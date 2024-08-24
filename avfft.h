/** AVFFT v1.71 C++ wrapper class written by Dmitry Boldyrev
 **
 **  File: avfft.h
 **  Main AVFFT class definition. You only need to include this file into your project
 **  for AVFFT function, and you also need to add the .S file into your XCode project 
 **  to be compiled-in.
 **
 **  This FFT wrapper is based on AV MPEG AUDIO source code, so I am not really
 **  going to officially claim any (C) to this, because honestly I just put it together
 **  into an easy to use C++ templated wrapper that supports SIMD vectorization via template.
 **
 **  I would appreciate credits in your app if you decide to use this convinient and well tuned
 **  FFT class, and send me email link to your app for me to appreciate what you have created.
 **
 **  GITHUB: https://github.com/mewza
 **  Email: subband@protonmail.com
 **/

#pragma once

#include "const1.h"

// USE_NEON must be defined as true or false for NEON asm optimizations
// Note: for now it only works if zfloat is defined as float not double

#define USE_NEON false

typedef cmplxT<zfloat> FFTComplex;

#define DEF_ALIGNED(n,t,v) extern t __attribute__ ((aligned (n))) v
#define DEF_COSINTABLE(size) DEF_ALIGNED(32, zfloat, ff_cos_##size)[size/2]; \
                         DEF_ALIGNED(32, zfloat, ff_sin_##size)[size/2]

#ifdef __cplusplus
extern "C" {
#endif

// needed to be "C" for asm links

void ff_fft_permute_neon(void *s, FFTComplex *z);
void ff_fft_calc_neon(void *s, FFTComplex *z);

#ifdef __cplusplus
}
#endif

void ff_init_ff_cosin_tabs(int index);
zfloat *get_ff_cos_tab(int index);
zfloat *get_ff_sin_tab(int index);
void pass_neon(FFTComplex*, zfloat *, unsigned int);

/* cos(2*pi*x/n) for 0<=x<=n/4, followed by its reverse */

DEF_COSINTABLE(16);
DEF_COSINTABLE(32);
DEF_COSINTABLE(64);
DEF_COSINTABLE(128);
DEF_COSINTABLE(256);
DEF_COSINTABLE(512);
DEF_COSINTABLE(1024);
DEF_COSINTABLE(2048);
DEF_COSINTABLE(4096);
DEF_COSINTABLE(8192);
DEF_COSINTABLE(16384);
DEF_COSINTABLE(32768);
DEF_COSINTABLE(65536);
DEF_COSINTABLE(131072);

static inline constexpr unsigned int floorlog2(unsigned int x) {
    return (x == 1) ? 0 : 1 + floorlog2(x >> 1);
}

#define M_TWOPI (M_PI * 2.)
#define NEON_ASM_COND (USE_NEON && std::is_same_v<T, float>)

template <typename T>
class AVFFT
{
    using T1 = SimdBase<T>;
    typedef cmplxT<T> ComplexT;
    
public:
    
    void real_fft(const T* x, cmplxT<T>* y, int N, bool do_scale = false)
    {
        rfft(x, (T*)y, N, true, do_scale);
    }

    void real_ifft(const cmplxT<T>* x, T* y, int N, bool do_scale = false)
    {
        rfft((T*)x, y, N, false, do_scale);
    }

    void rfft(const T* x, T* y, int N, bool forward, bool do_scale = false)
    {
        ComplexT *yp = (ComplexT*)y;
        const int N2 = N/2, N4 = N2/2, N8 = N4/2;
        zfloat *ctab = get_ff_cos_tab(N);
        zfloat *stab = get_ff_sin_tab(N);
        
        cmplxT<T> *p, *q, tw, sum, diff;
        T tw1, tw2;
        
        if (x != y) {
            memcpy(y, x, (forward)? (N * sizeof(T)) : ((N2+1) * sizeof(cmplxT<T>)));
        }
        if (forward) {
            cmplx_fft(yp, yp, N2, true);
            T1 scv = 0.5;
            cmplxT<T> *yb = yp, *ye = &yp[N];
            if (do_scale) scv *= 1./(T1)N2;
            while (yb != ye) *yb++ *= scv;
        }
        
        // Shamelessly borrowed from the Source: http://www.katjaas.nl/realFFT/realFFT2.html
        
        for (int i = 1; i < N4; ++i)
        {
            p = &yp[i];
            q = &yp[N2-i];
            
            if (i < N8) {
                tw = ComplexT(ctab[i],stab[i]);
            } else if (i > N8) {
                tw = ComplexT(ctab[N4-i],stab[N4-i]);
            } else {
                tw = ComplexT(M_SQRT1_2, M_SQRT1_2);
            }
            if (forward) tw.re = -tw.re;
            
            sum  = *p + *q;
            diff = *p - *q;
            
            tw1 = tw.re * sum.im + tw.im * diff.re;
            tw2 = tw.im * sum.im - tw.re * diff.re;
            
            p->re = sum.re - tw1;
            p->im = diff.im - tw2;
            q->re = sum.re + tw1;
            q->im = -(diff.im + tw2);
        }
        
        p = &yp[N4];
        p->re *=  2;
        p->im *= -2;
        
        if (!forward)  {
            cmplx_fft(yp, yp, N2, false);
            T1 scv = 0.5;
            T *ye = &y[N];
            if (do_scale) scv *= 1./(T1)N;
            while (y != ye) *y++ *= scv;
        }
    }

    void cmplx_fft(const cmplxT<T> *in, cmplxT<T> *out, int NC, bool forward, bool do_scale = false)
    {
        if (forward)
        {
            if (!fwd.initialized()) {
                fwd.init(NC, false);
            }
            
            memcpy(out, in, sizeof(T) * (NC));
            
            if constexpr(NEON_ASM_COND)  {
                ff_fft_permute_neon(&fwd, (FFTComplex*)out);
                ff_fft_calc_neon(&fwd, (FFTComplex*)out);
            } else {
                fwd.permute(out);
                fwd.calc(out, true);
            }
        } else
        {
            if (!rev.initialized()) {
                rev.init(NC, true);
            }
            
            memcpy(out, in, sizeof(cmplxT<T>) * (NC+1));
            
            if constexpr(NEON_ASM_COND)  {
                ff_fft_permute_neon(&rev, (FFTComplex*)out);
                ff_fft_calc_neon(&rev, (FFTComplex*)out);
            } else {
                rev.permute(out);
                rev.calc(out, false);
            }
        }
    }
    
protected:
    
    class FFTContext
    {
        enum fft_permutation_type
        {
            FF_FFT_PERM_DEFAULT,
            FF_FFT_PERM_SWAP_LSBS,
            FF_FFT_PERM_AVX,
        };

    public:
        
        FFTContext()
        {
            _revtab = NULL;
            _inverse = false;
            _initialized = false;
            _tmpbuf = NULL;
            _nbits = 0;
        }
        
        ~FFTContext()
        {
            if (_revtab) callocT_free(_revtab);
            _revtab = NULL;
            if (_tmpbuf) callocT_free(_tmpbuf);
            _tmpbuf = NULL;
        }
        
        uint16_t *get_perm_table() const { return _revtab; }
        
        bool initialized() const { return _initialized; }
        
        int split_radix_permutation(int i, int n, int inv)
        {
            if (n <= 2) return i & 1;
            int m = n >> 1;
            if (!(i&m)) return split_radix_permutation(i, m, inv) * 2;
            m >>= 1;
            if (inv == !(i&m))
                return split_radix_permutation(i, m, inv) * 4 + 1;
            else
                return split_radix_permutation(i, m, inv) * 4 - 1;
        }

        int is_second_half_of_fft32(int i, int n)
        {
            if (n <= 32)
                return i >= 16;
            else if (i < n/2)
                return is_second_half_of_fft32(i, n/2);
            else if (i < 3*n/4)
                return is_second_half_of_fft32(i - n/2, n/4);
            else
                return is_second_half_of_fft32(i - 3*n/4, n/4);
        }

        void perm_avx()
        {
            int n = 1 << _nbits;
            for (int i=0; i < n; i+=16)
            {
                if (is_second_half_of_fft32(i, n))
                {
                    for (int k=0; k < 16; k++)
                        _revtab[-split_radix_permutation(i + k, n, _inverse) & (n - 1)] = i + avx_tab[k];
                } else
                {
                    for (int k=0; k < 16; k++)
                    {
                        int j = i + k;
                        j = (j & ~7) | ((j >> 1) & 3) | ((j << 2) & 4);
                        _revtab[-split_radix_permutation(i + k, n, _inverse) & (n - 1)] = j;
                    }
                }
            }
        }
        
        int init(int N, bool inverse)
        {
            int i, j, n;

            _mdct_bits = floorlog2(N);
            _mdct_size = N;
        
            _revtab = NULL;

            if (_mdct_bits < 2 || _mdct_bits > 17)
                goto fail;
            _nbits = _mdct_bits;
            n = 1 << _nbits;
            _inverse = inverse;
           
            _tmpbuf = callocT<ComplexT>(n);
            if (!_tmpbuf)
                goto fail;
            
            _revtab = callocT<int>(n);
            if (!_revtab)
                goto fail;
           
            fft_permutation = FF_FFT_PERM_DEFAULT; //FF_FFT_PERM_AVX; //FF_FFT_PERM_SWAP_LSBS; //FF_FFT_PERM_DEFAULT;

            for(j = 4; j <= _nbits + 1; j++)
                ff_init_ff_cosin_tabs(j);

            if (fft_permutation == FF_FFT_PERM_AVX) {
                perm_avx();
            } else {
#define PROCESS_FFT_PERM_SWAP_LSBS() do { \
            for(i = 0; i < n; i++) {\
                int k;\
                j = i;\
                j = (j & ~3) | ((j >> 1) & 1) | ((j << 1) & 2);\
                k = -split_radix_permutation(i, n, _inverse) & (n - 1);\
                _revtab[k] = j;\
            } \
            } while (0);
                
#define PROCESS_FFT_PERM_DEFAULT(num) do { \
            for(i = 0; i < n; i++) {\
            int k;\
            j = i;\
            k = -split_radix_permutation(i, n, _inverse) & (n - 1);\
            _revtab[k] = j;\
            } \
            } while (0);
                
#define SPLIT_RADIX_PERMUTATION(num) do { \
            if (fft_permutation == FF_FFT_PERM_SWAP_LSBS) {\
                PROCESS_FFT_PERM_SWAP_LSBS(num) \
            } else {\
                PROCESS_FFT_PERM_DEFAULT(num) \
            } \
            } while (0);
                 
            if (_revtab) SPLIT_RADIX_PERMUTATION()
                        
#undef PROCESS_FFT_PERM_DEFAULT
#undef PROCESS_FFT_PERM_SWAP_LSBS
#undef SPLIT_RADIX_PERMUTATION
                        }
            _initialized = true;
            return 0;
         fail:
            if (_revtab)
                callocT_free(_revtab);
            _revtab = NULL;
            if (_tmpbuf)
                callocT_free(_tmpbuf);
            _tmpbuf = NULL;
            return -1;
        }

        void permute(ComplexT *z)
        {
            int np = 1 << _nbits;
            const int *rt = _revtab;
            ComplexT *tmp = (ComplexT*)_tmpbuf;
            
            // TODO: handle split-radix permute in a more optimal way, probably in-place
            if (rt) {
                for(int j=0;j<np;j++) tmp[rt[j]] = z[j];
            } else {
                fprintf(stderr, "WARNING: revtab nor revtab32 are allocated. Abort.\n");
                return;
            }
            for(int j=0;j<np;j++) z[j] = tmp[j];
        }
        
        void calc(ComplexT *z, bool forward)
        {
            if (initialized()) {
                if (forward)
                    fft_dispatch[_nbits-2](z);
                else
                    fft_dispatch[_nbits-2](z);
            } else
                fprintf(stderr, "WARNING: init() was not called before calc() as required. Abort.\n");
        }
        
        void calc2(ComplexT *z)
        {
            if (initialized())
                fft_dispatch[_nbits-2](z);
            else
                fprintf(stderr, "WARNING: init() was not called before calc() as required. Abort.\n");
        }
        
        
    protected:
        
        int         _nbits;         // +0
        int         _inverse;       // +4  forward / inverse transform
        int        *_revtab;       // +8
        void       *_tmpbuf;       // +16
        int         _mdct_size;     // +24 size of MDCT (i.e. number of input data * 2)
        int         _mdct_bits;     // +28 n = 2^nbits
        bool        _initialized;

        enum fft_permutation_type fft_permutation;
        
        static constexpr int avx_tab[] = {
            0, 4, 1, 5, 8, 12, 9, 13, 2, 6, 3, 7, 10, 14, 11, 15
        };

#define BF(x, y, a, b) \
        x = a - b; \
        y = a + b;

#define BUTTERFLIES(a0,a1,a2,a3) { \
        BF(t3, t5, t5, t1); \
        BF(a2.re, a0.re, a0.re, t5); \
        BF(a3.im, a1.im, a1.im, t3); \
        BF(t4, t6, t2, t6); \
        BF(a3.re, a1.re, a1.re, t4); \
        BF(a2.im, a0.im, a0.im, t6); \
        }

    // force loading all the inputs before storing any.
    // this is slightly slower for small data, but avoids
    // store->load aliasing for addresses separated by large
    // powers of 2

#define BUTTERFLIES_BIG(a0,a1,a2,a3) { \
        T r0=a0.re, i0=a0.im, r1=a1.re, i1=a1.im; \
        BF(t3, t5, t5, t1); \
        BF(a2.re, a0.re, r0, t5); \
        BF(a3.im, a1.im, i1, t3); \
        BF(t4, t6, t2, t6); \
        BF(a3.re, a1.re, r1, t4); \
        BF(a2.im, a0.im, i0, t6); \
        }

#define sqrthalf   M_SQRT1_2

#define CMUL(dre, dim, are, aim, bre, bim) \
        (dre) = (are) * (bre) - (aim) * (bim); \
        (dim) = (are) * (bim) + (aim) * (bre);

#define TRANSFORM(a0,a1,a2,a3,wre,wim) { \
        CMUL(t1, t2, a2.re, a2.im, wre, -wim); \
        CMUL(t5, t6, a3.re, a3.im, wre,  wim); \
        BUTTERFLIES(a0,a1,a2,a3) \
        }

#define TRANSFORM_ZERO(a0,a1,a2,a3) { \
        t1 = a2.re; \
        t2 = a2.im; \
        t5 = a3.re; \
        t6 = a3.im; \
        BUTTERFLIES(a0,a1,a2,a3) \
        }

    // z[0...8n-1], w[1...2n-1]

#define PASS(name) \
        static void name(ComplexT *z, const zfloat *wre, unsigned int n)\
        {\
            T t1, t2, t3, t4, t5, t6; \
            int o1 = 2*n; \
            int o2 = 4*n; \
            int o3 = 6*n; \
            const zfloat *wim = wre+o1; \
            n--; \
            \
            TRANSFORM_ZERO(z[0],z[o1],z[o2],z[o3]); \
            TRANSFORM(z[1],z[o1+1],z[o2+1],z[o3+1],wre[1],wim[-1]); \
            do { \
                z += 2; \
                wre += 2; \
                wim -= 2; \
                TRANSFORM(z[0],z[o1],z[o2],z[o3],wre[0],wim[0]); \
                TRANSFORM(z[1],z[o1+1],z[o2+1],z[o3+1],wre[1],wim[-1]); \
            } while(--n); \
        }
        PASS(pass)

#undef BUTTERFLIES
        
#define BUTTERFLIES BUTTERFLIES_BIG
        PASS(pass_big)

#define DECL_FFT(n,n2,n4) \
        static void fft##n(ComplexT *z) \
        { \
            fft##n2(z); \
            fft##n4(z+n4*2); \
            fft##n4(z+n4*3); \
            pass(z,ff_cos_##n,n4/2); \
        }
      
        static inline void fft4(ComplexT *z)
        {
            T t1, t2, t3, t4, t5, t6, t7, t8;
            
            BF(t3, t1, z[0].re, z[1].re);
            BF(t8, t6, z[3].re, z[2].re);
            BF(z[2].re, z[0].re, t1, t6);
            BF(t4, t2, z[0].im, z[1].im);
            BF(t7, t5, z[2].im, z[3].im);
            BF(z[3].im, z[1].im, t4, t8);
            BF(z[3].re, z[1].re, t3, t7);
            BF(z[2].im, z[0].im, t2, t5);
        }

        static inline void fft8(ComplexT *z)
        {
           T t1, t2, t3, t4, t5, t6;
                
            fft4(z);
            
            BF(t1, z[5].re, z[4].re, -z[5].re);
            BF(t2, z[5].im, z[4].im, -z[5].im);
            BF(t5, z[7].re, z[6].re, -z[7].re);
            BF(t6, z[7].im, z[6].im, -z[7].im);
            
            BUTTERFLIES(z[0],z[2],z[4],z[6]);
            TRANSFORM(z[1],z[3],z[5],z[7],sqrthalf,sqrthalf);
        }

        static inline void fft16(ComplexT *z)
        {
            T t1, t2, t3, t4, t5, t6;
            zfloat cos_16_1 = ff_cos_16[1];
            zfloat cos_16_3 = ff_cos_16[3];
            
            fft8(z);
            fft4(z+8);
            fft4(z+12);
            
            TRANSFORM_ZERO(z[0],z[4],z[8],z[12]);
            TRANSFORM(z[2],z[6],z[10],z[14],sqrthalf,sqrthalf);
            TRANSFORM(z[1],z[5],z[9],z[13],cos_16_1,cos_16_3);
            TRANSFORM(z[3],z[7],z[11],z[15],cos_16_3,cos_16_1);
        }
        
        DECL_FFT(32,16,8)
        DECL_FFT(64,32,16)
        DECL_FFT(128,64,32)
        DECL_FFT(256,128,64)
        DECL_FFT(512,256,128)
        
#define pass pass_big
        
        DECL_FFT(1024,512,256)
        DECL_FFT(2048,1024,512)
        DECL_FFT(4096,2048,1024)
        DECL_FFT(8192,4096,2048)
        DECL_FFT(16384,8192,4096)
        DECL_FFT(32768,16384,8192)
        DECL_FFT(65536,32768,16384)
        DECL_FFT(131072,65536,32768)

        static constexpr void (* const fft_dispatch[])(ComplexT*) = {
            fft4, fft8, fft16, fft32, fft64, fft128, fft256, fft512, fft1024,
            fft2048, fft4096, fft8192, fft16384, fft32768, fft65536, fft131072
        };
    };
    
    FFTContext  fwd, rev;
};
