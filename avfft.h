/** AVFFT v1.6 C++ wrapper class written by Dmitry Boldyrev
 **
 **  GITHUB: https://github.com/mewza
 **  Email: subband@protonmail.com
 **
 **  This FFT wrapper is based on AV MPEG AUDIO source code, so I am not really
 **  going to officially claim any (C) to it, because I honestly,
 **  just put it together into an aeasy to use C++ templated wrapper,
 **  the great benefit I tested it, it works well with 32-bit and 64-bit
 **  floating point single types or as a form of intrinsic SIMD vectors.
 **
 **  This is a much more robust version of AVFFT class. I made it 100% compatible
 **  with PFFFT, so you can just drop replace it and it will produce almost identical output,
 **  plus the intrinsic vector support via template T parameter.
 **
 **  I would appreciate credits in the app if you use this fancy c++ wrapper, and I might try to
 **  add assembler optimizations next.
 **/

#pragma once

#include "const1.h"

// USE_NEON must be defined as true or false for NEON asm optimizations
// Note: for now it only works if zfloat is defined as float not double

#define USE_NEON true

#ifdef __cplusplus
extern "C" {
#endif

typedef cmplxT<zfloat> FFTComplex;

typedef struct CosTabsInitOnce {
    void (*func)(void);
    bool control;
} CosTabsInitOnce;

#ifdef __cplusplus
}
#endif
 
template <typename TT> class FFTContext;

#define DECLARE_ALIGNED(n,t,v) t __attribute__ ((aligned (n))) v
#define COSTABLE(size) DECLARE_ALIGNED(32, zfloat, ff_cos_##size)[size/2]

#ifdef __cplusplus
extern "C" {
#endif

void ff_fft_permute_neon(void *s, FFTComplex *z);
void ff_fft_calc_neon(void *s, FFTComplex *z);

static inline void ff_init_ff_cos_tabs(int index);
static inline void pass_neon(FFTComplex*, zfloat *, unsigned int);

/* cos(2*pi*x/n) for 0<=x<=n/4, followed by its reverse */

COSTABLE(16);
COSTABLE(32);
COSTABLE(64);
COSTABLE(128);
COSTABLE(256);
COSTABLE(512);
COSTABLE(1024);
COSTABLE(2048);
COSTABLE(4096);
COSTABLE(8192);
COSTABLE(16384);
COSTABLE(32768);
COSTABLE(65536);
COSTABLE(131072);

#ifdef __cplusplus
}
#endif

#define M_TWOPI (M_PI * 2.)
#define NEON_ASM_COND (USE_NEON && std::is_same_v<T, float>)

template <typename T>
class AVFFT
{
    using FT =  std::conditional_t<std::is_same_v<T, float8v>, float,
                std::conditional_t<std::is_same_v<T, float4v>, float,
                std::conditional_t<std::is_same_v<T, float2v>, float,
                std::conditional_t<std::is_same_v<T, float>, float,
                std::conditional_t<std::is_same_v<T, double8v>, double,
                std::conditional_t<std::is_same_v<T, double4v>, double,
                std::conditional_t<std::is_same_v<T, double2v>, double,
                std::conditional_t<std::is_same_v<T, double>, double, T >>>>>>>>;

    typedef cmplxT<T> ComplexT;

public:
    
    AVFFT() {}

    void real_fft(const T* x, T* y, int N, bool forward, bool do_scale = false)
    {
        FT wpr, wpi, theta;
        T xr, xi, c1, c2, h1r, h1i, h2r, h2i, wr, wi, temp;
        int N2 = N >> 1, N4 = N >> 2;
        int i, i1, i2, i3, i4, N2p1;

        if (forward) {
            memcpy(y, x, sizeof(T)*N);
            y[N] = y[N+1] = 0.0;
        } else
            memcpy(y, x, sizeof(T)*(N+2));
        
        theta = M_PI / (FT)(N2);
        wr = 1.;
        wi = 0.;
        c1 = 0.5;

        if (forward) {
            c2 = -0.5;
            cmplx_fft((cmplxT<T>*)y, (cmplxT<T>*)y, N2, forward, do_scale);
            xr = y[0];
            xi = y[1];
        }
        else {
            c2 = 0.5;
            theta = -theta;
            xr = y[1];
            xi = 0.;
            y[1] = 0.;
        }

        wpr = (-2.*F_POW(F_SIN(0.5 * theta), 2.));
        wpi = F_SIN(theta);
        N2p1 = N + 1;

        for (i = 0; i <= N4; i++) {
            i1 = i << 1;
            i2 = i1 + 1;
            i3 = N2p1 - i2;
            i4 = i3 + 1;
            if (i == 0) 
            {
                h1r =  c1*(y[i1] + xr);
                h1i =  c1*(y[i2] - xi);
                h2r = -c2*(y[i2] + xi);
                h2i =  c2*(y[i1] - xr);
                y[i1] = h1r + wr*h2r - wi*h2i;
                y[i2] = h1i + wr*h2i + wi*h2r;
                xr =  h1r - wr*h2r - wi*h2i;
                xi = -h1i + wr*h2i + wi*h2r;
            } else
            {
                h1r = c1*(y[i1] + y[i3]);
                h1i = c1*(y[i2] - y[i4]);
                h2r = -c2*(y[i2] + y[i4]);
                h2i = c2*(y[i1] - y[i3]);
                y[i1] = h1r + wr*h2r - wi*h2i;
                y[i2] = h1i + wr*h2i + wi*h2r;
                y[i3] = h1r - wr*h2r + wi*h2i;
                y[i4] = -h1i + wr*h2i + wi*h2r;
            }

            wr = (temp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + temp*wpi + wi;
        }

        if (forward) {
          //  y[1] = xr;
            y[1] = y[N];
            y[N] = 0.0;
        } else {
            cmplx_fft((cmplxT<T>*)y, (cmplxT<T>*)y, N2, forward, do_scale);
        }
    }

    void cmplx_fft(const cmplxT<T> *in, cmplxT<T> *out, int NC, bool forward, bool scale = false)
    {
        int ND = NC << 1;
        if (forward) {
            if (!fwd.initialized()) {
                fwd.init(NC, false);
            }
            memcpy(out, in, sizeof(cmplxT<T>) * (NC+1));
           
            if constexpr(NEON_ASM_COND)  {
                ff_fft_permute_neon(&fwd, (FFTComplex*)out);
                ff_fft_calc_neon(&fwd, (FFTComplex*)out);
            } else {
                fwd.permute(out);
                fwd.calc(out);
            }
            zfloat scv = 1.0/(zfloat)ND;
            cmplxT<T> *xi = out, *xe = &out[NC];
            if (scale) {
                while (xi < xe) *xi++ *= scv;
            }
        } else {
            if (!rev.initialized()) {
                rev.init(NC, true);
            }
            memcpy(out, in, sizeof(cmplxT<T>) * (NC+1));
          
            if constexpr(NEON_ASM_COND)  {
                ff_fft_permute_neon(&rev, (FFTComplex*)out);
                ff_fft_calc_neon(&rev, (FFTComplex*)out);
            } else {
                rev.permute(out);
                rev.calc(out);
            }
            zfloat scv = 2.0;
            cmplxT<T> *xi = out, *xe = &out[NC];
            if (scale) {
                while (xi < xe) *xi++ *= scv;
            }
        }
    }
   
protected:
    
    template <typename TT>
    class FFTContext
    {
        typedef cmplxT<TT> ComplexT;
        
        enum fft_permutation_type {
            FF_FFT_PERM_DEFAULT,
            FF_FFT_PERM_SWAP_LSBS,
            FF_FFT_PERM_AVX,
        };

    public:
        FFTContext()
        {
            _revtab = NULL;
            _revtab32 = NULL;
            _inverse = false;
            _initialized = false;
            _tmpbuf = NULL;
            _nbits = 0;
        }
        
        ~FFTContext()
        {
            if (_revtab) callocT_free(_revtab);
            _revtab = NULL;
            if (_revtab32) callocT_free(_revtab32);
            _revtab32 = NULL;
            if (_tmpbuf) callocT_free(_tmpbuf);
            _tmpbuf = NULL;
        }
        
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
                    for (int k=0; k < 16; k++) {
                        _revtab[-split_radix_permutation(i + k, n, _inverse) & (n - 1)] = i + avx_tab[k];
                    }
                } else {
                    for (int k=0; k < 16; k++) {
                        int j = i + k;
                        j = (j & ~7) | ((j >> 1) & 3) | ((j << 2) & 4);
                        _revtab[-split_radix_permutation(i + k, n, _inverse) & (n - 1)] = j;
                    }
                }
            }
        }
        
        static inline constexpr unsigned int floorlog2(unsigned int x) {
            return (x == 1) ? 0 : 1 + floorlog2(x >> 1);
        }

        int init(int N, bool inverse)
        {
            int i, j, n;

            _mdct_bits = floorlog2(N);
            _mdct_size = N;
            
            _revtab = NULL;
            _revtab32 = NULL;

            if (_mdct_bits < 2 || _mdct_bits > 17)
                goto fail;
            _nbits = _mdct_bits;
            n = 1 << _nbits;
            _inverse = inverse;
           
            _tmpbuf = callocT<ComplexT>(n);
            if (!_tmpbuf)
                goto fail;
            
            if (_nbits <= 16) {
                _revtab = callocT<uint16_t>(n);
                if (!_revtab)
                    goto fail;
            } else {
                _revtab32 = callocT<uint32_t>(n);
                if (!_revtab32)
                    goto fail;
            }
            fft_permutation = FF_FFT_PERM_DEFAULT; //FF_FFT_PERM_AVX; //FF_FFT_PERM_SWAP_LSBS; //FF_FFT_PERM_DEFAULT;

            for(j=4; j <= _nbits; j++) {
                ff_init_ff_cos_tabs(j);
            }

            if (fft_permutation == FF_FFT_PERM_AVX) {
                perm_avx();
            } else {
#define PROCESS_FFT_PERM_SWAP_LSBS(num) do { \
            for(i = 0; i < n; i++) {\
                int k;\
                j = i;\
                j = (j & ~3) | ((j >> 1) & 1) | ((j << 1) & 2);\
                k = -split_radix_permutation(i, n, _inverse) & (n - 1);\
                _revtab##num[k] = j;\
            } \
            } while (0);
                
#define PROCESS_FFT_PERM_DEFAULT(num) do { \
            for(i = 0; i < n; i++) {\
            int k;\
            j = i;\
            k = -split_radix_permutation(i, n, _inverse) & (n - 1);\
            _revtab##num[k] = j;\
            } \
            } while (0);
                
#define SPLIT_RADIX_PERMUTATION(num) do { \
            if (fft_permutation == FF_FFT_PERM_SWAP_LSBS) {\
                PROCESS_FFT_PERM_SWAP_LSBS(num) \
            } else {\
                PROCESS_FFT_PERM_DEFAULT(num) \
            } \
            } while (0);
                 
            if (_revtab)
                SPLIT_RADIX_PERMUTATION()
                
            if (_revtab32)
                SPLIT_RADIX_PERMUTATION(32)
                        
#undef PROCESS_FFT_PERM_DEFAULT
#undef PROCESS_FFT_PERM_SWAP_LSBS
#undef SPLIT_RADIX_PERMUTATION
                        }
            _initialized = true;
            return 0;
         fail:
            if (_revtab) callocT_free(_revtab);
            _revtab = NULL;
            if (_revtab32) callocT_free(_revtab32);
            _revtab32 = NULL;
            if (_tmpbuf) callocT_free(_tmpbuf);
            _tmpbuf = NULL;
            return -1;
        }

        void permute(ComplexT *z)
        {
            int np = 1 << _nbits;
            const uint16_t *rt = _revtab;
            const uint32_t *rt32 = _revtab32;
            ComplexT *tmp = (ComplexT*)_tmpbuf;
            
            // TODO: handle split-radix permute in a more optimal way, probably in-place
            if (rt) {
                for(int j=0;j<np;j++) tmp[rt[j]] = z[j];
            } else if (rt32) {
                for(int j=0;j<np;j++) tmp[rt32[j]] = z[j];
            } else {
                fprintf(stderr, "WARNING: revtab nor revtab32 are allocated. Abort.\n");
                return;
            }
            for(int j=0;j<np;j++) z[j] = tmp[j];
        }
        
        void calc(ComplexT *z)
        {
            if (initialized())
                fft_dispatch[_nbits-2](z);
            else
                fprintf(stderr, "WARNING: init() was not called before calc() as required. Abort.\n");
        }

    protected:
        
        int         _nbits;         // +0
        int         _inverse;       // +4  forward / inverse transform
        uint16_t    *_revtab;       // +8
        void        *_tmpbuf;       // +16
        int         _mdct_size;     // +24 size of MDCT (i.e. number of input data * 2)
        int         _mdct_bits;     // +28 n = 2^nbits
        uint32_t    *_revtab32;
        bool        _initialized;

        enum fft_permutation_type fft_permutation;
        
        static constexpr int avx_tab[] = {
            0, 4, 1, 5, 8, 12, 9, 13, 2, 6, 3, 7, 10, 14, 11, 15
        };

#define BF(x, y, a, b) \
        x = a - b; \
        y = a + b;

#define BUTTERFLIES(a0,a1,a2,a3) {\
        BF(t3, t5, t5, t1);\
        BF(a2.re, a0.re, a0.re, t5);\
        BF(a3.im, a1.im, a1.im, t3);\
        BF(t4, t6, t2, t6);\
        BF(a3.re, a1.re, a1.re, t4);\
        BF(a2.im, a0.im, a0.im, t6);\
        }

    // force loading all the inputs before storing any.
    // this is slightly slower for small data, but avoids
    // store->load aliasing for addresses separated by large
    // powers of 2

#define BUTTERFLIES_BIG(a0,a1,a2,a3) {\
        TT r0=a0.re, i0=a0.im, r1=a1.re, i1=a1.im;\
        BF(t3, t5, t5, t1);\
        BF(a2.re, a0.re, r0, t5);\
        BF(a3.im, a1.im, i1, t3);\
        BF(t4, t6, t2, t6);\
        BF(a3.re, a1.re, r1, t4);\
        BF(a2.im, a0.im, i0, t6);\
        }

#define sqrthalf   M_SQRT1_2

#define CMUL(dre, dim, are, aim, bre, bim) \
        (dre) = (are) * (bre) - (aim) * (bim);  \
        (dim) = (are) * (bim) + (aim) * (bre);

#define TRANSFORM(a0,a1,a2,a3,wre,wim) {\
        CMUL(t1, t2, a2.re, a2.im, wre, -wim);\
        CMUL(t5, t6, a3.re, a3.im, wre,  wim);\
        BUTTERFLIES(a0,a1,a2,a3)\
        }

#define TRANSFORM_ZERO(a0,a1,a2,a3) {\
        t1 = a2.re;\
        t2 = a2.im;\
        t5 = a3.re;\
        t6 = a3.im;\
        BUTTERFLIES(a0,a1,a2,a3)\
        }

    // z[0...8n-1], w[1...2n-1]

#define PASS(name) \
        static void name(ComplexT *z, const zfloat *wre, unsigned int n)\
        {\
            TT t1, t2, t3, t4, t5, t6;\
            int o1 = 2*n;\
            int o2 = 4*n;\
            int o3 = 6*n;\
            const zfloat *wim = wre+o1;\
            n--;\
            \
            TRANSFORM_ZERO(z[0],z[o1],z[o2],z[o3]);\
            TRANSFORM(z[1],z[o1+1],z[o2+1],z[o3+1],wre[1],wim[-1]);\
            do {\
                z += 2;\
                wre += 2;\
                wim -= 2;\
                TRANSFORM(z[0],z[o1],z[o2],z[o3],wre[0],wim[0]);\
                TRANSFORM(z[1],z[o1+1],z[o2+1],z[o3+1],wre[1],wim[-1]);\
            } while(--n);\
        }
        PASS(pass)

#undef BUTTERFLIES
        
#define BUTTERFLIES BUTTERFLIES_BIG
        PASS(pass_big)

#define DECL_FFT(n,n2,n4)\
        static void fft##n(ComplexT *z)\
        {\
            fft##n2(z);\
            fft##n4(z+n4*2);\
            fft##n4(z+n4*3);\
            if constexpr(NEON_ASM_COND) \
                pass_neon((FFTComplex*)z,ff_cos_##n,n4/2);\
            else \
               pass(z,ff_cos_##n,n4/2);\
        }
      
        static inline void fft4(ComplexT *z)
        {
            if constexpr(NEON_ASM_COND)
                fft4_neon(z);
                else {
                    TT t1, t2, t3, t4, t5, t6, t7, t8;
                    
                    BF(t3, t1, z[0].re, z[1].re);
                    BF(t8, t6, z[3].re, z[2].re);
                    BF(z[2].re, z[0].re, t1, t6);
                    BF(t4, t2, z[0].im, z[1].im);
                    BF(t7, t5, z[2].im, z[3].im);
                    BF(z[3].im, z[1].im, t4, t8);
                    BF(z[3].re, z[1].re, t3, t7);
                    BF(z[2].im, z[0].im, t2, t5);
                }
        }

        static inline void fft8(ComplexT *z)
        {
            if constexpr(NEON_ASM_COND)
                fft8_neon(z);
            else {
                TT t1, t2, t3, t4, t5, t6;
                
                fft4(z);
                
                BF(t1, z[5].re, z[4].re, -z[5].re);
                BF(t2, z[5].im, z[4].im, -z[5].im);
                BF(t5, z[7].re, z[6].re, -z[7].re);
                BF(t6, z[7].im, z[6].im, -z[7].im);
                
                BUTTERFLIES(z[0],z[2],z[4],z[6]);
                TRANSFORM(z[1],z[3],z[5],z[7],sqrthalf,sqrthalf);
            }
        }

        static inline void fft16(ComplexT *z)
        {
            if constexpr(NEON_ASM_COND)
                fft16_neon(z);
            else {
                TT t1, t2, t3, t4, t5, t6;
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
    } __attribute__ ((aligned (32)));
                     
    FFTContext<T>  fwd, rev;
  
};

#ifdef __cplusplus
extern "C" {
#endif

static constexpr zfloat * ff_cos_tabs[18] = {
    NULL, NULL, NULL, NULL,
    ff_cos_16,
    ff_cos_32,
    ff_cos_64,
    ff_cos_128,
    ff_cos_256,
    ff_cos_512,
    ff_cos_1024,
    ff_cos_2048,
    ff_cos_4096,
    ff_cos_8192,
    ff_cos_16384,
    ff_cos_32768,
    ff_cos_65536,
    ff_cos_131072,
};
static void init_ff_cos_tabs(int index)
{
    int i;
    int m = 1<<index;
    zfloat freq = 2.*M_PI/(zfloat)m;
    zfloat *tab = ff_cos_tabs[index];
    for(i=0; i<=m/4; i++)
        tab[i] = cos(i*freq);
    for(i=1; i<m/4; i++)
        tab[m/2-i] = tab[i];
}

#define INIT_FF_COS_TABS_FUNC(index, size)   \
static void init_ff_cos_tabs_ ## size (void) \
{                                            \
init_ff_cos_tabs(index);                 \
}

INIT_FF_COS_TABS_FUNC(4, 16)
INIT_FF_COS_TABS_FUNC(5, 32)
INIT_FF_COS_TABS_FUNC(6, 64)
INIT_FF_COS_TABS_FUNC(7, 128)
INIT_FF_COS_TABS_FUNC(8, 256)
INIT_FF_COS_TABS_FUNC(9, 512)
INIT_FF_COS_TABS_FUNC(10, 1024)
INIT_FF_COS_TABS_FUNC(11, 2048)
INIT_FF_COS_TABS_FUNC(12, 4096)
INIT_FF_COS_TABS_FUNC(13, 8192)
INIT_FF_COS_TABS_FUNC(14, 16384)
INIT_FF_COS_TABS_FUNC(15, 32768)
INIT_FF_COS_TABS_FUNC(16, 65536)
INIT_FF_COS_TABS_FUNC(17, 131072)

static CosTabsInitOnce cos_tabs_init_once[18] = {
    { NULL },
    { NULL },
    { NULL },
    { NULL },
    { init_ff_cos_tabs_16, false },
    { init_ff_cos_tabs_32, false },
    { init_ff_cos_tabs_64, false },
    { init_ff_cos_tabs_128, false },
    { init_ff_cos_tabs_256, false },
    { init_ff_cos_tabs_512, false },
    { init_ff_cos_tabs_1024, false },
    { init_ff_cos_tabs_2048, false },
    { init_ff_cos_tabs_4096, false },
    { init_ff_cos_tabs_8192, false },
    { init_ff_cos_tabs_16384, false },
    { init_ff_cos_tabs_32768, false },
    { init_ff_cos_tabs_65536, false },
    { init_ff_cos_tabs_131072, false },
};

static inline void ff_init_ff_cos_tabs(int index)
{
    if (!cos_tabs_init_once[index].control) {
        cos_tabs_init_once[index].func();
        cos_tabs_init_once[index].control = true;
    }
}

#ifdef __cplusplus
}
#endif
