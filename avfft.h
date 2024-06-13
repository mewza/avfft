#pragma once

/** AVFFT v1.3 C++ wrapper class written by Dmitry Boldyrev
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

// you can configure FFT to run as float or double by changing zfloat
// preicison floating point (default = double)

// for intrinsic vector support on iOS SDK uncomment
//#define SUPPORT_IVECTOR

#ifdef SUPPORT_IVECTOR
#include <simd/simd.h>

typedef simd_double2 zfloat2;
typedef simd_double4 zfloat4;
typedef simd_double8 zfloat8;

typedef simd_float2 float2v;
typedef simd_float4 float4v;
typedef simd_float8 float8v;

typedef simd_double2 double2v;
typedef simd_double4 double4v;
typedef simd_double8 double8v;

static __inline float8v F_SQRT(float8v s) { return simd::sqrt(s); }
static __inline double8v F_SQRT(double8v s) { return simd::sqrt(s); }
static __inline float4v F_SQRT(float4v s) { return simd::sqrt(s); }
static __inline double4v F_SQRT(double4v s) { return simd::sqrt(s); }
static __inline float2v F_SQRT(float2v s) { return simd::sqrt(s); }
static __inline double2v F_SQRT(double2v s) { return simd::sqrt(s); }
static __inline double F_SQRT(double s) { return __builtin_sqrt(s); }
static __inline float F_SQRT(float s) { return __builtin_sqrt(s); }

#endif // SUPPORT_IVECTOR


typedef double zfloat;

#ifndef D_CALLOCT
#define D_CALLOCT

#define MALLOC_ALIGNMENT 64

static inline void *aligned_malloc( size_t nb_bytes )
{
    void *p, *p0 = malloc(nb_bytes + MALLOC_ALIGNMENT);
    if (!p0) return (void *) 0;
   
    p = (void *) (((size_t) p0 + MALLOC_ALIGNMENT) & (~((size_t) (MALLOC_ALIGNMENT-1))));
    *((void **) p - 1) = p0;
   
    return p;
}

static inline void aligned_free( void *p )
{
    if (p) free(*((void **) p - 1));
}

template <typename T>
static inline T *callocT( int64_t k )
{
    T *p = (T*)aligned_malloc (k * sizeof(T));
    if (!p)
        return NULL;
    memset (p, 0, k * sizeof (T));
    return p;
}

static inline void callocT_free( void *p )
{
    aligned_free(p);
}

#endif // D_CALLOCT

#ifndef D_CMPLXT
#define D_CMPLXT

template <typename T>
struct cmplxT {
public:
    cmplxT(T r, T i) { re = r; im = i; }
    cmplxT(const cmplxT& v) {
        re = v.re; im = v.im;
    }
    cmplxT(long double v) {
        re = v; im = v;
    }
    cmplxT() { re = 0.0; im = 0.0; }
    T mag() const {
        return F_SQRT(re * re + im * im);
    }
    
    inline cmplxT<T> operator * (const double d) {
        return cmplxT(re * d, im * d);
    }
    inline constexpr cmplxT<T> operator *= (long double d) const {
        return cmplxT(re * d, im * d);
    }
    inline constexpr cmplxT<T>& operator += (const cmplxT<T>& d) {
        return cmplxT(re + d.re, im + d.im);
    }
    T re;
    T im;
};
//__attribute__((packed));

#endif // D_CMPLXT

typedef struct CosTabsInitOnce {
    void (*func)(void);
    bool control;
} CosTabsInitOnce;

#define DECLARE_ALIGNED(n,t,v) t __attribute__ ((aligned (n))) v
#define COSTABLE(size) static DECLARE_ALIGNED(32, zfloat, ff_cos_##size)[size/2]

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

#define M_TWOPI (M_PI * 2.)

template <typename T>
class AVFFT
{
    static constexpr int floorlog2(int x) {
        return (x == 1) ? 0 : 1 + floorlog2(x >> 1);
    }
    typedef cmplxT<T> FFTComplex;
    
public:
    
    AVFFT() : fwd(this), rev(this) { }
    
    void real_fft(const T* in, T *out, int N, bool forward, bool scale = false)
    {
        alignas(64) cmplxT<T> fft[N];
        cmplxT<T> x, h1, h2;
        zfloat tmp, c2, theta;
        cmplxT<zfloat> wp, w;
        const zfloat c1 = 0.5;

        zfloat scv = 1.;
        if (scale) scv = 1./(zfloat)N;
        
        w = cmplxT<zfloat>(1,0);
        theta = M_PI/(zfloat)N;
        
        if (forward)
        {
            for (int i=0; i<N; i++) {
                fft[i] = cmplxT<T>(in[i] * scv, 0.0);
            }
            cmplx_fft(fft, fft, N, forward);
           
            c2 = -0.5;
            x = fft[0];
        } else 
        {
            memcpy(fft, in, sizeof(cmplxT<T>) * N);
            c2 = 0.5;
            theta = -theta;
            x = cmplxT<T>(fft[0].im, 0);
            fft[0].im = 0;
        }
        const zfloat xx = sin(0.5*theta);
        wp = cmplxT<zfloat>(-2. * xx * xx, sin(theta));
        
        for (int i=0; i <= N/2; i++) {
            int ir = N-i;
            if (i) {
                h1 = cmplxT<T>( fft[i].re + fft[ir].re, fft[i].im - fft[ir].im) * c1;
                h2 = cmplxT<T>(-fft[i].im - fft[ir].im, fft[i].re - fft[ir].re) * c2;
                fft[i] =  cmplxT<T>(h1.re + w.re * h2.re - w.im * h2.im,  h1.im + w.re * h2.im + w.im * h2.re);
                fft[ir] = cmplxT<T>(h1.re - w.re * h2.re + w.im * h2.im, -h1.im + w.re * h2.im + w.im * h2.re);
            } else {
                h1 = cmplxT<T>( fft[i].re + x.re, fft[i].im - x.im) * c1;
                h2 = cmplxT<T>(-fft[i].im - x.im, fft[i].re - x.re) * c2;
                fft[i] = cmplxT<T>(h1.re + w.re * h2.re - w.im * h2.im,  h1.im + w.re * h2.im + w.im * h2.re);
                x      = cmplxT<T>(h1.re - w.re * h2.re - w.im * h2.im, -h1.im + w.re * h2.im + w.im * h2.re);
            }
            w.re = (tmp = w.re) * wp.re - w.im * wp.im + w.re;
            w.im = w.im * wp.re + tmp * wp.im + w.im;
        }
        
        if (forward) {
            //fft[0].im = x.re;
            fft[0].im = fft[N/2].re;
            fft[N/2].re = 0.0;
            memcpy(out, fft, sizeof(cmplxT<T>)*N);
        } else
        {
            cmplx_fft(fft, fft, N, forward);
            for (int i=0; i<N; i++)
               out[i] = fft[i].re * scv;
        }
    }

    void cmplx_fft(const cmplxT<T> *in, cmplxT<T> *out, int N, bool forward)
    {
        if (forward) {
            if (!fwd.initialized()) {
                fwd.init(floorlog2(N), false);
            }
            memcpy(out, in, sizeof(cmplxT<T>) * N);
            fwd.permute(out);
            fwd.calc(out);
        } else {
            if (!rev.initialized()) {
                rev.init(floorlog2(N), true);
            }
            memcpy(out, in, sizeof(cmplxT<T>) * N);
            rev.permute(out);
            rev.calc(out);
        }
    }
   
protected:
    
    enum fft_permutation_type {
        FF_FFT_PERM_DEFAULT,
        FF_FFT_PERM_SWAP_LSBS,
        FF_FFT_PERM_AVX,
    };
    
    class FFTContext
    {
    protected:
        int         _nbits;
        bool        _inverse;
        uint16_t    *_revtab;
        uint32_t    *_revtab32;
        FFTComplex  *_tmpbuf;
        bool        _initialized;

    public:

        FFTContext(AVFFT *owner)
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
            if (_revtab) free(_revtab);
            _revtab = NULL;
            if (_revtab32) free(_revtab32);
            _revtab32 = NULL;
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

        int init(int nb, bool inverse)
        {
            int i, j, n;

            _revtab = NULL;
            _revtab32 = NULL;

            if (nb < 2 || nb > 17)
                goto fail;
            _nbits = nb;
            n = 1 << _nbits;
            _inverse = inverse;
           
            _tmpbuf = callocT<FFTComplex>(n);
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

        void permute(FFTComplex *z)
        {
            int j, np = 1 << _nbits;
            const uint16_t *rt = _revtab;
            const uint32_t *rt32 = _revtab32;
            FFTComplex *tmp = _tmpbuf;
            
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
        
        void calc(FFTComplex *z)
        {
            if (initialized())
                fft_dispatch[_nbits-2](z);
            else
                fprintf(stderr, "WARNING: init() was not called before calc() as required. Abort.\n");
        }

    protected:
        
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
    T r0=a0.re, i0=a0.im, r1=a1.re, i1=a1.im;\
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
static void name(FFTComplex *z, const zfloat *wre, unsigned int n)\
{\
    T t1, t2, t3, t4, t5, t6;\
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
        static void fft##n(FFTComplex *z)\
        {\
            fft##n2(z);\
            fft##n4(z+n4*2);\
            fft##n4(z+n4*3);\
            pass(z,ff_cos_##n,n4/2);\
        }
      
        static void fft4(FFTComplex *z)
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

        static void fft8(FFTComplex *z)
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

        static void fft16(FFTComplex *z)
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

        static constexpr void (* const fft_dispatch[])(FFTComplex*) = {
            fft4, fft8, fft16, fft32, fft64, fft128, fft256, fft512, fft1024,
            fft2048, fft4096, fft8192, fft16384, fft32768, fft65536, fft131072
        };
        enum fft_permutation_type fft_permutation;
    };

protected:
    
    FFTContext  fwd, rev;
  
};





