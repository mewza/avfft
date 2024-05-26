#pragma once

/**     AVFFT v1.0 C++ wrapper class written by Dmitry Boldyrev 
***     GITHUB: https://github.com/mewza
***     Email: subband@protonmail.com
***
***     It is based on AV MPEG AUDIO source code, so I am not really
***     going to officially claim any (C) to it, because I honestly,
***     just put it together into an aeasy to use C++ templated class,
***    the great benefit I tested it, it works well with 32-bit and 64-bit
***    floating point single types or as a form of intrinsic SIMD vectors. 
*** 
***     NOTE: This one unlike WDLFFT actually worked and produced a spectrum
***     I could recognize, I think WDL has some sort of bug in permutation
***     but AVFFT seems to be working in a more similar way to PFFFT.
**/

#include "const1.h"
#include "avfft_tab.h"

template <typename T>
class AVFFT
{
    static constexpr int floorlog2(int x) {
        return (x == 1) ? 0 : 1 + floorlog2(x >> 1);
    }
    
    typedef cmplxT<T> FFTComplex;

public:
    
    AVFFT() {
        init_forward = false;
        init_reverse = false;
        
        if (!cos_tabs_init_once_initialized) {
            for (int i=0; i<18; i++) {
                cos_tabs_init_once[i] = cos_tabs_init_once_const[i];
            }
            cos_tabs_init_once_initialized = true;
        }
    }
    
    ~AVFFT() {
        if (init_forward)
            ff_fft_end(&ctx_fwd);
        if (init_reverse)
            ff_fft_end(&ctx_rev);
    }
    
    void real_fft_forward(const void *in, void *out, int size, bool scale = false)
    {
        const int bitSize = floorlog2(size);
        cmplxT<T> fft[size];
        const T *tin = (const T *)in;
        cmplxT<T> *cout = (cmplxT<T> *)out;
        T scv = 1.;
        
        if (!init_forward) {
            ff_fft_init(&ctx_fwd, bitSize, 0);
            init_forward = true;
        }
        
        if (scale) scv = 1./(double)size;
        for (int i=0; i<size; i++)
            cout[i] = cmplxT<T>(tin[i] * scv, 0.0);
        
        ctx_fwd.fft_permute(&ctx_fwd, cout);
        ctx_fwd.fft_calc(&ctx_fwd, cout);
        
        // memcpy(out, fft, size * sizeof(cmplxT<T>));
        /* DEBUG */
        //  if constexpr( std::is_same_v<T, mssFloat8> )
        //      PRINT_D8(fft[2].re);
        //  if constexpr( std::is_same_v<T, mssFloat4> )
        //      PRINT_D4(fft[2].re);
        //  if constexpr( std::is_same_v<T, mssFloat2> )
        //      PRINT_D2(fft[2].re);
        //  if constexpr( std::is_same_v<T, mssFloat> )
        //      PRINT_D1(fft[2].re);
    }
    
    void real_fft_reverse(const void *in, void *out, int size)
    {
        const int bitSize = floorlog2(size);
        const cmplxT<T> *cin = (const cmplxT<T> *)in;
        T *tout = (T *)out;
        cmplxT<T> fft[size];
        
        memcpy(fft, cin, sizeof(cmplxT<T>) * size);
     
        if (!init_reverse) {
            ff_fft_init(&ctx_rev, bitSize, 1);
            init_reverse = true;
        }
        ctx_fwd.fft_permute(&ctx_rev, fft);
        ctx_rev.fft_calc(&ctx_rev, fft);
        
        for (int i=0; i<size; i++)
            tout[i] = fft[i].re;
        
        /* DEBUG */
        //  if constexpr( std::is_same_v<T, mssFloat8> )
        //     PRINT_D8(out[2] );
        // if constexpr( std::is_same_v<T, mssFloat4> )
        //     PRINT_D4(out[2]);
        // if constexpr( std::is_same_v<T, mssFloat2> )
        //    PRINT_D2(out[2]);
        // if constexpr( std::is_same_v<T, mssFloat> )
        //    PRINT_D1(out[2]);
        
    }
    
    void cmplx_fft_forward(const cmplxT<T> *in, cmplxT<T> *out, int size)
    {
        const int bitSize = floorlog2(size);
        
        if (!init_forward) {
            ff_fft_init(&ctx_fwd, bitSize, 0);
            init_forward = true;
        }
        memcpy(out, in, sizeof(cmplxT<T>));
        
        ctx_fwd.fft_permute(&ctx_fwd, out);
        ctx_fwd.fft_calc(&ctx_fwd, out);
    }
    void cmplx_fft_reverse(const cmplxT<T> *in, cmplxT<T> *out,  int size)
    {
        const int bitSize = floorlog2(size);
        
        if (!init_reverse) {
            ff_fft_init(&ctx_rev, bitSize, 1);
            init_reverse = true;
        }
        memcpy(out, in, sizeof(cmplxT<T>));
        ctx_rev.fft_permute(&ctx_rev, out);
        ctx_rev.fft_calc(&ctx_rev, out);
    }
    
    
protected:
    
    enum fft_permutation_type {
        FF_FFT_PERM_DEFAULT,
        FF_FFT_PERM_SWAP_LSBS,
        FF_FFT_PERM_AVX,
    };

    enum mdct_permutation_type {
        FF_MDCT_PERM_NONE,
        FF_MDCT_PERM_INTERLEAVE,
    };

    struct FFTContext {
        int         nbits;
        int         inverse;
        uint16_t    *revtab;
        FFTComplex  *tmp_buf;
        int         mdct_size;  /* size of MDCT (i.e. number of input data * 2) */
        int         mdct_bits;  /* n = 2^nbits */
        uint32_t    *revtab32;
       
        /* pre/post rotation tables */
        T           *tcos;
        T           *tsin;
       
        /**
         * Do the permutation needed BEFORE calling fft_calc().
         */
        void (*fft_permute)(struct FFTContext *s, FFTComplex *z);
        /**
         * Do a complex FFT with the parameters defined in ff_fft_init(). The
         * input data must be permuted before. No 1.0/sqrt(n) normalization is done.
         */
       
        void (*fft_calc)(struct FFTContext *s, FFTComplex *z);
        void (*imdct_calc)(struct FFTContext *s, T *output, const T *input);
        void (*imdct_half)(struct FFTContext *s, T *output, const T *input);
        void (*mdct_calc)(struct FFTContext *s, T *output, const T *input);
        
        enum fft_permutation_type fft_permutation;
        enum mdct_permutation_type mdct_permutation;
    };
    
    static constexpr int avx_tab[] = {
        0, 4, 1, 5, 8, 12, 9, 13, 2, 6, 3, 7, 10, 14, 11, 15
    };

    FFTContext  ctx_fwd, ctx_rev;
    bool init_forward, init_reverse;
    
    static int split_radix_permutation(int i, int n, int inverse)
    {
        int m;
        if(n <= 2) return i&1;
        m = n >> 1;
        if(!(i&m))            return split_radix_permutation(i, m, inverse)*2;
        m >>= 1;
        if(inverse == !(i&m)) return split_radix_permutation(i, m, inverse)*4 + 1;
        else                  return split_radix_permutation(i, m, inverse)*4 - 1;
    }


    static int is_second_half_of_fft32(int i, int n)
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

    static void fft_perm_avx(FFTContext *s)
    {
        int i;
        int n = 1 << s->nbits;

        for (i = 0; i < n; i += 16) {
            int k;
            if (is_second_half_of_fft32(i, n)) {
                for (k = 0; k < 16; k++)
                    s->revtab[-split_radix_permutation(i + k, n, s->inverse) & (n - 1)] =
                        i + avx_tab[k];

            } else {
                for (k = 0; k < 16; k++) {
                    int j = i + k;
                    j = (j & ~7) | ((j >> 1) & 3) | ((j << 2) & 4);
                    s->revtab[-split_radix_permutation(i + k, n, s->inverse) & (n - 1)] = j;
                }
            }
        }
    }

     int ff_fft_init(FFTContext *s, int nbits, int inverse)
    {
        int i, j, n;

        s->revtab = NULL;
        s->revtab32 = NULL;

        if (nbits < 2 || nbits > 17)
            goto fail;
        s->nbits = nbits;
        n = 1 << nbits;

        if (nbits <= 16) {
            s->revtab = (uint16_t*)malloc(n * sizeof(uint16_t));
            if (!s->revtab)
                goto fail;
        } else {
            s->revtab32 = (uint32_t*)malloc(n * sizeof(uint32_t));
            if (!s->revtab32)
                goto fail;
        }
        s->tmp_buf = (FFTComplex*)malloc(n * sizeof(FFTComplex));
        if (!s->tmp_buf)
            goto fail;
        s->inverse = inverse;
         s->fft_permutation = FF_FFT_PERM_DEFAULT; //FF_FFT_PERM_SWAP_LSBS; //FF_FFT_PERM_DEFAULT;

        s->fft_permute = fft_permute_c;
        s->fft_calc    = fft_calc_c;
    #if CONFIG_MDCT
        s->imdct_calc  = ff_imdct_calc_c;
        s->imdct_half  = ff_imdct_half_c;
        s->mdct_calc   = ff_mdct_calc_c;
    #endif

    #if ARCH_AARCH64
        ff_fft_init_aarch64(s);
    #elif ARCH_ARM
        ff_fft_init_arm(s);
    #elif ARCH_PPC
        ff_fft_init_ppc(s);
    #elif ARCH_X86
        ff_fft_init_x86(s);
    #endif
         
    #if HAVE_MIPSFPU
        ff_fft_init_mips(s);
    #endif
        for(j=4; j <= nbits; j++) {
            ff_init_ff_cos_tabs(j);
        }


    #define PROCESS_FFT_PERM_SWAP_LSBS(num) do {\
        for(i = 0; i < n; i++) {\
            int k;\
            j = i;\
            j = (j & ~3) | ((j >> 1) & 1) | ((j << 1) & 2);\
            k = -split_radix_permutation(i, n, s->inverse) & (n - 1);\
            s->revtab##num[k] = j;\
        } \
    } while(0);

    #define PROCESS_FFT_PERM_DEFAULT(num) do {\
        for(i = 0; i < n; i++) {\
            int k;\
            j = i;\
            k = -split_radix_permutation(i, n, s->inverse) & (n - 1);\
            s->revtab##num[k] = j;\
        } \
    } while(0);

    #define SPLIT_RADIX_PERMUTATION(num) do { \
        if (s->fft_permutation == FF_FFT_PERM_SWAP_LSBS) {\
            PROCESS_FFT_PERM_SWAP_LSBS(num) \
        } else {\
            PROCESS_FFT_PERM_DEFAULT(num) \
        }\
    } while(0);

        if (s->revtab)
            SPLIT_RADIX_PERMUTATION()
        if (s->revtab32)
            SPLIT_RADIX_PERMUTATION(32)

    #undef PROCESS_FFT_PERM_DEFAULT
    #undef PROCESS_FFT_PERM_SWAP_LSBS
    #undef SPLIT_RADIX_PERMUTATION
        

        return 0;
     fail:
        free(s->revtab);
        free(s->revtab32);
        free(s->tmp_buf);
        return -1;
    }

    static void fft_permute_c(FFTContext *s, FFTComplex *z)
    {
        int j, np;
        const uint16_t *revtab = s->revtab;
        const uint32_t *revtab32 = s->revtab32;
        np = 1 << s->nbits;
        /* TODO: handle split-radix permute in a more optimal way, probably in-place */
        if (revtab) {
            for(j=0;j<np;j++) s->tmp_buf[revtab[j]] = z[j];
        } else
            for(j=0;j<np;j++) s->tmp_buf[revtab32[j]] = z[j];

        for(j=0;j<np;j++) z[j] = s->tmp_buf[j];
//        memcpy(z, s->tmp_buf, np * sizeof(FFTComplex));
    }

     void ff_fft_end(FFTContext *s)
    {
        free(s->revtab);
        free(s->revtab32);
        free(s->tmp_buf);
    }

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
    // this is slightly slower for small data, but avoids store->load aliasing
    // for addresses separated by large powers of 2.
    #define BUTTERFLIES_BIG(a0,a1,a2,a3) {\
        T r0=a0.re, i0=a0.im, r1=a1.re, i1=a1.im;\
        BF(t3, t5, t5, t1);\
        BF(a2.re, a0.re, r0, t5);\
        BF(a3.im, a1.im, i1, t3);\
        BF(t4, t6, t2, t6);\
        BF(a3.re, a1.re, r1, t4);\
        BF(a2.im, a0.im, i0, t6);\
    }

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

    /* z[0...8n-1], w[1...2n-1] */
    #define PASS(name)\
    static void name(FFTComplex *z, const double *wre, unsigned int n)\
    {\
        T t1, t2, t3, t4, t5, t6;\
        int o1 = 2*n;\
        int o2 = 4*n;\
        int o3 = 6*n;\
        const double *wim = wre+o1;\
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
//#if !CONFIG_SMALL
    #undef BUTTERFLIES
    #define BUTTERFLIES BUTTERFLIES_BIG
    PASS(pass_big)
//#endif

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
        TRANSFORM(z[1],z[3],z[5],z[7],_sqrthalf,_sqrthalf);
    }

//#if !CONFIG_SMALL
    static void fft16(FFTComplex *z)
    {
        T t1, t2, t3, t4, t5, t6;
        T cos_16_1 = ff_cos_16[1];
        T cos_16_3 = ff_cos_16[3];

        fft8(z);
        fft4(z+8);
        fft4(z+12);

        TRANSFORM_ZERO(z[0],z[4],z[8],z[12]);
        TRANSFORM(z[2],z[6],z[10],z[14],_sqrthalf,_sqrthalf);
        TRANSFORM(z[1],z[5],z[9],z[13],cos_16_1,cos_16_3);
        TRANSFORM(z[3],z[7],z[11],z[15],cos_16_3,cos_16_1);
    }
//#else
//  DECL_FFT(16,8,4)
//#endif
    DECL_FFT(32,16,8)
    DECL_FFT(64,32,16)
    DECL_FFT(128,64,32)
    DECL_FFT(256,128,64)
    DECL_FFT(512,256,128)
//#if !CONFIG_SMALL
    #define pass pass_big
//#endif
    DECL_FFT(1024,512,256)
    DECL_FFT(2048,1024,512)
    DECL_FFT(4096,2048,1024)
    DECL_FFT(8192,4096,2048)
    DECL_FFT(16384,8192,4096)
    DECL_FFT(32768,16384,8192)
    DECL_FFT(65536,32768,16384)
    DECL_FFT(131072,65536,32768)

    constexpr static void (* const fft_dispatch[])(FFTComplex*) = {
        fft4, fft8, fft16, fft32, fft64, fft128, fft256, fft512, fft1024,
        fft2048, fft4096, fft8192, fft16384, fft32768, fft65536, fft131072
    };

    static void fft_calc_c(FFTContext *s, FFTComplex *z)
    {
        fft_dispatch[s->nbits-2](z);
    }
};
