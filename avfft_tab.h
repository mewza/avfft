#pragma once

/**     AVFFT v1.0 C++ wrapper class written by Dmitry Boldyrev 
***     GITHUB: https://github.com/mewza
***     Email: subband@protonmail.com
***
***     It is based on AV MPEG AUDIO source code, so I am not really
***     going to officially claim any (C) to it, because I honestly,
***     just put it together into an aeasy to use C++ templated wrapper,
***     the great benefit I tested it, it works well with 32-bit and 64-bit
***     floating point single types or as a form of intrinsic SIMD vectors. 
**/

static bool cos_tabs_init_once_initialized = false;

typedef struct CosTabsInitOnce {
    void (*func)(void);
    bool control;
} CosTabsInitOnce;

#define DECLARE_ALIGNED(n,t,v) t __attribute__ ((aligned (n))) v
#define COSTABLE(size) static DECLARE_ALIGNED(32, double, ff_cos_##size)[size/2]

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

static constexpr double * ff_cos_tabs[18] = {
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
#define _sqrthalf   M_SQRT1_2

#define CMUL(dre, dim, are, aim, bre, bim) do { \
    (dre) = (are) * (bre) - (aim) * (bim);  \
    (dim) = (are) * (bim) + (aim) * (bre);  \
} while (0)

static void init_ff_cos_tabs(int index)
{
    int i;
    int m = 1<<index;
    double freq = 2.*M_PI/(double)m;
    double *tab = ff_cos_tabs[index];
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

static constexpr CosTabsInitOnce cos_tabs_init_once_const[18] = {
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

static CosTabsInitOnce cos_tabs_init_once[18];

void ff_init_ff_cos_tabs(int index)
{
    if (!cos_tabs_init_once[index].control) {
        cos_tabs_init_once[index].func();
        cos_tabs_init_once[index].control = true;
    }
}
