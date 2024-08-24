/** AVFFT v1.71 C++ wrapper class written by Dmitry Boldyrev
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

#include "avfft.h"

zfloat * ff_cos_tabs[18] = {
    NULL,
    NULL,
    NULL,
    NULL,
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

zfloat * ff_sin_tabs[18] = {
    NULL,
    NULL,
    NULL,
    NULL,
    ff_sin_16,
    ff_sin_32,
    ff_sin_64,
    ff_sin_128,
    ff_sin_256,
    ff_sin_512,
    ff_sin_1024,
    ff_sin_2048,
    ff_sin_4096,
    ff_sin_8192,
    ff_sin_16384,
    ff_sin_32768,
    ff_sin_65536,
    ff_sin_131072,
};

zfloat * get_ff_cos_tab(int size) {
    int index = floorlog2(size);
    return ff_cos_tabs[index];
}

zfloat *get_ff_sin_tab(int size) {
    int index = floorlog2(size);
    return ff_sin_tabs[index];
}

void init_ff_cosin_tabs(int index)
{
    int i;
    int m = 1<<index;
    zfloat freq = 2*M_PI/(zfloat)m;
    zfloat *ctab = ff_cos_tabs[index];
    zfloat *stab = ff_sin_tabs[index];
    for(i=0; i<=m/4; i++) {
        ctab[i] = cos(i*freq);
        stab[i] = sin(i*freq);
    }
    for(i=1; i<m/4; i++) {
        ctab[m/2-i] = ctab[i];
        stab[m/2-i] = stab[i];
    }
}


#define INIT_FF_COSIN_TABS_FUNC(index, size)   \
void init_ff_cosin_tabs_ ##size (void)  \
{                                              \
init_ff_cosin_tabs(index);                 \
}

INIT_FF_COSIN_TABS_FUNC(4, 16)
INIT_FF_COSIN_TABS_FUNC(5, 32)
INIT_FF_COSIN_TABS_FUNC(6, 64)
INIT_FF_COSIN_TABS_FUNC(7, 128)
INIT_FF_COSIN_TABS_FUNC(8, 256)
INIT_FF_COSIN_TABS_FUNC(9, 512)
INIT_FF_COSIN_TABS_FUNC(10, 1024)
INIT_FF_COSIN_TABS_FUNC(11, 2048)
INIT_FF_COSIN_TABS_FUNC(12, 4096)
INIT_FF_COSIN_TABS_FUNC(13, 8192)
INIT_FF_COSIN_TABS_FUNC(14, 16384)
INIT_FF_COSIN_TABS_FUNC(15, 32768)
INIT_FF_COSIN_TABS_FUNC(16, 65536)
INIT_FF_COSIN_TABS_FUNC(17, 131072)

typedef struct CoSinTabsInitOnce {
    void (*func)(void);
    bool loaded;
} CoSinTabsInitOnce;

CoSinTabsInitOnce cosin_tabs_init_once[18] = {
    { NULL, false },
    { NULL, false },
    { NULL, false },
    { NULL, false },
    { init_ff_cosin_tabs_16, false },
    { init_ff_cosin_tabs_32, false },
    { init_ff_cosin_tabs_64, false },
    { init_ff_cosin_tabs_128, false },
    { init_ff_cosin_tabs_256, false },
    { init_ff_cosin_tabs_512, false },
    { init_ff_cosin_tabs_1024, false },
    { init_ff_cosin_tabs_2048, false },
    { init_ff_cosin_tabs_4096, false },
    { init_ff_cosin_tabs_8192, false },
    { init_ff_cosin_tabs_16384, false },
    { init_ff_cosin_tabs_32768, false },
    { init_ff_cosin_tabs_65536, false },
    { init_ff_cosin_tabs_131072, false },
};

void ff_init_ff_cosin_tabs(int index)
{
    if (!cosin_tabs_init_once[index].loaded) {
        //LOG("generating sin + cos tabs for %d", 1<<index);
        cosin_tabs_init_once[index].func();
        cosin_tabs_init_once[index].loaded = true;
    }
}

#define DCLR_ALIGNED(n,t,v) t __attribute__ ((aligned (n))) v
#define DCLR_COSINTABLE(size) DCLR_ALIGNED(32, zfloat, ff_cos_##size)[size/2]; \
DCLR_ALIGNED(32, zfloat, ff_sin_##size)[size/2]


DCLR_COSINTABLE(16);
DCLR_COSINTABLE(32);
DCLR_COSINTABLE(64);
DCLR_COSINTABLE(128);
DCLR_COSINTABLE(256);
DCLR_COSINTABLE(512);
DCLR_COSINTABLE(1024);
DCLR_COSINTABLE(2048);
DCLR_COSINTABLE(4096);
DCLR_COSINTABLE(8192);
DCLR_COSINTABLE(16384);
DCLR_COSINTABLE(32768);
DCLR_COSINTABLE(65536);
DCLR_COSINTABLE(131072);
