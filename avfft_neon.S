/**
 **  ARM NEON optimizations of AVFFT
 **
 **  Copyright (c) 2009 Mans Rullgard <mans@mansr.com>
 **  Copyright (c) 2009 Naotoshi Nojiri
 **  Copyright (c) 2014 Janne Grunau <janne-libav@jannau.net>
 **/

#define HAVE_AS_FUNC 0
#define CONFIG_PIC 1
#define HAVE_SECTION_DATA_REL_RO 0
#define EXTERN_ASM _

#ifdef __ELF__
#   define ELF
#else
#   define ELF #
#endif

#if HAVE_AS_FUNC
#   define FUNC
#else
#   define FUNC #
#endif

#ifndef __has_feature
#   define __has_feature(x) 0
#endif

/* Support macros for
 *   - Armv8.3-A Pointer Authentication and
 *   - Armv8.5-A Branch Target Identification
 * features which require emitting a .note.gnu.property section with the
 * appropriate architecture-dependent feature bits set.
 *
 * |AARCH64_SIGN_LINK_REGISTER| and |AARCH64_VALIDATE_LINK_REGISTER| expand to
 * PACIxSP and AUTIxSP, respectively. |AARCH64_SIGN_LINK_REGISTER| should be
 * used immediately before saving the LR register (x30) to the stack.
 * |AARCH64_VALIDATE_LINK_REGISTER| should be used immediately after restoring
 * it. Note |AARCH64_SIGN_LINK_REGISTER|'s modifications to LR must be undone
 * with |AARCH64_VALIDATE_LINK_REGISTER| before RET. The SP register must also
 * have the same value at the two points. For example:
 *
 *   .global f
 *   f:
 *     AARCH64_SIGN_LINK_REGISTER
 *     stp x29, x30, [sp, #-96]!
 *     mov x29, sp
 *     ...
 *     ldp x29, x30, [sp], #96
 *     AARCH64_VALIDATE_LINK_REGISTER
 *     ret
 *
 * |AARCH64_VALID_CALL_TARGET| expands to BTI 'c'. Either it, or
 * |AARCH64_SIGN_LINK_REGISTER|, must be used at every point that may be an
 * indirect call target. In particular, all symbols exported from a file must
 * begin with one of these macros. For example, a leaf function that does not
 * save LR can instead use |AARCH64_VALID_CALL_TARGET|:
 *
 *   .globl return_zero
 *   return_zero:
 *     AARCH64_VALID_CALL_TARGET
 *     mov x0, #0
 *     ret
 *
 * A non-leaf function which does not immediately save LR may need both macros
 * because |AARCH64_SIGN_LINK_REGISTER| appears late. For example, the function
 * may jump to an alternate implementation before setting up the stack:
 *
 *   .globl with_early_jump
 *   with_early_jump:
 *     AARCH64_VALID_CALL_TARGET
 *     cmp x0, #128
 *     b.lt .Lwith_early_jump_128
 *     AARCH64_SIGN_LINK_REGISTER
 *     stp x29, x30, [sp, #-96]!
 *     mov x29, sp
 *     ...
 *     ldp x29, x30, [sp], #96
 *     AARCH64_VALIDATE_LINK_REGISTER
 *     ret
 *
 *  .Lwith_early_jump_128:
 *     ...
 *     ret
 *
 * These annotations are only required with indirect calls. Private symbols that
 * are only the target of direct calls do not require annotations. Also note
 * that |AARCH64_VALID_CALL_TARGET| is only valid for indirect calls (BLR), not
 * indirect jumps (BR). Indirect jumps in assembly are supported through
 * |AARCH64_VALID_JUMP_TARGET|. Landing Pads which shall serve for jumps and
 * calls can be created using |AARCH64_VALID_JUMP_CALL_TARGET|.
 *
 * Although not necessary, it is safe to use these macros in 32-bit ARM
 * assembly. This may be used to simplify dual 32-bit and 64-bit files.
 *
 * References:
 * - "ELF for the Arm® 64-bit Architecture"
 *   https: *github.com/ARM-software/abi-aa/blob/master/aaelf64/aaelf64.rst
 * - "Providing protection for complex software"
 *   https://developer.arm.com/architectures/learn-the-architecture/providing-protection-for-complex-software
 */
#if defined(__ARM_FEATURE_BTI_DEFAULT) && (__ARM_FEATURE_BTI_DEFAULT == 1)
#   define GNU_PROPERTY_AARCH64_BTI (1 << 0)   // Has BTI
#   define AARCH64_VALID_CALL_TARGET hint #34  // BTI 'c'
#   define AARCH64_VALID_JUMP_TARGET hint #38  // BTI 'j'
#else
#   define GNU_PROPERTY_AARCH64_BTI 0          // No BTI
#   define AARCH64_VALID_CALL_TARGET
#   define AARCH64_VALID_JUMP_TARGET
#endif

#if defined(__ARM_FEATURE_PAC_DEFAULT)
#   if ((__ARM_FEATURE_PAC_DEFAULT & (1 << 0)) != 0) // authentication using key A
#       define AARCH64_SIGN_LINK_REGISTER      paciasp
#       define AARCH64_VALIDATE_LINK_REGISTER  autiasp
#   elif ((__ARM_FEATURE_PAC_DEFAULT & (1 << 1)) != 0) // authentication using key B
#       define AARCH64_SIGN_LINK_REGISTER      pacibsp
#       define AARCH64_VALIDATE_LINK_REGISTER  autibsp
#   else
#       error Pointer authentication defines no valid key!
#   endif
#   if ((__ARM_FEATURE_PAC_DEFAULT & (1 << 2)) != 0)
#       error Authentication of leaf functions is enabled but not supported in FFmpeg!
#   endif
#   define GNU_PROPERTY_AARCH64_PAC (1 << 1)
#else
#   define GNU_PROPERTY_AARCH64_PAC 0
#   define AARCH64_SIGN_LINK_REGISTER
#   define AARCH64_VALIDATE_LINK_REGISTER
#endif


#if (GNU_PROPERTY_AARCH64_BTI != 0 || GNU_PROPERTY_AARCH64_PAC != 0) && defined(__ELF__)
        .pushsection .note.gnu.property, "a"
        .balign 8
        .long 4
        .long 0x10
        .long 0x5
        .asciz "GNU"
        .long 0xc0000000 /* GNU_PROPERTY_AARCH64_FEATURE_1_AND */
        .long 4
        .long (GNU_PROPERTY_AARCH64_BTI | GNU_PROPERTY_AARCH64_PAC)
        .long 0
        .popsection
#endif

.macro  function name, export=0, align=2
    .macro endfunc
ELF     .size   \name, . - \name
FUNC    .endfunc
        .purgem endfunc
    .endm
        .text
        .align          \align
    .if \export
        .global EXTERN_ASM\name
ELF     .type   EXTERN_ASM\name, %function
FUNC    .func   EXTERN_ASM\name
EXTERN_ASM\name:
        AARCH64_VALID_CALL_TARGET
    .else
ELF     .type   \name, %function
FUNC    .func   \name
\name:
    .endif
.endm

.macro  const   name, align=2, relocate=0
    .macro endconst
ELF     .size   \name, . - \name
        .purgem endconst
    .endm
#if HAVE_SECTION_DATA_REL_RO
.if \relocate
        .section        .data.rel.ro
.else
        .section        .rodata
.endif
#elif defined(_WIN32)
        .section        .rdata
#elif !defined(__MACH__)
        .section        .rodata
#else
        .const_data
#endif
        .align          \align
\name:
.endm

.macro  movrel rd, val, offset=0
#if CONFIG_PIC && defined(__APPLE__)
    .if \offset < 0
        adrp            \rd, \val@PAGE
        add             \rd, \rd, \val@PAGEOFF
        sub             \rd, \rd, -(\offset)
    .else
        adrp            \rd, \val+(\offset)@PAGE
        add             \rd, \rd, \val+(\offset)@PAGEOFF
    .endif
#elif CONFIG_PIC && defined(_WIN32)
    .if \offset < 0
        adrp            \rd, \val
        add             \rd, \rd, :lo12:\val
        sub             \rd, \rd, -(\offset)
    .else
        adrp            \rd, \val+(\offset)
        add             \rd, \rd, :lo12:\val+(\offset)
    .endif
#elif CONFIG_PIC
#   if __has_feature(hwaddress_sanitizer)
        adrp            \rd, :pg_hi21_nc:\val+(\offset)
#   else
        adrp            \rd, \val+(\offset)
#   endif
        add             \rd, \rd, :lo12:\val+(\offset)
#else
        ldr             \rd, =\val+\offset
#endif
.endm

#define GLUE(a, b) a ## b
#define JOIN(a, b) GLUE(a, b)
#define X(s) JOIN(EXTERN_ASM, s)

#define x18 do_not_use_x18
#define w18 do_not_use_w18


#define M_SQRT1_2 0.70710678118654752440

.macro transpose d0, d1, s0, s1
        trn1            \d0, \s0, \s1
        trn2            \d1, \s0, \s1
.endm


function fft4_neon, export=1
        AARCH64_VALID_JUMP_TARGET
        ld1             {v0.2s,v1.2s,v2.2s,v3.2s}, [x0]

        fadd            v4.2s,  v0.2s,  v1.2s   // r0+r1,i0+i1
        fsub            v6.2s,  v0.2s,  v1.2s   // r0-r1,i0-i1

        ext             v16.8b, v2.8b,  v3.8b,  #4
        ext             v17.8b, v3.8b,  v2.8b,  #4

        fadd            v5.2s,  v2.2s,  v3.2s   // i2+i3,r2+r3
        fsub            v7.2s,  v16.2s, v17.2s  // r3-r2,i2-i3

        fadd            v0.2s,  v4.2s,  v5.2s
        fsub            v2.2s,  v4.2s,  v5.2s
        fadd            v1.2s,  v6.2s,  v7.2s
        fsub            v3.2s,  v6.2s,  v7.2s

        st1             {v0.2s,v1.2s,v2.2s,v3.2s}, [x0]

        ret
endfunc

function fft8_neon, export=1
        AARCH64_VALID_JUMP_TARGET
        mov             x1,  x0
        ld1             {v0.2s, v1.2s, v2.2s, v3.2s},  [x0], #32
        ld1             {v16.2s,v17.2s,v18.2s,v19.2s}, [x0]
        ext             v22.8b, v2.8b,  v3.8b,  #4
        ext             v23.8b, v3.8b,  v2.8b,  #4
        fadd            v4.2s,  v16.2s, v17.2s           // r4+r5,i4+i5
        fadd            v5.2s,  v18.2s, v19.2s           // r6+r7,i6+i7
        fsub            v17.2s, v16.2s, v17.2s           // r4-r5,i4-i5
        fsub            v19.2s, v18.2s, v19.2s           // r6-r7,i6-i7
        rev64           v27.2s, v28.2s  // ???
        fadd            v20.2s, v0.2s,  v1.2s            // r0+r1,i0+i1
        fadd            v21.2s, v2.2s,  v3.2s            // r2+r3,i2+i3
        fmul            v26.2s, v17.2s, v28.2s           // -a2r*w,a2i*w
        ext             v6.8b,  v4.8b,  v5.8b,  #4
        ext             v7.8b,  v5.8b,  v4.8b,  #4
        fmul            v27.2s, v19.2s, v27.2s           // a3r*w,-a3i*w
        fsub            v23.2s, v22.2s, v23.2s           // i2-i3,r3-r2
        fsub            v22.2s, v0.2s,  v1.2s            // r0-r1,i0-i1
        fmul            v24.2s, v17.2s, v28.s[1]         // a2r*w,a2i*w
        fmul            v25.2s, v19.2s, v28.s[1]         // a3r*w,a3i*w
        fadd            v0.2s,  v20.2s, v21.2s
        fsub            v2.2s,  v20.2s, v21.2s
        fadd            v1.2s,  v22.2s, v23.2s
        rev64           v26.2s, v26.2s
        rev64           v27.2s, v27.2s
        fsub            v3.2s,  v22.2s, v23.2s
        fsub            v6.2s,  v6.2s,  v7.2s
        fadd            v24.2s, v24.2s, v26.2s  // a2r+a2i,a2i-a2r   t1,t2
        fadd            v25.2s, v25.2s, v27.2s  // a3r-a3i,a3i+a3r   t5,t6
        fadd            v7.2s,  v4.2s,  v5.2s
        fsub            v18.2s, v2.2s,  v6.2s
        ext             v26.8b, v24.8b, v25.8b, #4
        ext             v27.8b, v25.8b, v24.8b, #4
        fadd            v2.2s,  v2.2s,  v6.2s
        fsub            v16.2s, v0.2s,  v7.2s
        fadd            v5.2s,  v25.2s, v24.2s
        fsub            v4.2s,  v26.2s, v27.2s
        fadd            v0.2s,  v0.2s,  v7.2s
        fsub            v17.2s, v1.2s,  v5.2s
        fsub            v19.2s, v3.2s,  v4.2s
        fadd            v3.2s,  v3.2s,  v4.2s
        fadd            v1.2s,  v1.2s,  v5.2s

        st1             {v16.2s,v17.2s,v18.2s,v19.2s}, [x0]
        st1             {v0.2s, v1.2s, v2.2s, v3.2s},  [x1]

        ret
endfunc

function fft16_neon, export=1
        AARCH64_VALID_JUMP_TARGET
        mov             x1,  x0
        ld1             {v0.2s, v1.2s, v2.2s, v3.2s},  [x0], #32
        ld1             {v16.2s,v17.2s,v18.2s,v19.2s}, [x0], #32
        ext             v22.8b, v2.8b,  v3.8b,  #4
        ext             v23.8b, v3.8b,  v2.8b,  #4
        fadd            v4.2s,  v16.2s, v17.2s           // r4+r5,i4+i5
        fadd            v5.2s,  v18.2s, v19.2s           // r6+r7,i6+i7
        fsub            v17.2s, v16.2s, v17.2s           // r4-r5,i4-i5
        fsub            v19.2s, v18.2s, v19.2s           // r6-r7,i6-i7
        rev64           v27.2s, v28.2s  // ???
        fadd            v20.2s, v0.2s,  v1.2s            // r0+r1,i0+i1
        fadd            v21.2s, v2.2s,  v3.2s            // r2+r3,i2+i3
        fmul            v26.2s, v17.2s, v28.2s           // -a2r*w,a2i*w
        ext             v6.8b,  v4.8b,  v5.8b,  #4
        ext             v7.8b,  v5.8b,  v4.8b,  #4
        fmul            v27.2s, v19.2s, v27.2s           // a3r*w,-a3i*w
        fsub            v23.2s, v22.2s, v23.2s           // i2-i3,r3-r2
        fsub            v22.2s, v0.2s,  v1.2s            // r0-r1,i0-i1
        fmul            v24.2s, v17.2s, v28.s[1]         // a2r*w,a2i*w
        fmul            v25.2s, v19.2s, v28.s[1]         // a3r*w,a3i*w
        fadd            v0.2s,  v20.2s, v21.2s
        fsub            v2.2s,  v20.2s, v21.2s
        fadd            v1.2s,  v22.2s, v23.2s
        rev64           v26.2s, v26.2s
        rev64           v27.2s, v27.2s
        fsub            v3.2s,  v22.2s, v23.2s
        fsub            v6.2s,  v6.2s,  v7.2s
        fadd            v24.2s, v24.2s, v26.2s  // a2r+a2i,a2i-a2r   t1,t2
        fadd            v25.2s, v25.2s, v27.2s  // a3r-a3i,a3i+a3r   t5,t6
        fadd            v7.2s,  v4.2s,  v5.2s
        fsub            v18.2s, v2.2s,  v6.2s
        ld1             {v20.4s,v21.4s}, [x0], #32
        ld1             {v22.4s,v23.4s}, [x0], #32
        ext             v26.8b, v24.8b, v25.8b, #4
        ext             v27.8b, v25.8b, v24.8b, #4
        fadd            v2.2s,  v2.2s,  v6.2s
        fsub            v16.2s, v0.2s,  v7.2s
        fadd            v5.2s,  v25.2s, v24.2s
        fsub            v4.2s,  v26.2s, v27.2s
        transpose       v24.2d, v25.2d, v20.2d, v22.2d
        transpose       v26.2d, v27.2d, v21.2d, v23.2d
        fadd            v0.2s,  v0.2s,  v7.2s
        fsub            v17.2s, v1.2s,  v5.2s
        fsub            v19.2s, v3.2s,  v4.2s
        fadd            v3.2s,  v3.2s,  v4.2s
        fadd            v1.2s,  v1.2s,  v5.2s
        ext             v20.16b, v21.16b, v21.16b,  #4
        ext             v21.16b, v23.16b, v23.16b,  #4

        zip1            v0.2d,  v0.2d,  v1.2d   // {z[0],   z[1]}
        zip1            v1.2d,  v2.2d,  v3.2d   // {z[2],   z[3]}
        zip1            v2.2d,  v16.2d, v17.2d  // {z[o1],  z[o1+1]}
        zip1            v3.2d,  v18.2d, v19.2d  // {z[o1+2],z[o1+3]}

        // 2 x fft4
        transpose       v22.2d, v23.2d, v20.2d, v21.2d

        fadd            v4.4s,  v24.4s, v25.4s
        fadd            v5.4s,  v26.4s, v27.4s
        fsub            v6.4s,  v24.4s, v25.4s
        fsub            v7.4s,  v22.4s, v23.4s

        ld1             {v23.4s},  [x14]

        fadd            v24.4s, v4.4s,  v5.4s   // {z[o2+0],z[o2+1]}
        fsub            v26.4s, v4.4s,  v5.4s   // {z[o2+2],z[o2+3]}
        fadd            v25.4s, v6.4s,  v7.4s   // {z[o3+0],z[o3+1]}
        fsub            v27.4s, v6.4s,  v7.4s   // {z[o3+2],z[o3+3]}

        //fft_pass_neon_16
        rev64           v7.4s,  v25.4s
        fmul            v25.4s, v25.4s, v23.s[1]
        fmul            v7.4s,  v7.4s,  v29.4s
        fmla            v25.4s, v7.4s,  v23.s[3] // {t1a,t2a,t5a,t6a}

        zip1            v20.4s, v24.4s, v25.4s
        zip2            v21.4s, v24.4s, v25.4s
        fneg            v22.4s, v20.4s
        fadd            v4.4s,  v21.4s, v20.4s
        fsub            v6.4s,  v20.4s, v21.4s  // just the second half
        fadd            v5.4s,  v21.4s, v22.4s  // just the first half

        tbl             v4.16b, {v4.16b},        v30.16b // trans4_float
        tbl             v5.16b, {v5.16b,v6.16b}, v31.16b // trans8_float

        fsub            v20.4s, v0.4s,  v4.4s   // {z[o2],z[o2+1]}
        fadd            v16.4s, v0.4s,  v4.4s   // {z[0], z[1]}
        fsub            v22.4s, v2.4s,  v5.4s   // {z[o3],z[o3+1]}
        fadd            v18.4s, v2.4s,  v5.4s   // {z[o1],z[o1+1]}

//second half
        rev64           v6.4s,  v26.4s
        fmul            v26.4s, v26.4s, v23.s[2]
        rev64           v7.4s,  v27.4s
        fmul            v27.4s, v27.4s, v23.s[3]
        fmul            v6.4s,  v6.4s,  v29.4s
        fmul            v7.4s,  v7.4s,  v29.4s
        fmla            v26.4s, v6.4s,  v23.s[2] // {t1,t2,t5,t6}
        fmla            v27.4s, v7.4s,  v23.s[1] // {t1a,t2a,t5a,t6a}

        zip1            v24.4s, v26.4s, v27.4s
        zip2            v25.4s, v26.4s, v27.4s
        fneg            v26.4s, v24.4s
        fadd            v4.4s,  v25.4s, v24.4s
        fsub            v6.4s,  v24.4s, v25.4s  // just the second half
        fadd            v5.4s,  v25.4s, v26.4s  // just the first half

        tbl             v4.16b, {v4.16b},        v30.16b // trans4_float
        tbl             v5.16b, {v5.16b,v6.16b}, v31.16b // trans8_float

        fadd            v17.4s, v1.4s, v4.4s    // {z[2], z[3]}
        fsub            v21.4s, v1.4s, v4.4s    // {z[o2+2],z[o2+3]}
        fadd            v19.4s, v3.4s, v5.4s    // {z[o1+2],z[o1+3]}
        fsub            v23.4s, v3.4s, v5.4s    // {z[o3+2],z[o3+3]}

        st1             {v16.4s,v17.4s}, [x1], #32
        st1             {v18.4s,v19.4s}, [x1], #32
        st1             {v20.4s,v21.4s}, [x1], #32
        st1             {v22.4s,v23.4s}, [x1], #32

        ret
endfunc


const  trans4_float, align=4
        .byte    0,  1,  2,  3
        .byte    8,  9, 10, 11
        .byte    4,  5,  6,  7
        .byte   12, 13, 14, 15
endconst

const  trans8_float, align=4
        .byte   24, 25, 26, 27
        .byte    0,  1,  2,  3
        .byte   28, 29, 30, 31
        .byte    4,  5,  6,  7
endconst

function pass_neon, export=1
        sub             x6,  x2,  #1            // n - 1, loop counter
        lsl             x5,  x2,  #3            // 2 * n * sizeof FFTSample
        lsl             x1,  x2,  #4            // 2 * n * sizeof FFTComplex
        add             x5,  x4,  x5            // wim
        add             x3,  x1,  x2,  lsl #5   // 4 * n * sizeof FFTComplex
        add             x2,  x0,  x2,  lsl #5   // &z[o2]
        add             x3,  x0,  x3            // &z[o3]
        add             x1,  x0,  x1            // &z[o1]
        ld1             {v20.4s},[x2]           // {z[o2],z[o2+1]}
        ld1             {v22.4s},[x3]           // {z[o3],z[o3+1]}
        ld1             {v4.2s},  [x4], #8      // {wre[0],wre[1]}
        trn2            v25.2d, v20.2d, v22.2d
        sub             x5,  x5,  #4            // wim--
        trn1            v24.2d, v20.2d, v22.2d
        ld1             {v5.s}[0],  [x5], x7    // d5[0] = wim[-1]
        rev64           v7.4s,  v25.4s
        fmul            v25.4s, v25.4s, v4.s[1]
        ld1             {v16.4s}, [x0]          // {z[0],z[1]}
        fmul            v7.4s,  v7.4s,  v29.4s
        ld1             {v17.4s}, [x1]          // {z[o1],z[o1+1]}
        prfm            pldl1keep, [x2, #16]
        prfm            pldl1keep, [x3, #16]
        fmla            v25.4s, v7.4s,  v5.s[0] // {t1a,t2a,t5a,t6a}
        prfm            pldl1keep, [x0, #16]
        prfm            pldl1keep, [x1, #16]

        zip1            v20.4s, v24.4s, v25.4s
        zip2            v21.4s, v24.4s, v25.4s
        fneg            v22.4s, v20.4s
        fadd            v4.4s,  v21.4s, v20.4s
        fsub            v6.4s,  v20.4s, v21.4s  // just the second half
        fadd            v5.4s,  v21.4s, v22.4s  // just the first half

        tbl             v4.16b, {v4.16b},        v30.16b // trans4_float
        tbl             v5.16b, {v5.16b,v6.16b}, v31.16b // trans8_float

        fadd            v20.4s, v16.4s, v4.4s
        fsub            v22.4s, v16.4s, v4.4s
        fadd            v21.4s, v17.4s, v5.4s
        st1             {v20.4s}, [x0], #16     // {z[0], z[1]}
        fsub            v23.4s, v17.4s, v5.4s

        st1             {v21.4s}, [x1], #16     // {z[o1],z[o1+1]}
        st1             {v22.4s}, [x2], #16     // {z[o2],z[o2+1]}
        st1             {v23.4s}, [x3], #16     // {z[o3],z[o3+1]}
1:
        ld1             {v20.4s},[x2]    // {z[o2],z[o2+1]}
        ld1             {v22.4s},[x3]    // {z[o3],z[o3+1]}
        ld1             {v4.2s}, [x4], #8       // {wre[0],wre[1]}
        transpose       v26.2d, v27.2d, v20.2d, v22.2d
        ld1             {v5.2s}, [x5], x7       // {wim[-1],wim[0]}
        rev64           v6.4s,  v26.4s
        fmul            v26.4s, v26.4s, v4.s[0]
        rev64           v7.4s,  v27.4s
        fmul            v27.4s, v27.4s, v4.s[1]
        fmul            v6.4s,  v6.4s,  v29.4s
        fmul            v7.4s,  v7.4s,  v29.4s
        ld1             {v16.4s},[x0]           // {z[0],z[1]}
        fmla            v26.4s, v6.4s,  v5.s[1] // {t1,t2,t5,t6}
        fmla            v27.4s, v7.4s,  v5.s[0] // {t1a,t2a,t5a,t6a}
        ld1             {v17.4s},[x1]           // {z[o1],z[o1+1]}

        subs            x6,  x6,  #1            // n--

        zip1            v20.4s, v26.4s, v27.4s
        zip2            v21.4s, v26.4s, v27.4s
        fneg            v22.4s, v20.4s
        fadd            v4.4s,  v21.4s, v20.4s
        fsub            v6.4s,  v20.4s, v21.4s  // just the second half
        fadd            v5.4s,  v21.4s, v22.4s  // just the first half

        tbl             v4.16b, {v4.16b},        v30.16b // trans4_float
        tbl             v5.16b, {v5.16b,v6.16b}, v31.16b // trans8_float

        fadd            v20.4s, v16.4s, v4.4s
        fsub            v22.4s, v16.4s, v4.4s
        fadd            v21.4s, v17.4s, v5.4s
        st1             {v20.4s}, [x0], #16     // {z[0], z[1]}
        fsub            v23.4s, v17.4s, v5.4s

        st1             {v21.4s}, [x1], #16     // {z[o1],z[o1+1]}
        st1             {v22.4s}, [x2], #16     // {z[o2],z[o2+1]}
        st1             {v23.4s}, [x3], #16     // {z[o3],z[o3+1]}
        b.ne            1b

        ret
endfunc

.macro  def_fft n, n2, n4
function _fft\n\()_neon, align=6
        AARCH64_VALID_JUMP_TARGET
        AARCH64_SIGN_LINK_REGISTER
        stp             x28, x30, [sp, #-16]!
        add             x28, x0,  #\n4*2*8
        bl              _fft\n2\()_neon
        mov             x0,  x28
        bl              _fft\n4\()_neon
        add             x0,  x28, #\n4*1*8
        bl              _fft\n4\()_neon
        sub             x0,  x28, #\n4*2*8
        ldp             x28, x30, [sp], #16
        AARCH64_VALIDATE_LINK_REGISTER
        movrel          x4,  X(ff_cos_\n)
        mov             x2,  #\n4>>1
        b               _pass_neon
endfunc
.endm

        def_fft    32,    16,     8
        def_fft    64,    32,    16
        def_fft   128,    64,    32
        def_fft   256,   128,    64
        def_fft   512,   256,   128
        def_fft  1024,   512,   256
        def_fft  2048,  1024,   512
        def_fft  4096,  2048,  1024
        def_fft  8192,  4096,  2048
        def_fft 16384,  8192,  4096
        def_fft 32768, 16384,  8192
        def_fft 65536, 32768, 16384

function ff_fft_calc_neon, export=1
        prfm            pldl1keep, [x1]
        movrel          x10, trans4_float
        ldr             w2,  [x0]
        movrel          x11, trans8_float
        sub             w2,  w2,  #2
        movrel          x3,  fft_tab_neon
        ld1             {v30.16b}, [x10]
        mov             x7,  #-8
        movrel          x12, pmmp
        ldr             x3,  [x3, x2, lsl #3]
        movrel          x13, mppm
        movrel          x14, X(ff_cos_16)
        ld1             {v31.16b}, [x11]
        mov             x0,  x1
        ld1             {v29.4s},  [x12]         // pmmp
        ld1             {v28.4s},  [x13]
        br              x3
endfunc

function ff_fft_permute_neon, export=1
        mov             x6,  #1
        ldr             w2,  [x0]       // nbits
        ldr             x3,  [x0, #16]  // tmp_buf
        ldr             x0,  [x0, #8]   // revtab
        lsl             x6,  x6, x2
        mov             x2,  x6
1:
        ld1             {v0.2s,v1.2s}, [x1], #16
        ldr             w4,  [x0], #4
        uxth            w5,  w4
        lsr             w4,  w4,  #16
        add             x5,  x3,  x5,  lsl #3
        add             x4,  x3,  x4,  lsl #3
        st1             {v0.2s}, [x5]
        st1             {v1.2s}, [x4]
        subs            x6,  x6, #2
        b.gt            1b

        sub             x1,  x1,  x2,  lsl #3
1:
        ld1             {v0.4s,v1.4s}, [x3], #32
        st1             {v0.4s,v1.4s}, [x1], #32
        subs            x2,  x2,  #4
        b.gt            1b

        ret
endfunc

const   fft_tab_neon, relocate=1
        .quad _fft4_neon
        .quad _fft8_neon
        .quad _fft16_neon
        .quad _fft32_neon
        .quad _fft64_neon
        .quad _fft128_neon
        .quad _fft256_neon
        .quad _fft512_neon
        .quad _fft1024_neon
        .quad _fft2048_neon
        .quad _fft4096_neon
        .quad _fft8192_neon
        .quad _fft16384_neon
        .quad _fft32768_neon
        .quad _fft65536_neon
endconst

const   pmmp, align=4
        .float          +1.0, -1.0, -1.0, +1.0
endconst

const   mppm, align=4
        .float          -M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, -M_SQRT1_2
endconst
