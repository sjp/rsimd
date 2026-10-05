/* Shared prelude of every tier translation unit (src/tier_<tier>.c).
 *
 * A tier TU defines RSIMD_TIER as one of none, sse2, avx2, avx512, neon, sve
 * or sve2, includes this file once and then the kernel sources. The TU is
 * compiled with that tier's compiler flags (-mavx2 -mfma for avx2 ...), so the
 * SIMDe headers included here select native code generation for the tier.
 *
 * This file provides:
 *   - tier identification: RSIMD_TIER_IS(t), RSIMD_TIER_STRING;
 *   - name mangling: RSIMD_KERNEL(sum_f64) expands to rsimd_sum_f64_<tier>;
 *   - a guard that fails the build if the TU lacks its tier's flags;
 *   - the vector layer, so a kernel is written once for every tier.
 *
 * Vector layer
 * ------------
 * Types, one per element type: rsimd_vf64 (double), rsimd_vi32 (int32_t),
 * rsimd_vi64 (int64_t); comparison masks rsimd_mf64, rsimd_mi32, rsimd_mi64;
 * lane predicates rsimd_p64 (for f64 and i64 lanes) and rsimd_p32.
 *
 * Lane counts: RSIMD_LANES_64 and RSIMD_LANES_32 are valid on every tier but
 * are run-time values on the SVE tiers, where RSIMD_VLA is defined. The
 * fixed-width tiers also define the constants RSIMD_WIDTH_F64,
 * RSIMD_WIDTH_I64 and RSIMD_WIDTH_I32.
 *
 * Operations are static inline functions named rsimd_<type>_<op>, each with
 * an upper-case alias (rsimd_vf64_add is also RSIMD_VF64_ADD):
 *
 *   all vector types  loadu storeu loadu_p storeu_p set1 zero add sub mul
 *                     min max and or xor andnot cmp_eq cmp_gt cmp_lt blend
 *                     is_na reduce_add reduce_min reduce_max
 *   rsimd_vf64 also   div fma abs neg sqrt floor ceil trunc rint cmp_ne cmp_le cmp_ge
 *                     is_nan as_vi64 loadu_i32 loadu_i32_p storeu_i32
 *                     storeu_i32_p uzp_even uzp_odd recip_approx rsqrt_approx
 *   rsimd_vi32 also   mulhi sll srl sra loadu_u8 loadu_u8_p storeu_u8
 *                     storeu_u8_p inc
 *   rsimd_vi64 also   as_vf64 loadu_i32 loadu_i32_p storeu_i32 storeu_i32_p
 *                     inc
 *   mask types        and or not andnot any all count; rsimd_mf64_to_mi64
 *                     and rsimd_mi64_to_mf64 convert between the f64 and
 *                     i64 masks
 *   predicate types   while true count
 *
 * Semantics shared by all tiers:
 *   - Loads and stores are unaligned (R vectors are only 8-byte aligned).
 *   - fma(a, b, c) is a * b + c with a single rounding, on every tier.
 *   - recip_approx(a) and rsqrt_approx(a) are 1 / a and 1 / sqrt(a)
 *     within a relative error of 2^-22, for 2^-1022 <= |a| < 2^1022 (and
 *     a > 0 for rsqrt_approx); other lanes are unspecified. Each tier uses
 *     a hardware estimate refined by Newton steps (x86 below AVX-512
 *     estimates in float, from the significand) or the exact value.
 *   - floor, ceil, trunc and rint (half to even) are exact for every
 *     double, keep the sign of zero and pass NaN payloads through.
 *   - min(a, b) is a < b ? a : b and max(a, b) is a > b ? a : b, as on x86:
 *     when either operand is NaN the result is b.
 *   - andnot(a, b) is (~a) & b, as on x86.
 *   - blend(a, b, m) takes b in the lanes where m is set and a elsewhere.
 *   - cmp_ne is true when either operand is NaN; the other comparisons are
 *     false then. is_na tests for R's NA_real_ (a NaN whose low 32 bits are
 *     1954), is_nan for any NaN including NA. The integer is_na tests for
 *     INT32_MIN or INT64_MIN.
 *   - rsimd_vi32_inc(acc, m) and rsimd_vi64_inc(acc, m) add 1 to acc in the
 *     lanes where m is set (wrapping), for counting a mask without a
 *     horizontal step per vector.
 *   - Integer add, sub and mul wrap (two's complement). rsimd_vi32_mulhi
 *     is the high 32 bits of the signed 64-bit product.
 *   - loadu_i32 and storeu_i32 convert between int32 elements in memory and
 *     64-bit lanes, so they move RSIMD_LANES_64 elements (the _p forms take
 *     a 64-bit predicate). rsimd_vi64 sign-extends on load and keeps the low
 *     32 bits on store; rsimd_vf64 converts exactly on load and truncates
 *     toward zero on store, where every lane must be within the int32 range
 *     (the result for other lanes differs between tiers).
 *     rsimd_vi32_reduce_add returns the exact sum as int64_t;
 *     rsimd_vi64_reduce_add wraps.
 *   - rsimd_vi32_sll, _srl and _sra shift every lane left, right
 *     logically and right arithmetically by the same count k, a run-time
 *     value in [0, 31].
 *   - rsimd_vi32_loadu_u8 reads RSIMD_LANES_32 bytes into the 32-bit lanes,
 *     zero-extended; rsimd_vi32_storeu_u8 writes the low 8 bits of every
 *     lane as bytes (the _p forms take a 32-bit predicate).
 *   - reduce_add of f64 lanes uses a tier-specific order; reduce_min and
 *     reduce_max are unspecified when a lane is NaN.
 *   - Mask any/all/count look at every lane, so the lanes a predicated load
 *     filled with its fill value count too: pick a neutral fill value.
 *   - uzp_even(a, b) and uzp_odd(a, b) are the even and the odd lanes of
 *     the concatenation a:b (on the none tier a and b), which deinterleave
 *     the real and imaginary parts of complex numbers loaded as doubles.
 *   - Prefix scans (sse2, avx2, avx512, neon only; not none or the SVE
 *     tiers): rsimd_vf64_shift_up(v, k, fill) moves lane j to lane j + k
 *     for a power of two k < RSIMD_WIDTH_F64 and takes lanes 0 .. k - 1
 *     from fill; rsimd_vf64_bcast_last(v) copies the last lane to every
 *     lane; rsimd_vf64_first(v) is lane 0.
 *
 * Loops
 * -----
 * A predicate selects the active lanes: rsimd_p64_while(i, n) is the lanes
 * i, i + 1, ... that are < n, and rsimd_p64_true() is all lanes.
 * loadu_p(pg, p, fill) reads only the active lanes and sets the others to
 * fill; storeu_p(pg, p, v) writes only the active lanes. A kernel processes
 * full vectors and then the tail with the same code:
 *
 *   ptrdiff_t i = 0;
 *   for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
 *     rsimd_vf64_storeu(out + i, rsimd_vf64_add(rsimd_vf64_loadu(x + i), k));
 *   }
 *   if (i < n) {
 *     rsimd_p64 pg = rsimd_p64_while(i, n);
 *     rsimd_vf64 v = rsimd_vf64_loadu_p(pg, x + i, 0.0);
 *     rsimd_vf64_storeu_p(pg, out + i, rsimd_vf64_add(v, k));
 *   }
 *
 * The SVE tiers implement predicates with svwhilelt, AVX-512 with mask
 * registers, the other fixed-width tiers with a lane count and a small stack
 * buffer, and the none tier with a single lane.
 *
 * Ops that a tier has no instruction for (64-bit integer multiply on SSE2,
 * NEON and AVX2, for example) use SIMDe's portable code or a lane loop: they
 * are correct everywhere but may be slow. Kernels needing anything else use
 * #if RSIMD_TIER_IS(...) blocks, always with a generic fallback.
 *
 * On 32-bit ARM the neon tier has no f64 vector layer: RSIMD_NO_F64_SIMD is
 * defined and f64 kernels are left to the none tier.
 */

#ifndef RSIMD_KERNELS_COMMON_INC_H
#define RSIMD_KERNELS_COMMON_INC_H

#ifndef RSIMD_TIER
#error "define RSIMD_TIER before including kernels/common.inc.h"
#endif

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "tiers.h"

/* Tier identification and name mangling. */
#define RSIMD_PASTE_(a, b) a##b
#define RSIMD_PASTE(a, b) RSIMD_PASTE_(a, b)
#define RSIMD_CAT_(a, b) a##_##b
#define RSIMD_CAT(a, b) RSIMD_CAT_(a, b)
#define RSIMD_STR_(x) #x
#define RSIMD_STR(x) RSIMD_STR_(x)

#define RSIMD_TIERNUM_none 1
#define RSIMD_TIERNUM_sse2 2
#define RSIMD_TIERNUM_avx2 3
#define RSIMD_TIERNUM_avx512 4
#define RSIMD_TIERNUM_neon 5
#define RSIMD_TIERNUM_sve 6
#define RSIMD_TIERNUM_sve2 7

#if !RSIMD_PASTE(RSIMD_TIERNUM_, RSIMD_TIER)
#error "RSIMD_TIER is not a tier with kernels"
#endif

#define RSIMD_TIER_IS(t) \
  (RSIMD_PASTE(RSIMD_TIERNUM_, t) == RSIMD_PASTE(RSIMD_TIERNUM_, RSIMD_TIER))
#define RSIMD_TIER_STRING RSIMD_STR(RSIMD_TIER)
#define RSIMD_KERNEL(name) RSIMD_CAT(RSIMD_CAT(rsimd, name), RSIMD_TIER)

#define RSIMD_INLINE static inline

/* Fail loudly if the TU was compiled without its tier's flags, instead of
   letting SIMDe silently emulate the instructions. */
#if RSIMD_TIER_IS(sse2) && !(defined(__SSE2__) || defined(_M_X64))
#error "tier_sse2.c compiled without SSE2 flags"
#endif
#if RSIMD_TIER_IS(avx2) && !(defined(__AVX2__) && defined(__FMA__))
#error "tier_avx2.c compiled without AVX2 and FMA flags"
#endif
#if RSIMD_TIER_IS(avx512) && !(defined(__AVX512F__) && defined(__AVX512BW__) && \
                               defined(__AVX512DQ__) && defined(__AVX512VL__))
#error "tier_avx512.c compiled without AVX-512 F, BW, DQ and VL flags"
#endif
#if RSIMD_TIER_IS(neon) && !defined(__ARM_NEON)
#error "tier_neon.c compiled without NEON flags"
#endif
#if RSIMD_TIER_IS(sve) && !defined(__ARM_FEATURE_SVE)
#error "tier_sve.c compiled without SVE flags"
#endif
#if RSIMD_TIER_IS(sve2) && !defined(__ARM_FEATURE_SVE2)
#error "tier_sve2.c compiled without SVE2 flags"
#endif

#if RSIMD_TIER_IS(neon) && !defined(__aarch64__) && !defined(_M_ARM64)
#define RSIMD_NO_F64_SIMD 1
#endif

/* Population count of a lane bitmask, for the mask count operations. */
RSIMD_INLINE int rsimd_popcount32(uint32_t x) {
  x = x - ((x >> 1) & 0x55555555u);
  x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
  x = (x + (x >> 4)) & 0x0F0F0F0Fu;
  return (int) ((x * 0x01010101u) >> 24);
}

/* R's NA_real_ is a NaN whose low 32 bits are 1954. */
#define RSIMD_NA_LOW_WORD 1954

#if RSIMD_TIER_IS(none)
#include "vec_none.inc.h"
#elif RSIMD_TIER_IS(sse2) || RSIMD_TIER_IS(neon)
#include "vec_simde128.inc.h"
#elif RSIMD_TIER_IS(avx2)
#include "vec_simde256.inc.h"
#elif RSIMD_TIER_IS(avx512)
#include "vec_simde512.inc.h"
#elif RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
#include "vec_sve.inc.h"
#endif

/* Upper-case aliases. */
#define RSIMD_PRED64 rsimd_p64_while
#define RSIMD_PRED32 rsimd_p32_while
#define RSIMD_P64_WHILE rsimd_p64_while
#define RSIMD_P64_TRUE rsimd_p64_true
#define RSIMD_P64_COUNT rsimd_p64_count
#define RSIMD_P32_WHILE rsimd_p32_while
#define RSIMD_P32_TRUE rsimd_p32_true
#define RSIMD_P32_COUNT rsimd_p32_count

#ifndef RSIMD_NO_F64_SIMD
#define RSIMD_VF64_LOADU rsimd_vf64_loadu
#define RSIMD_VF64_STOREU rsimd_vf64_storeu
#define RSIMD_VF64_LOADU_P rsimd_vf64_loadu_p
#define RSIMD_VF64_STOREU_P rsimd_vf64_storeu_p
#define RSIMD_VF64_SET1 rsimd_vf64_set1
#define RSIMD_VF64_ZERO rsimd_vf64_zero
#define RSIMD_VF64_ADD rsimd_vf64_add
#define RSIMD_VF64_SUB rsimd_vf64_sub
#define RSIMD_VF64_MUL rsimd_vf64_mul
#define RSIMD_VF64_DIV rsimd_vf64_div
#define RSIMD_VF64_FMA rsimd_vf64_fma
#define RSIMD_VF64_MIN rsimd_vf64_min
#define RSIMD_VF64_MAX rsimd_vf64_max
#define RSIMD_VF64_ABS rsimd_vf64_abs
#define RSIMD_VF64_NEG rsimd_vf64_neg
#define RSIMD_VF64_SQRT rsimd_vf64_sqrt
#define RSIMD_VF64_FLOOR rsimd_vf64_floor
#define RSIMD_VF64_CEIL rsimd_vf64_ceil
#define RSIMD_VF64_TRUNC rsimd_vf64_trunc
#define RSIMD_VF64_RINT rsimd_vf64_rint
#define RSIMD_VF64_LOADU_I32 rsimd_vf64_loadu_i32
#define RSIMD_VF64_LOADU_I32_P rsimd_vf64_loadu_i32_p
#define RSIMD_VF64_STOREU_I32 rsimd_vf64_storeu_i32
#define RSIMD_VF64_STOREU_I32_P rsimd_vf64_storeu_i32_p
#define RSIMD_VF64_AND rsimd_vf64_and
#define RSIMD_VF64_OR rsimd_vf64_or
#define RSIMD_VF64_XOR rsimd_vf64_xor
#define RSIMD_VF64_ANDNOT rsimd_vf64_andnot
#define RSIMD_VF64_CMP_EQ rsimd_vf64_cmp_eq
#define RSIMD_VF64_CMP_NE rsimd_vf64_cmp_ne
#define RSIMD_VF64_CMP_LT rsimd_vf64_cmp_lt
#define RSIMD_VF64_CMP_LE rsimd_vf64_cmp_le
#define RSIMD_VF64_CMP_GT rsimd_vf64_cmp_gt
#define RSIMD_VF64_CMP_GE rsimd_vf64_cmp_ge
#define RSIMD_VF64_IS_NAN rsimd_vf64_is_nan
#define RSIMD_VF64_IS_NA rsimd_vf64_is_na
#define RSIMD_VF64_BLEND rsimd_vf64_blend
#define RSIMD_VF64_REDUCE_ADD rsimd_vf64_reduce_add
#define RSIMD_VF64_REDUCE_MIN rsimd_vf64_reduce_min
#define RSIMD_VF64_REDUCE_MAX rsimd_vf64_reduce_max
#define RSIMD_VF64_AS_VI64 rsimd_vf64_as_vi64
#define RSIMD_VF64_UZP_EVEN rsimd_vf64_uzp_even
#define RSIMD_VF64_UZP_ODD rsimd_vf64_uzp_odd
#if !RSIMD_TIER_IS(none) && !RSIMD_TIER_IS(sve) && !RSIMD_TIER_IS(sve2)
#define RSIMD_VF64_SHIFT_UP rsimd_vf64_shift_up
#define RSIMD_VF64_BCAST_LAST rsimd_vf64_bcast_last
#define RSIMD_VF64_FIRST rsimd_vf64_first
#endif
#define RSIMD_VI64_AS_VF64 rsimd_vi64_as_vf64
#define RSIMD_MF64_AND rsimd_mf64_and
#define RSIMD_MF64_OR rsimd_mf64_or
#define RSIMD_MF64_NOT rsimd_mf64_not
#define RSIMD_MF64_ANDNOT rsimd_mf64_andnot
#define RSIMD_MF64_ANY rsimd_mf64_any
#define RSIMD_MF64_ALL rsimd_mf64_all
#define RSIMD_MF64_COUNT rsimd_mf64_count
#define RSIMD_MF64_TO_MI64 rsimd_mf64_to_mi64
#define RSIMD_MI64_TO_MF64 rsimd_mi64_to_mf64
#endif

#define RSIMD_VI32_LOADU rsimd_vi32_loadu
#define RSIMD_VI32_STOREU rsimd_vi32_storeu
#define RSIMD_VI32_LOADU_P rsimd_vi32_loadu_p
#define RSIMD_VI32_STOREU_P rsimd_vi32_storeu_p
#define RSIMD_VI32_SET1 rsimd_vi32_set1
#define RSIMD_VI32_ZERO rsimd_vi32_zero
#define RSIMD_VI32_ADD rsimd_vi32_add
#define RSIMD_VI32_SUB rsimd_vi32_sub
#define RSIMD_VI32_MUL rsimd_vi32_mul
#define RSIMD_VI32_MULHI rsimd_vi32_mulhi
#define RSIMD_VI32_SLL rsimd_vi32_sll
#define RSIMD_VI32_SRL rsimd_vi32_srl
#define RSIMD_VI32_SRA rsimd_vi32_sra
#define RSIMD_VI32_LOADU_U8 rsimd_vi32_loadu_u8
#define RSIMD_VI32_LOADU_U8_P rsimd_vi32_loadu_u8_p
#define RSIMD_VI32_STOREU_U8 rsimd_vi32_storeu_u8
#define RSIMD_VI32_STOREU_U8_P rsimd_vi32_storeu_u8_p
#define RSIMD_VI32_MIN rsimd_vi32_min
#define RSIMD_VI32_MAX rsimd_vi32_max
#define RSIMD_VI32_AND rsimd_vi32_and
#define RSIMD_VI32_OR rsimd_vi32_or
#define RSIMD_VI32_XOR rsimd_vi32_xor
#define RSIMD_VI32_ANDNOT rsimd_vi32_andnot
#define RSIMD_VI32_CMP_EQ rsimd_vi32_cmp_eq
#define RSIMD_VI32_CMP_GT rsimd_vi32_cmp_gt
#define RSIMD_VI32_CMP_LT rsimd_vi32_cmp_lt
#define RSIMD_VI32_BLEND rsimd_vi32_blend
#define RSIMD_VI32_IS_NA rsimd_vi32_is_na
#define RSIMD_VI32_INC rsimd_vi32_inc
#define RSIMD_VI32_REDUCE_ADD rsimd_vi32_reduce_add
#define RSIMD_VI32_REDUCE_MIN rsimd_vi32_reduce_min
#define RSIMD_VI32_REDUCE_MAX rsimd_vi32_reduce_max
#define RSIMD_MI32_AND rsimd_mi32_and
#define RSIMD_MI32_OR rsimd_mi32_or
#define RSIMD_MI32_NOT rsimd_mi32_not
#define RSIMD_MI32_ANDNOT rsimd_mi32_andnot
#define RSIMD_MI32_ANY rsimd_mi32_any
#define RSIMD_MI32_ALL rsimd_mi32_all
#define RSIMD_MI32_COUNT rsimd_mi32_count

#define RSIMD_VI64_LOADU rsimd_vi64_loadu
#define RSIMD_VI64_STOREU rsimd_vi64_storeu
#define RSIMD_VI64_LOADU_P rsimd_vi64_loadu_p
#define RSIMD_VI64_STOREU_P rsimd_vi64_storeu_p
#define RSIMD_VI64_SET1 rsimd_vi64_set1
#define RSIMD_VI64_ZERO rsimd_vi64_zero
#define RSIMD_VI64_ADD rsimd_vi64_add
#define RSIMD_VI64_SUB rsimd_vi64_sub
#define RSIMD_VI64_MUL rsimd_vi64_mul
#define RSIMD_VI64_LOADU_I32 rsimd_vi64_loadu_i32
#define RSIMD_VI64_LOADU_I32_P rsimd_vi64_loadu_i32_p
#define RSIMD_VI64_STOREU_I32 rsimd_vi64_storeu_i32
#define RSIMD_VI64_STOREU_I32_P rsimd_vi64_storeu_i32_p
#define RSIMD_VI64_MIN rsimd_vi64_min
#define RSIMD_VI64_MAX rsimd_vi64_max
#define RSIMD_VI64_AND rsimd_vi64_and
#define RSIMD_VI64_OR rsimd_vi64_or
#define RSIMD_VI64_XOR rsimd_vi64_xor
#define RSIMD_VI64_ANDNOT rsimd_vi64_andnot
#define RSIMD_VI64_CMP_EQ rsimd_vi64_cmp_eq
#define RSIMD_VI64_CMP_GT rsimd_vi64_cmp_gt
#define RSIMD_VI64_CMP_LT rsimd_vi64_cmp_lt
#define RSIMD_VI64_BLEND rsimd_vi64_blend
#define RSIMD_VI64_IS_NA rsimd_vi64_is_na
#define RSIMD_VI64_INC rsimd_vi64_inc
#define RSIMD_VI64_REDUCE_ADD rsimd_vi64_reduce_add
#define RSIMD_VI64_REDUCE_MIN rsimd_vi64_reduce_min
#define RSIMD_VI64_REDUCE_MAX rsimd_vi64_reduce_max
#define RSIMD_MI64_AND rsimd_mi64_and
#define RSIMD_MI64_OR rsimd_mi64_or
#define RSIMD_MI64_NOT rsimd_mi64_not
#define RSIMD_MI64_ANDNOT rsimd_mi64_andnot
#define RSIMD_MI64_ANY rsimd_mi64_any
#define RSIMD_MI64_ALL rsimd_mi64_all
#define RSIMD_MI64_COUNT rsimd_mi64_count

/* 1 on the tiers with a native fused multiply-add, which the fast-mode
   folds use for product terms (na.h). SSE2 has none (its rsimd_vf64_fma
   calls libm's fma per lane), and the none tier adds rounded products. */
#if RSIMD_TIER_IS(avx2) || RSIMD_TIER_IS(avx512) || RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2) || \
  (RSIMD_TIER_IS(neon) && !defined(RSIMD_NO_F64_SIMD))
#define RSIMD_NATIVE_FMA 1
#else
#define RSIMD_NATIVE_FMA 0
#endif

/* NA, precision and overflow helpers (scalar and vector forms). */
#include "na.h"

/* Vectorised elementary functions (SLEEF)
 * ---------------------------------------
 * configure compiles a tier with -DRSIMD_HAVE_SLEEF_<TIER>=1 when the
 * tier's SLEEF inline header (src/vendor/sleef, named in tools/tiers.txt;
 * sve2 uses the sve one) compiles with the tier's flags. RSIMD_HAVE_SLEEF
 * is then 1 and these wrappers of rsimd_vf64 exist, each also under its
 * upper-case name (rsimd_sleef_exp is RSIMD_SLEEF_EXP):
 *
 *   rsimd_sleef_<f>(v)   exp exp2 exp10 expm1 log log2 log10 log1p cbrt
 *                        sin cos tan asin acos atan sinh cosh tanh asinh
 *                        acosh atanh, and sinpi cospi
 *   rsimd_sleef_<f>(a, b)  pow atan2, and hypot; nextafter and remainder,
 *                        which are exact
 *   rsimd_sleef_sincos(v, &s, &c), rsimd_sleef_sincospi(v, &s, &c)
 *                        store the sine and cosine through pointers (SVE
 *                        vectors cannot be struct members)
 *   rsimd_sleef_<f>_fast  the same functions in fast mode (option
 *                        rsimd.math_accuracy = "fast")
 *
 * Accuracy policy: every function uses SLEEF's 1.0-ULP (_u10) variant,
 * except sinpi, cospi, sincospi and hypot, which use the 0.5-ULP (_u05)
 * one. The _fast wrappers use the 3.5-ULP (_u35) variant where SLEEF has
 * one and it is faster: sin cos tan asin acos atan atan2 log log2 cbrt
 * sinh cosh tanh hypot sincos sincospi, and sinpi and cospi as one half of
 * sincospi_u35. The others (exp expm1 log10 log1p pow asinh acosh atanh,
 * which have no _u35, and exp2 and exp10, whose _u35 is no faster) are
 * the accurate wrappers under the _fast name. SLEEF's functions
 * handle the special values (NaN, infinities, signed zeros) as C99 does but
 * do not carry R's NA payload through, so kernels blend NA back in. The none
 * tier never uses SLEEF: its math kernels call C99 libm, as base R does, and
 * are the reference the other tiers are tested against. A tier without
 * RSIMD_HAVE_SLEEF leaves its math slots empty, so the dispatcher takes them
 * from the next lower tier that has them, ultimately none.
 *
 * SLEEF's functions and helpers are static, named with an _<isa>_sleef
 * suffix and only compiled when used; a tier includes one header. The
 * headers need -ffp-contract=off, which configure gives every tier. */
#if (RSIMD_TIER_IS(sse2) && defined(RSIMD_HAVE_SLEEF_SSE2)) ||           \
  (RSIMD_TIER_IS(avx2) && defined(RSIMD_HAVE_SLEEF_AVX2)) ||             \
  (RSIMD_TIER_IS(avx512) && defined(RSIMD_HAVE_SLEEF_AVX512)) ||         \
  (RSIMD_TIER_IS(neon) && defined(RSIMD_HAVE_SLEEF_NEON) &&              \
   !defined(RSIMD_NO_F64_SIMD)) ||                                       \
  (RSIMD_TIER_IS(sve) && defined(RSIMD_HAVE_SLEEF_SVE)) ||               \
  (RSIMD_TIER_IS(sve2) && defined(RSIMD_HAVE_SLEEF_SVE2))
#define RSIMD_HAVE_SLEEF 1
#endif

#ifdef RSIMD_HAVE_SLEEF
/* RSIMD_SLEEF_FN(exp, u10) is the tier's Sleef_expd<W>_u10<isa>;
   RSIMD_SLEEF_IN and _OUT convert between rsimd_vf64 and SLEEF's vector
   type; RSIMD_SLEEF_PAIR is the type of a two-result value, read with
   RSIMD_SLEEF_FIRST and _SECOND. */
#if RSIMD_TIER_IS(sse2)
#include <emmintrin.h>
#include "sleefinline_sse2.h"
#define RSIMD_SLEEF_FN(f, acc) Sleef_##f##d2_##acc##sse2
#define RSIMD_SLEEF_PAIR vdouble2_sse2_sleef
#elif RSIMD_TIER_IS(avx2)
#include <immintrin.h>
#include "sleefinline_avx2.h"
#define RSIMD_SLEEF_FN(f, acc) Sleef_##f##d4_##acc##avx2
#define RSIMD_SLEEF_PAIR vdouble2_avx2_sleef
#elif RSIMD_TIER_IS(avx512)
#include <immintrin.h>
#include "sleefinline_avx512f.h"
#define RSIMD_SLEEF_FN(f, acc) Sleef_##f##d8_##acc##avx512f
#define RSIMD_SLEEF_PAIR vdouble2_avx512f_sleef
#elif RSIMD_TIER_IS(neon)
#include <arm_neon.h>
#include "sleefinline_advsimd.h"
#define RSIMD_SLEEF_FN(f, acc) Sleef_##f##d2_##acc##advsimd
#define RSIMD_SLEEF_PAIR vdouble2_advsimd_sleef
#else /* sve, sve2 */
#include <arm_sve.h>
#include "sleefinline_sve.h"
#define RSIMD_SLEEF_FN(f, acc) Sleef_##f##dx_##acc##sve
#define RSIMD_SLEEF_PAIR vdouble2_sve_sleef
#endif

#if RSIMD_TIER_IS(neon)
#define RSIMD_SLEEF_IN(v) simde__m128d_to_neon_f64(v)
#define RSIMD_SLEEF_OUT(v) simde__m128d_from_neon_f64(v)
#else
#define RSIMD_SLEEF_IN(v) (v)
#define RSIMD_SLEEF_OUT(v) (v)
#endif
#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
#define RSIMD_SLEEF_FIRST(r) svget2_f64((r), 0)
#define RSIMD_SLEEF_SECOND(r) svget2_f64((r), 1)
#else
#define RSIMD_SLEEF_FIRST(r) ((r).x)
#define RSIMD_SLEEF_SECOND(r) ((r).y)
#endif

#define RSIMD_SLEEF_DEF1(f, acc)                                                           \
  RSIMD_INLINE rsimd_vf64 rsimd_sleef_##f(rsimd_vf64 a) {                                  \
    return RSIMD_SLEEF_OUT(RSIMD_SLEEF_FN(f, acc)(RSIMD_SLEEF_IN(a)));                     \
  }

#define RSIMD_SLEEF_DEF2(f, acc)                                                           \
  RSIMD_INLINE rsimd_vf64 rsimd_sleef_##f(rsimd_vf64 a, rsimd_vf64 b) {                    \
    return RSIMD_SLEEF_OUT(RSIMD_SLEEF_FN(f, acc)(RSIMD_SLEEF_IN(a), RSIMD_SLEEF_IN(b)));  \
  }

#define RSIMD_SLEEF_DEFSC(f, acc)                                                          \
  RSIMD_INLINE void rsimd_sleef_##f(rsimd_vf64 a, rsimd_vf64 *s, rsimd_vf64 *c) {          \
    RSIMD_SLEEF_PAIR r = RSIMD_SLEEF_FN(f, acc)(RSIMD_SLEEF_IN(a));                        \
    *s = RSIMD_SLEEF_OUT(RSIMD_SLEEF_FIRST(r));                                            \
    *c = RSIMD_SLEEF_OUT(RSIMD_SLEEF_SECOND(r));                                           \
  }

RSIMD_SLEEF_DEF1(exp, u10)
RSIMD_SLEEF_DEF1(exp2, u10)
RSIMD_SLEEF_DEF1(exp10, u10)
RSIMD_SLEEF_DEF1(expm1, u10)
RSIMD_SLEEF_DEF1(log, u10)
RSIMD_SLEEF_DEF1(log2, u10)
RSIMD_SLEEF_DEF1(log10, u10)
RSIMD_SLEEF_DEF1(log1p, u10)
RSIMD_SLEEF_DEF1(cbrt, u10)
RSIMD_SLEEF_DEF1(sin, u10)
RSIMD_SLEEF_DEF1(cos, u10)
RSIMD_SLEEF_DEF1(tan, u10)
RSIMD_SLEEF_DEF1(asin, u10)
RSIMD_SLEEF_DEF1(acos, u10)
RSIMD_SLEEF_DEF1(atan, u10)
RSIMD_SLEEF_DEF1(sinh, u10)
RSIMD_SLEEF_DEF1(cosh, u10)
RSIMD_SLEEF_DEF1(tanh, u10)
RSIMD_SLEEF_DEF1(asinh, u10)
RSIMD_SLEEF_DEF1(acosh, u10)
RSIMD_SLEEF_DEF1(atanh, u10)
RSIMD_SLEEF_DEF1(sinpi, u05)
RSIMD_SLEEF_DEF1(cospi, u05)
RSIMD_SLEEF_DEF2(pow, u10)
RSIMD_SLEEF_DEF2(atan2, u10)
RSIMD_SLEEF_DEF2(hypot, u05)
RSIMD_SLEEF_DEF2(nextafter, )
RSIMD_SLEEF_DEF2(remainder, )
RSIMD_SLEEF_DEFSC(sincos, u10)
RSIMD_SLEEF_DEFSC(sincospi, u05)

/* Fast mode. RSIMD_SLEEF_DEF*_FAST(f) defines rsimd_sleef_f_fast with
   SLEEF's _u35 variant. */
#define RSIMD_SLEEF_DEF1_FAST(f)                                                           \
  RSIMD_INLINE rsimd_vf64 rsimd_sleef_##f##_fast(rsimd_vf64 a) {                           \
    return RSIMD_SLEEF_OUT(RSIMD_SLEEF_FN(f, u35)(RSIMD_SLEEF_IN(a)));                     \
  }
#define RSIMD_SLEEF_DEF2_FAST(f)                                                           \
  RSIMD_INLINE rsimd_vf64 rsimd_sleef_##f##_fast(rsimd_vf64 a, rsimd_vf64 b) {             \
    return RSIMD_SLEEF_OUT(RSIMD_SLEEF_FN(f, u35)(RSIMD_SLEEF_IN(a), RSIMD_SLEEF_IN(b)));  \
  }
#define RSIMD_SLEEF_DEFSC_FAST(f)                                                          \
  RSIMD_INLINE void rsimd_sleef_##f##_fast(rsimd_vf64 a, rsimd_vf64 *s, rsimd_vf64 *c) {   \
    RSIMD_SLEEF_PAIR r = RSIMD_SLEEF_FN(f, u35)(RSIMD_SLEEF_IN(a));                        \
    *s = RSIMD_SLEEF_OUT(RSIMD_SLEEF_FIRST(r));                                            \
    *c = RSIMD_SLEEF_OUT(RSIMD_SLEEF_SECOND(r));                                           \
  }

RSIMD_SLEEF_DEF1_FAST(log)
RSIMD_SLEEF_DEF1_FAST(log2)
RSIMD_SLEEF_DEF1_FAST(cbrt)
RSIMD_SLEEF_DEF1_FAST(sin)
RSIMD_SLEEF_DEF1_FAST(cos)
RSIMD_SLEEF_DEF1_FAST(tan)
RSIMD_SLEEF_DEF1_FAST(asin)
RSIMD_SLEEF_DEF1_FAST(acos)
RSIMD_SLEEF_DEF1_FAST(atan)
RSIMD_SLEEF_DEF1_FAST(sinh)
RSIMD_SLEEF_DEF1_FAST(cosh)
RSIMD_SLEEF_DEF1_FAST(tanh)
RSIMD_SLEEF_DEF2_FAST(atan2)
RSIMD_SLEEF_DEF2_FAST(hypot)
RSIMD_SLEEF_DEFSC_FAST(sincos)
RSIMD_SLEEF_DEFSC_FAST(sincospi)

/* SLEEF has no sinpi_u35 or cospi_u35. */
RSIMD_INLINE rsimd_vf64 rsimd_sleef_sinpi_fast(rsimd_vf64 a) {
  return RSIMD_SLEEF_OUT(RSIMD_SLEEF_FIRST(RSIMD_SLEEF_FN(sincospi, u35)(RSIMD_SLEEF_IN(a))));
}
RSIMD_INLINE rsimd_vf64 rsimd_sleef_cospi_fast(rsimd_vf64 a) {
  return RSIMD_SLEEF_OUT(RSIMD_SLEEF_SECOND(RSIMD_SLEEF_FN(sincospi, u35)(RSIMD_SLEEF_IN(a))));
}

#define rsimd_sleef_exp_fast rsimd_sleef_exp
#define rsimd_sleef_exp2_fast rsimd_sleef_exp2
#define rsimd_sleef_exp10_fast rsimd_sleef_exp10
#define rsimd_sleef_expm1_fast rsimd_sleef_expm1
#define rsimd_sleef_log10_fast rsimd_sleef_log10
#define rsimd_sleef_log1p_fast rsimd_sleef_log1p
#define rsimd_sleef_asinh_fast rsimd_sleef_asinh
#define rsimd_sleef_acosh_fast rsimd_sleef_acosh
#define rsimd_sleef_atanh_fast rsimd_sleef_atanh
#define rsimd_sleef_pow_fast rsimd_sleef_pow

#define RSIMD_SLEEF_EXP rsimd_sleef_exp
#define RSIMD_SLEEF_EXP2 rsimd_sleef_exp2
#define RSIMD_SLEEF_EXP10 rsimd_sleef_exp10
#define RSIMD_SLEEF_EXPM1 rsimd_sleef_expm1
#define RSIMD_SLEEF_LOG rsimd_sleef_log
#define RSIMD_SLEEF_LOG2 rsimd_sleef_log2
#define RSIMD_SLEEF_LOG10 rsimd_sleef_log10
#define RSIMD_SLEEF_LOG1P rsimd_sleef_log1p
#define RSIMD_SLEEF_CBRT rsimd_sleef_cbrt
#define RSIMD_SLEEF_SIN rsimd_sleef_sin
#define RSIMD_SLEEF_COS rsimd_sleef_cos
#define RSIMD_SLEEF_TAN rsimd_sleef_tan
#define RSIMD_SLEEF_ASIN rsimd_sleef_asin
#define RSIMD_SLEEF_ACOS rsimd_sleef_acos
#define RSIMD_SLEEF_ATAN rsimd_sleef_atan
#define RSIMD_SLEEF_SINH rsimd_sleef_sinh
#define RSIMD_SLEEF_COSH rsimd_sleef_cosh
#define RSIMD_SLEEF_TANH rsimd_sleef_tanh
#define RSIMD_SLEEF_ASINH rsimd_sleef_asinh
#define RSIMD_SLEEF_ACOSH rsimd_sleef_acosh
#define RSIMD_SLEEF_ATANH rsimd_sleef_atanh
#define RSIMD_SLEEF_SINPI rsimd_sleef_sinpi
#define RSIMD_SLEEF_COSPI rsimd_sleef_cospi
#define RSIMD_SLEEF_POW rsimd_sleef_pow
#define RSIMD_SLEEF_ATAN2 rsimd_sleef_atan2
#define RSIMD_SLEEF_HYPOT rsimd_sleef_hypot
#define RSIMD_SLEEF_NEXTAFTER rsimd_sleef_nextafter
#define RSIMD_SLEEF_REMAINDER rsimd_sleef_remainder
#define RSIMD_SLEEF_SINCOS rsimd_sleef_sincos
#define RSIMD_SLEEF_SINCOSPI rsimd_sleef_sincospi
#define RSIMD_SLEEF_EXP_FAST rsimd_sleef_exp_fast
#define RSIMD_SLEEF_EXP2_FAST rsimd_sleef_exp2_fast
#define RSIMD_SLEEF_EXP10_FAST rsimd_sleef_exp10_fast
#define RSIMD_SLEEF_EXPM1_FAST rsimd_sleef_expm1_fast
#define RSIMD_SLEEF_LOG_FAST rsimd_sleef_log_fast
#define RSIMD_SLEEF_LOG2_FAST rsimd_sleef_log2_fast
#define RSIMD_SLEEF_LOG10_FAST rsimd_sleef_log10_fast
#define RSIMD_SLEEF_LOG1P_FAST rsimd_sleef_log1p_fast
#define RSIMD_SLEEF_CBRT_FAST rsimd_sleef_cbrt_fast
#define RSIMD_SLEEF_SIN_FAST rsimd_sleef_sin_fast
#define RSIMD_SLEEF_COS_FAST rsimd_sleef_cos_fast
#define RSIMD_SLEEF_TAN_FAST rsimd_sleef_tan_fast
#define RSIMD_SLEEF_ASIN_FAST rsimd_sleef_asin_fast
#define RSIMD_SLEEF_ACOS_FAST rsimd_sleef_acos_fast
#define RSIMD_SLEEF_ATAN_FAST rsimd_sleef_atan_fast
#define RSIMD_SLEEF_SINH_FAST rsimd_sleef_sinh_fast
#define RSIMD_SLEEF_COSH_FAST rsimd_sleef_cosh_fast
#define RSIMD_SLEEF_TANH_FAST rsimd_sleef_tanh_fast
#define RSIMD_SLEEF_ASINH_FAST rsimd_sleef_asinh_fast
#define RSIMD_SLEEF_ACOSH_FAST rsimd_sleef_acosh_fast
#define RSIMD_SLEEF_ATANH_FAST rsimd_sleef_atanh_fast
#define RSIMD_SLEEF_SINPI_FAST rsimd_sleef_sinpi_fast
#define RSIMD_SLEEF_COSPI_FAST rsimd_sleef_cospi_fast
#define RSIMD_SLEEF_POW_FAST rsimd_sleef_pow_fast
#define RSIMD_SLEEF_ATAN2_FAST rsimd_sleef_atan2_fast
#define RSIMD_SLEEF_HYPOT_FAST rsimd_sleef_hypot_fast
#define RSIMD_SLEEF_SINCOS_FAST rsimd_sleef_sincos_fast
#define RSIMD_SLEEF_SINCOSPI_FAST rsimd_sleef_sincospi_fast
#endif /* RSIMD_HAVE_SLEEF */

#endif /* RSIMD_KERNELS_COMMON_INC_H */
