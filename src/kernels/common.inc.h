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
 *   rsimd_vf64 also   div fma abs neg sqrt cmp_ne cmp_le cmp_ge is_nan
 *                     as_vi64
 *   rsimd_vi64 also   as_vf64
 *   mask types        and or not andnot any all count; rsimd_mf64_to_mi64
 *                     and rsimd_mi64_to_mf64 convert between the f64 and
 *                     i64 masks
 *   predicate types   while true count
 *
 * Semantics shared by all tiers:
 *   - Loads and stores are unaligned (R vectors are only 8-byte aligned).
 *   - fma(a, b, c) is a * b + c with a single rounding, on every tier.
 *   - min(a, b) is a < b ? a : b and max(a, b) is a > b ? a : b, as on x86:
 *     when either operand is NaN the result is b.
 *   - andnot(a, b) is (~a) & b, as on x86.
 *   - blend(a, b, m) takes b in the lanes where m is set and a elsewhere.
 *   - cmp_ne is true when either operand is NaN; the other comparisons are
 *     false then. is_na tests for R's NA_real_ (a NaN whose low 32 bits are
 *     1954), is_nan for any NaN including NA. The integer is_na tests for
 *     INT32_MIN or INT64_MIN.
 *   - Integer add, sub and mul wrap (two's complement).
 *     rsimd_vi32_reduce_add returns the exact sum as int64_t;
 *     rsimd_vi64_reduce_add wraps.
 *   - reduce_add of f64 lanes uses a tier-specific order; reduce_min and
 *     reduce_max are unspecified when a lane is NaN.
 *   - Mask any/all/count look at every lane, so the lanes a predicated load
 *     filled with its fill value count too: pick a neutral fill value.
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

/* Vectorised elementary functions (SLEEF) are included here per tier once
   they are bundled; RSIMD_HAVE_SLEEF_<TIER> then says whether the tier has
   them. */

#endif /* RSIMD_KERNELS_COMMON_INC_H */
