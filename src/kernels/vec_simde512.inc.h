/* 512-bit vector layer over SIMDe's AVX-512 F+BW+DQ+VL intrinsics, for the
   avx512 tier. Comparisons return mask registers, and predicates are mask
   registers too, so predicated loads and stores are native masked
   instructions (they never touch inactive lanes). Included by
   kernels/common.inc.h only. */

#include "x86/avx2.h"
#include "x86/fma.h"
#include "x86/avx512/types.h"
#include "x86/avx512/loadu.h"
#include "x86/avx512/storeu.h"
#include "x86/avx512/set1.h"
#include "x86/avx512/setzero.h"
#include "x86/avx512/add.h"
#include "x86/avx512/sub.h"
#include "x86/avx512/mul.h"
#include "x86/avx512/mullo.h"
#include "x86/avx512/div.h"
#include "x86/avx512/fmadd.h"
#include "x86/avx512/max.h"
#include "x86/avx512/min.h"
#include "x86/avx512/abs.h"
#include "x86/avx512/sqrt.h"
#include "x86/avx512/and.h"
#include "x86/avx512/andnot.h"
#include "x86/avx512/or.h"
#include "x86/avx512/xor.h"
#include "x86/avx512/cmp.h"
#include "x86/avx512/cmpeq.h"
#include "x86/avx512/cmpgt.h"
#include "x86/avx512/blend.h"
#include "x86/avx512/cast.h"
#include "x86/avx512/reduce.h"
#include "x86/avx512/srli.h"
#include "x86/avx512/cvt.h"
#include "x86/avx512/extract.h"
#include "x86/avx512/roundscale.h"
#include "x86/avx512/set.h"
#include "x86/avx512/permutexvar.h"
#include "x86/avx512/sll.h"
#include "x86/avx512/srl.h"
#include "x86/avx512/srai.h"
#include "x86/avx512/shuffle.h"

#define RSIMD_WIDTH_F64 8
#define RSIMD_WIDTH_I64 8
#define RSIMD_WIDTH_I32 16
#define RSIMD_LANES_64 ((ptrdiff_t) 8)
#define RSIMD_LANES_32 ((ptrdiff_t) 16)

typedef simde__m512d rsimd_vf64;
typedef simde__m512i rsimd_vi32;
typedef simde__m512i rsimd_vi64;
typedef simde__mmask8 rsimd_mf64; /* one bit per lane */
typedef simde__mmask16 rsimd_mi32;
typedef simde__mmask8 rsimd_mi64;
typedef simde__mmask8 rsimd_p64;
typedef simde__mmask16 rsimd_p32;

/* Predicates */
RSIMD_INLINE rsimd_p64 rsimd_p64_while(ptrdiff_t i, ptrdiff_t n) {
  ptrdiff_t k = n - i;
  if (k >= 8) return (rsimd_p64) 0xFF;
  return (rsimd_p64) (k > 0 ? (1u << k) - 1u : 0u);
}
RSIMD_INLINE rsimd_p64 rsimd_p64_true(void) { return (rsimd_p64) 0xFF; }
RSIMD_INLINE int rsimd_p64_count(rsimd_p64 pg) { return rsimd_popcount32(pg); }
RSIMD_INLINE rsimd_p32 rsimd_p32_while(ptrdiff_t i, ptrdiff_t n) {
  ptrdiff_t k = n - i;
  if (k >= 16) return (rsimd_p32) 0xFFFF;
  return (rsimd_p32) (k > 0 ? (1u << k) - 1u : 0u);
}
RSIMD_INLINE rsimd_p32 rsimd_p32_true(void) { return (rsimd_p32) 0xFFFF; }
RSIMD_INLINE int rsimd_p32_count(rsimd_p32 pg) { return rsimd_popcount32(pg); }

/* Masks */
#define RSIMD_K_MASK(m, all)                                                      \
  RSIMD_INLINE rsimd_##m rsimd_##m##_and(rsimd_##m a, rsimd_##m b) {              \
    return (rsimd_##m) (a & b);                                                  \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_or(rsimd_##m a, rsimd_##m b) {               \
    return (rsimd_##m) (a | b);                                                  \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_not(rsimd_##m a) { return (rsimd_##m) ~a; }  \
  RSIMD_INLINE rsimd_##m rsimd_##m##_andnot(rsimd_##m a, rsimd_##m b) {           \
    return (rsimd_##m) (~a & b);                                                 \
  }                                                                              \
  RSIMD_INLINE int rsimd_##m##_any(rsimd_##m a) { return a != 0; }                \
  RSIMD_INLINE int rsimd_##m##_all(rsimd_##m a) { return a == (all); }            \
  RSIMD_INLINE int rsimd_##m##_count(rsimd_##m a) { return rsimd_popcount32(a); }
RSIMD_K_MASK(mf64, 0xFF)
RSIMD_K_MASK(mi32, 0xFFFF)
RSIMD_K_MASK(mi64, 0xFF)
#undef RSIMD_K_MASK

RSIMD_INLINE rsimd_mf64 rsimd_mi64_to_mf64(rsimd_mi64 m) { return m; }
RSIMD_INLINE rsimd_mi64 rsimd_mf64_to_mi64(rsimd_mf64 m) { return m; }

/* Doubles */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu(const double *p) { return simde_mm512_loadu_pd(p); }
RSIMD_INLINE void rsimd_vf64_storeu(double *p, rsimd_vf64 v) { simde_mm512_storeu_pd(p, v); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_p(rsimd_p64 pg, const double *p, double fill) {
  return simde_mm512_mask_loadu_pd(simde_mm512_set1_pd(fill), pg, p);
}
RSIMD_INLINE void rsimd_vf64_storeu_p(rsimd_p64 pg, double *p, rsimd_vf64 v) {
  simde_mm512_mask_storeu_pd(p, pg, v);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_set1(double x) { return simde_mm512_set1_pd(x); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_zero(void) { return simde_mm512_setzero_pd(); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_add(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_add_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sub(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_sub_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_mul(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_mul_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_div(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_div_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_fma(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c) {
  return simde_mm512_fmadd_pd(a, b, c);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_min(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_min_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_max(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_max_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_abs(rsimd_vf64 a) { return simde_mm512_abs_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_neg(rsimd_vf64 a) {
  return simde_mm512_xor_pd(simde_mm512_set1_pd(-0.0), a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sqrt(rsimd_vf64 a) { return simde_mm512_sqrt_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_and(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_and_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_or(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_or_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_xor(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm512_xor_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_andnot(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_andnot_pd(a, b);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_eq(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_EQ_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ne(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_NEQ_UQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_lt(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_LT_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_le(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_LE_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_gt(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_GT_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ge(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm512_cmp_pd_mask(a, b, SIMDE_CMP_GE_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_nan(rsimd_vf64 a) {
  return simde_mm512_cmp_pd_mask(a, a, SIMDE_CMP_UNORD_Q);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_na(rsimd_vf64 a) {
  simde__m512i lo = simde_mm512_and_si512(simde_mm512_castpd_si512(a),
                                          simde_mm512_set1_epi64(0xFFFFFFFF));
  return (rsimd_mf64) (simde_mm512_cmpeq_epi64_mask(lo, simde_mm512_set1_epi64(RSIMD_NA_LOW_WORD)) &
                       rsimd_vf64_is_nan(a));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_blend(rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m) {
  return simde_mm512_mask_blend_pd(m, a, b);
}
RSIMD_INLINE double rsimd_vf64_reduce_add(rsimd_vf64 a) { return simde_mm512_reduce_add_pd(a); }
RSIMD_INLINE double rsimd_vf64_reduce_min(rsimd_vf64 a) { return simde_mm512_reduce_min_pd(a); }
RSIMD_INLINE double rsimd_vf64_reduce_max(rsimd_vf64 a) { return simde_mm512_reduce_max_pd(a); }
/* Lane moves for prefix scans (see the 128-bit layer); k is 1, 2 or 4. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_shift_up(rsimd_vf64 v, int k, rsimd_vf64 fill) {
  simde__m512i idx = simde_mm512_set_epi64((7 - k) & 7, (6 - k) & 7, (5 - k) & 7, (4 - k) & 7,
                                           (3 - k) & 7, (2 - k) & 7, (1 - k) & 7, (0 - k) & 7);
  return simde_mm512_mask_blend_pd((simde__mmask8) ((1u << k) - 1u),
                                   simde_mm512_permutexvar_pd(idx, v), fill);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_bcast_last(rsimd_vf64 v) {
  return simde_mm512_permutexvar_pd(simde_mm512_set1_epi64(7), v);
}
RSIMD_INLINE double rsimd_vf64_first(rsimd_vf64 v) {
  return simde_mm_cvtsd_f64(simde_mm512_castpd512_pd128(v));
}
RSIMD_INLINE rsimd_vi64 rsimd_vf64_as_vi64(rsimd_vf64 a) { return simde_mm512_castpd_si512(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vi64_as_vf64(rsimd_vi64 a) { return simde_mm512_castsi512_pd(a); }

#define RSIMD_S512_BITWISE(v)                                                     \
  RSIMD_INLINE rsimd_##v rsimd_##v##_and(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm512_and_si512(a, b);                                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_or(rsimd_##v a, rsimd_##v b) {               \
    return simde_mm512_or_si512(a, b);                                           \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_xor(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm512_xor_si512(a, b);                                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_andnot(rsimd_##v a, rsimd_##v b) {           \
    return simde_mm512_andnot_si512(a, b);                                       \
  }

/* 32-bit integers */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu(const int32_t *p) { return simde_mm512_loadu_si512(p); }
RSIMD_INLINE void rsimd_vi32_storeu(int32_t *p, rsimd_vi32 v) { simde_mm512_storeu_si512(p, v); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_p(rsimd_p32 pg, const int32_t *p, int32_t fill) {
  return simde_mm512_mask_loadu_epi32(simde_mm512_set1_epi32(fill), pg, p);
}
RSIMD_INLINE void rsimd_vi32_storeu_p(rsimd_p32 pg, int32_t *p, rsimd_vi32 v) {
  simde_mm512_mask_storeu_epi32(p, pg, v);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_set1(int32_t x) { return simde_mm512_set1_epi32(x); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_zero(void) { return simde_mm512_setzero_si512(); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_add(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm512_add_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sub(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm512_sub_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mul(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm512_mullo_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_min(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm512_min_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_max(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm512_max_epi32(a, b); }
RSIMD_S512_BITWISE(vi32)
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_eq(rsimd_vi32 a, rsimd_vi32 b) {
  return simde_mm512_cmpeq_epi32_mask(a, b);
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_gt(rsimd_vi32 a, rsimd_vi32 b) {
  return simde_mm512_cmpgt_epi32_mask(a, b);
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_lt(rsimd_vi32 a, rsimd_vi32 b) {
  return simde_mm512_cmpgt_epi32_mask(b, a);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_blend(rsimd_vi32 a, rsimd_vi32 b, rsimd_mi32 m) {
  return simde_mm512_mask_blend_epi32(m, a, b);
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_is_na(rsimd_vi32 a) {
  return simde_mm512_cmpeq_epi32_mask(a, simde_mm512_set1_epi32(INT32_MIN));
}
RSIMD_INLINE int64_t rsimd_vi32_reduce_add(rsimd_vi32 a) {
  int32_t buf[16];
  int64_t r = 0;
  int j;
  simde_mm512_storeu_si512(buf, a);
  for (j = 0; j < 16; j++) r += buf[j];
  return r;
}
RSIMD_INLINE int32_t rsimd_vi32_reduce_min(rsimd_vi32 a) { return simde_mm512_reduce_min_epi32(a); }
RSIMD_INLINE int32_t rsimd_vi32_reduce_max(rsimd_vi32 a) { return simde_mm512_reduce_max_epi32(a); }

/* 64-bit integers */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu(const int64_t *p) { return simde_mm512_loadu_si512(p); }
RSIMD_INLINE void rsimd_vi64_storeu(int64_t *p, rsimd_vi64 v) { simde_mm512_storeu_si512(p, v); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_p(rsimd_p64 pg, const int64_t *p, int64_t fill) {
  return simde_mm512_mask_loadu_epi64(simde_mm512_set1_epi64(fill), pg, p);
}
RSIMD_INLINE void rsimd_vi64_storeu_p(rsimd_p64 pg, int64_t *p, rsimd_vi64 v) {
  simde_mm512_mask_storeu_epi64(p, pg, v);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_set1(int64_t x) { return simde_mm512_set1_epi64(x); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_zero(void) { return simde_mm512_setzero_si512(); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_add(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm512_add_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sub(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm512_sub_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_mul(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm512_mullo_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_min(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm512_min_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_max(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm512_max_epi64(a, b); }
RSIMD_S512_BITWISE(vi64)
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_eq(rsimd_vi64 a, rsimd_vi64 b) {
  return simde_mm512_cmpeq_epi64_mask(a, b);
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_gt(rsimd_vi64 a, rsimd_vi64 b) {
  return simde_mm512_cmpgt_epi64_mask(a, b);
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_lt(rsimd_vi64 a, rsimd_vi64 b) {
  return simde_mm512_cmpgt_epi64_mask(b, a);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_blend(rsimd_vi64 a, rsimd_vi64 b, rsimd_mi64 m) {
  return simde_mm512_mask_blend_epi64(m, a, b);
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_is_na(rsimd_vi64 a) {
  return simde_mm512_cmpeq_epi64_mask(a, simde_mm512_set1_epi64(INT64_MIN));
}
RSIMD_INLINE int64_t rsimd_vi64_reduce_add(rsimd_vi64 a) {
  int64_t buf[8];
  uint64_t r = 0;
  int j;
  simde_mm512_storeu_si512(buf, a);
  for (j = 0; j < 8; j++) r += (uint64_t) buf[j];
  return (int64_t) r;
}
RSIMD_INLINE int64_t rsimd_vi64_reduce_min(rsimd_vi64 a) { return simde_mm512_reduce_min_epi64(a); }
RSIMD_INLINE int64_t rsimd_vi64_reduce_max(rsimd_vi64 a) { return simde_mm512_reduce_max_epi64(a); }

#undef RSIMD_S512_BITWISE

/* Width conversions and extras. The 64-bit-lane conversions move
   RSIMD_WIDTH_I64 (8) int32 elements. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mulhi(rsimd_vi32 a, rsimd_vi32 b) {
  /* See the 128-bit layer. */
  simde__m512i even = simde_mm512_mul_epi32(a, b);
  simde__m512i odd =
    simde_mm512_mul_epi32(simde_mm512_srli_epi64(a, 32), simde_mm512_srli_epi64(b, 32));
  return simde_mm512_mask_blend_epi32((simde__mmask16) 0xAAAA, simde_mm512_srli_epi64(even, 32),
                                      odd);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sll(rsimd_vi32 a, int k) {
  return simde_mm512_sll_epi32(a, simde_mm_cvtsi32_si128(k));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_srl(rsimd_vi32 a, int k) {
  return simde_mm512_srl_epi32(a, simde_mm_cvtsi32_si128(k));
}
/* The vendored SIMDe subset has no 512-bit arithmetic shift by a count
   register: shift the complement of negative lanes logically, which is
   the same, and complement back. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sra(rsimd_vi32 a, int k) {
  simde__m512i s = simde_mm512_srai_epi32(a, 31);
  return simde_mm512_xor_si512(rsimd_vi32_srl(simde_mm512_xor_si512(a, s), k), s);
}
/* Sixteen bytes zero-extended to 32-bit lanes, and back (low 8 bits:
   byte 0 of each lane gathered into the low 128 bits). */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_u8(const uint8_t *p) {
  return simde_mm512_cvtepu8_epi32(simde_mm_loadu_si128((const simde__m128i *) (const void *) p));
}
RSIMD_INLINE void rsimd_vi32_storeu_u8(uint8_t *p, rsimd_vi32 v) {
  const int32_t z = (int32_t) 0x80808080u, b = 0x0C080400;
  simde__m512i g = simde_mm512_shuffle_epi8(
    v, simde_mm512_set_epi32(z, z, z, b, z, z, z, b, z, z, z, b, z, z, z, b));
  g = simde_mm512_permutexvar_epi32(
    simde_mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 8, 4, 0), g);
  simde_mm_storeu_si128((simde__m128i *) (void *) p, simde_mm512_castsi512_si128(g));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_u8_p(rsimd_p32 pg, const uint8_t *p, uint8_t fill) {
  uint8_t buf[16];
  int j;
  for (j = 0; j < 16; j++) buf[j] = (pg >> j) & 1 ? p[j] : fill;
  return rsimd_vi32_loadu_u8(buf);
}
RSIMD_INLINE void rsimd_vi32_storeu_u8_p(rsimd_p32 pg, uint8_t *p, rsimd_vi32 v) {
  uint8_t buf[16];
  int j;
  rsimd_vi32_storeu_u8(buf, v);
  for (j = 0; j < 16; j++) {
    if ((pg >> j) & 1) p[j] = buf[j];
  }
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32(const int32_t *p) {
  return simde_mm512_cvtepi32_epi64(simde_mm256_loadu_si256((const simde__m256i *) (const void *) p));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return simde_mm512_cvtepi32_epi64(simde_mm256_mask_loadu_epi32(simde_mm256_set1_epi32(fill), pg, p));
}
RSIMD_INLINE void rsimd_vi64_storeu_i32(int32_t *p, rsimd_vi64 v) {
  simde_mm256_storeu_si256((simde__m256i *) (void *) p, simde_mm512_cvtepi64_epi32(v));
}
RSIMD_INLINE void rsimd_vi64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vi64 v) {
  simde_mm256_mask_storeu_epi32(p, pg, simde_mm512_cvtepi64_epi32(v));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_floor(rsimd_vf64 a) {
  return simde_mm512_roundscale_pd(a, SIMDE_MM_FROUND_TO_NEG_INF | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_ceil(rsimd_vf64 a) {
  return simde_mm512_roundscale_pd(a, SIMDE_MM_FROUND_TO_POS_INF | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_trunc(rsimd_vf64 a) {
  return simde_mm512_roundscale_pd(a, SIMDE_MM_FROUND_TO_ZERO | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_rint(rsimd_vf64 a) {
  return simde_mm512_roundscale_pd(a, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32(const int32_t *p) {
  return simde_mm512_cvtepi64_pd(rsimd_vi64_loadu_i32(p));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return simde_mm512_cvtepi64_pd(rsimd_vi64_loadu_i32_p(pg, p, fill));
}
/* The vendored SIMDe subset has no 512-bit double -> int32 conversion, so
   convert the two 256-bit halves. */
RSIMD_INLINE simde__m256i rsimd_s512_cvttpd_epi32(rsimd_vf64 v) {
  return simde_mm256_set_m128i(simde_mm256_cvttpd_epi32(simde_mm512_extractf64x4_pd(v, 1)),
                               simde_mm256_cvttpd_epi32(simde_mm512_castpd512_pd256(v)));
}
RSIMD_INLINE void rsimd_vf64_storeu_i32(int32_t *p, rsimd_vf64 v) {
  simde_mm256_storeu_si256((simde__m256i *) (void *) p, rsimd_s512_cvttpd_epi32(v));
}
RSIMD_INLINE void rsimd_vf64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vf64 v) {
  simde_mm256_mask_storeu_epi32(p, pg, rsimd_s512_cvttpd_epi32(v));
}
