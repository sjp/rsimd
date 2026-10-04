/* 256-bit vector layer over SIMDe's AVX2 and FMA intrinsics, for the avx2
   tier. Included by kernels/common.inc.h only. */

#include "x86/avx2.h"
#include "x86/fma.h"

#define RSIMD_WIDTH_F64 4
#define RSIMD_WIDTH_I64 4
#define RSIMD_WIDTH_I32 8
#define RSIMD_LANES_64 ((ptrdiff_t) 4)
#define RSIMD_LANES_32 ((ptrdiff_t) 8)

typedef simde__m256d rsimd_vf64;
typedef simde__m256d rsimd_mf64; /* all-ones lanes where true */
typedef simde__m256i rsimd_vi32;
typedef simde__m256i rsimd_vi64;
typedef simde__m256i rsimd_mi32;
typedef simde__m256i rsimd_mi64;

RSIMD_INLINE simde__m256i rsimd_s256_ones(void) { return simde_mm256_set1_epi32(-1); }

/* Doubles */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu(const double *p) { return simde_mm256_loadu_pd(p); }
RSIMD_INLINE void rsimd_vf64_storeu(double *p, rsimd_vf64 v) { simde_mm256_storeu_pd(p, v); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_set1(double x) { return simde_mm256_set1_pd(x); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_zero(void) { return simde_mm256_setzero_pd(); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_add(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_add_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sub(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_sub_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_mul(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_mul_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_div(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_div_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_fma(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c) {
  return simde_mm256_fmadd_pd(a, b, c);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_min(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_min_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_max(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_max_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_abs(rsimd_vf64 a) {
  return simde_mm256_andnot_pd(simde_mm256_set1_pd(-0.0), a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_neg(rsimd_vf64 a) {
  return simde_mm256_xor_pd(simde_mm256_set1_pd(-0.0), a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sqrt(rsimd_vf64 a) { return simde_mm256_sqrt_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_and(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_and_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_or(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_or_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_xor(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm256_xor_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_andnot(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_andnot_pd(a, b);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_eq(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_EQ_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ne(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_NEQ_UQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_lt(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_LT_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_le(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_LE_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_gt(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_GT_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ge(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm256_cmp_pd(a, b, SIMDE_CMP_GE_OQ);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_nan(rsimd_vf64 a) {
  return simde_mm256_cmp_pd(a, a, SIMDE_CMP_UNORD_Q);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_na(rsimd_vf64 a) {
  simde__m256i lo = simde_mm256_cmpeq_epi64(
    simde_mm256_and_si256(simde_mm256_castpd_si256(a), simde_mm256_set1_epi64x(0xFFFFFFFF)),
    simde_mm256_set1_epi64x(RSIMD_NA_LOW_WORD));
  return simde_mm256_and_pd(simde_mm256_castsi256_pd(lo), rsimd_vf64_is_nan(a));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_blend(rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m) {
  return simde_mm256_blendv_pd(a, b, m);
}
/* Lane moves for prefix scans (see the 128-bit layer); k is 1 or 2. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_shift_up(rsimd_vf64 v, int k, rsimd_vf64 fill) {
  if (k == 1) {
    return simde_mm256_blend_pd(simde_mm256_permute4x64_pd(v, SIMDE_MM_SHUFFLE(2, 1, 0, 3)), fill,
                                0x1);
  }
  return simde_mm256_blend_pd(simde_mm256_permute4x64_pd(v, SIMDE_MM_SHUFFLE(1, 0, 3, 2)), fill,
                              0x3);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_bcast_last(rsimd_vf64 v) {
  return simde_mm256_permute4x64_pd(v, 0xFF);
}
RSIMD_INLINE double rsimd_vf64_first(rsimd_vf64 v) { return simde_mm256_cvtsd_f64(v); }
RSIMD_INLINE rsimd_vi64 rsimd_vf64_as_vi64(rsimd_vf64 a) { return simde_mm256_castpd_si256(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vi64_as_vf64(rsimd_vi64 a) { return simde_mm256_castsi256_pd(a); }

RSIMD_INLINE rsimd_mf64 rsimd_mf64_and(rsimd_mf64 a, rsimd_mf64 b) { return simde_mm256_and_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_mf64_or(rsimd_mf64 a, rsimd_mf64 b) { return simde_mm256_or_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_mf64_not(rsimd_mf64 a) {
  return simde_mm256_xor_pd(a, simde_mm256_castsi256_pd(rsimd_s256_ones()));
}
RSIMD_INLINE rsimd_mf64 rsimd_mf64_andnot(rsimd_mf64 a, rsimd_mf64 b) {
  return simde_mm256_andnot_pd(a, b);
}
RSIMD_INLINE int rsimd_mf64_any(rsimd_mf64 a) { return simde_mm256_movemask_pd(a) != 0; }
RSIMD_INLINE int rsimd_mf64_all(rsimd_mf64 a) { return simde_mm256_movemask_pd(a) == 0xF; }
RSIMD_INLINE int rsimd_mf64_count(rsimd_mf64 a) {
  return rsimd_popcount32((uint32_t) simde_mm256_movemask_pd(a));
}
RSIMD_INLINE rsimd_mf64 rsimd_mi64_to_mf64(rsimd_mi64 m) { return simde_mm256_castsi256_pd(m); }
RSIMD_INLINE rsimd_mi64 rsimd_mf64_to_mi64(rsimd_mf64 m) { return simde_mm256_castpd_si256(m); }

#define RSIMD_S256_BITWISE(v)                                                     \
  RSIMD_INLINE rsimd_##v rsimd_##v##_and(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm256_and_si256(a, b);                                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_or(rsimd_##v a, rsimd_##v b) {               \
    return simde_mm256_or_si256(a, b);                                           \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_xor(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm256_xor_si256(a, b);                                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_andnot(rsimd_##v a, rsimd_##v b) {           \
    return simde_mm256_andnot_si256(a, b);                                       \
  }
#define RSIMD_S256_MASK(m)                                                        \
  RSIMD_INLINE rsimd_##m rsimd_##m##_and(rsimd_##m a, rsimd_##m b) {              \
    return simde_mm256_and_si256(a, b);                                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_or(rsimd_##m a, rsimd_##m b) {               \
    return simde_mm256_or_si256(a, b);                                           \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_not(rsimd_##m a) {                           \
    return simde_mm256_xor_si256(a, rsimd_s256_ones());                          \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_andnot(rsimd_##m a, rsimd_##m b) {           \
    return simde_mm256_andnot_si256(a, b);                                       \
  }                                                                              \
  RSIMD_INLINE int rsimd_##m##_any(rsimd_##m a) {                                 \
    return simde_mm256_movemask_epi8(a) != 0;                                    \
  }                                                                              \
  RSIMD_INLINE int rsimd_##m##_all(rsimd_##m a) {                                 \
    return simde_mm256_movemask_epi8(a) == -1;                                   \
  }

/* 32-bit integers */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu(const int32_t *p) {
  return simde_mm256_loadu_si256((const simde__m256i *) (const void *) p);
}
RSIMD_INLINE void rsimd_vi32_storeu(int32_t *p, rsimd_vi32 v) {
  simde_mm256_storeu_si256((simde__m256i *) (void *) p, v);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_set1(int32_t x) { return simde_mm256_set1_epi32(x); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_zero(void) { return simde_mm256_setzero_si256(); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_add(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_add_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sub(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_sub_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mul(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_mullo_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_min(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_min_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_max(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_max_epi32(a, b); }
RSIMD_S256_BITWISE(vi32)
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_eq(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_cmpeq_epi32(a, b); }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_gt(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_cmpgt_epi32(a, b); }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_lt(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm256_cmpgt_epi32(b, a); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_blend(rsimd_vi32 a, rsimd_vi32 b, rsimd_mi32 m) {
  return simde_mm256_blendv_epi8(a, b, m);
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_is_na(rsimd_vi32 a) {
  return simde_mm256_cmpeq_epi32(a, simde_mm256_set1_epi32(INT32_MIN));
}
RSIMD_S256_MASK(mi32)
RSIMD_INLINE int rsimd_mi32_count(rsimd_mi32 a) {
  return rsimd_popcount32((uint32_t) simde_mm256_movemask_ps(simde_mm256_castsi256_ps(a)));
}

/* 64-bit integers */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu(const int64_t *p) {
  return simde_mm256_loadu_si256((const simde__m256i *) (const void *) p);
}
RSIMD_INLINE void rsimd_vi64_storeu(int64_t *p, rsimd_vi64 v) {
  simde_mm256_storeu_si256((simde__m256i *) (void *) p, v);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_set1(int64_t x) { return simde_mm256_set1_epi64x(x); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_zero(void) { return simde_mm256_setzero_si256(); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_add(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm256_add_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sub(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm256_sub_epi64(a, b); }
RSIMD_S256_BITWISE(vi64)
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_eq(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm256_cmpeq_epi64(a, b); }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_gt(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm256_cmpgt_epi64(a, b); }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_lt(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm256_cmpgt_epi64(b, a); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_blend(rsimd_vi64 a, rsimd_vi64 b, rsimd_mi64 m) {
  return simde_mm256_blendv_epi8(a, b, m);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_min(rsimd_vi64 a, rsimd_vi64 b) {
  return simde_mm256_blendv_epi8(b, a, simde_mm256_cmpgt_epi64(b, a));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_max(rsimd_vi64 a, rsimd_vi64 b) {
  return simde_mm256_blendv_epi8(b, a, simde_mm256_cmpgt_epi64(a, b));
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_is_na(rsimd_vi64 a) {
  return simde_mm256_cmpeq_epi64(a, simde_mm256_set1_epi64x(INT64_MIN));
}
RSIMD_S256_MASK(mi64)
RSIMD_INLINE int rsimd_mi64_count(rsimd_mi64 a) {
  return rsimd_popcount32((uint32_t) simde_mm256_movemask_pd(simde_mm256_castsi256_pd(a)));
}

#undef RSIMD_S256_BITWISE
#undef RSIMD_S256_MASK

/* Width conversions and extras. The 64-bit-lane conversions move
   RSIMD_WIDTH_I64 (4) int32 elements. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mulhi(rsimd_vi32 a, rsimd_vi32 b) {
  /* See the 128-bit layer. */
  simde__m256i even = simde_mm256_mul_epi32(a, b);
  simde__m256i odd =
    simde_mm256_mul_epi32(simde_mm256_srli_epi64(a, 32), simde_mm256_srli_epi64(b, 32));
  return simde_mm256_blend_epi32(simde_mm256_srli_epi64(even, 32), odd, 0xAA);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32(const int32_t *p) {
  return simde_mm256_cvtepi32_epi64(simde_mm_loadu_si128((const simde__m128i *) (const void *) p));
}
RSIMD_INLINE void rsimd_vi64_storeu_i32(int32_t *p, rsimd_vi64 v) {
  simde__m256i low = simde_mm256_permutevar8x32_epi32(v, simde_mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
  simde_mm_storeu_si128((simde__m128i *) (void *) p, simde_mm256_castsi256_si128(low));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_floor(rsimd_vf64 a) { return simde_mm256_floor_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_ceil(rsimd_vf64 a) { return simde_mm256_ceil_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_trunc(rsimd_vf64 a) {
  return simde_mm256_round_pd(a, SIMDE_MM_FROUND_TO_ZERO | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_rint(rsimd_vf64 a) {
  return simde_mm256_round_pd(a, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32(const int32_t *p) {
  return simde_mm256_cvtepi32_pd(simde_mm_loadu_si128((const simde__m128i *) (const void *) p));
}
RSIMD_INLINE void rsimd_vf64_storeu_i32(int32_t *p, rsimd_vf64 v) {
  simde_mm_storeu_si128((simde__m128i *) (void *) p, simde_mm256_cvttpd_epi32(v));
}

#include "vec_fixed.inc.h"
