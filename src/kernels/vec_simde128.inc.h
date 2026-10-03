/* 128-bit vector layer over SIMDe's SSE intrinsics, for the sse2 tier
   (native SSE2 on x86) and the neon tier (SIMDe compiles the same calls to
   NEON). Calls above SSE2 (SSE4.1 min/max, 64-bit compares ...) are native
   on AArch64 and use SIMDe's SSE2 implementations on x86. Included by
   kernels/common.inc.h only. */

#include "x86/sse2.h"
#include "x86/sse4.2.h"
#if defined(__aarch64__) || defined(_M_ARM64)
#include "x86/fma.h"
#endif

#define RSIMD_WIDTH_F64 2
#define RSIMD_WIDTH_I64 2
#define RSIMD_WIDTH_I32 4
#define RSIMD_LANES_64 ((ptrdiff_t) 2)
#define RSIMD_LANES_32 ((ptrdiff_t) 4)

typedef simde__m128i rsimd_vi32;
typedef simde__m128i rsimd_vi64;
typedef simde__m128i rsimd_mi32; /* all-ones lanes where true */
typedef simde__m128i rsimd_mi64;

RSIMD_INLINE simde__m128i rsimd_s128_ones(void) { return simde_mm_set1_epi32(-1); }

/* Doubles */
#ifndef RSIMD_NO_F64_SIMD
typedef simde__m128d rsimd_vf64;
typedef simde__m128d rsimd_mf64;

RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu(const double *p) { return simde_mm_loadu_pd(p); }
RSIMD_INLINE void rsimd_vf64_storeu(double *p, rsimd_vf64 v) { simde_mm_storeu_pd(p, v); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_set1(double x) { return simde_mm_set1_pd(x); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_zero(void) { return simde_mm_setzero_pd(); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_add(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_add_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sub(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_sub_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_mul(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_mul_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_div(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_div_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_fma(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c) {
#if defined(__aarch64__) || defined(_M_ARM64)
  return simde_mm_fmadd_pd(a, b, c);
#else
  /* SSE2 has no fused multiply-add and SIMDe's fallback rounds twice. */
  double x[2], y[2], z[2];
  simde_mm_storeu_pd(x, a);
  simde_mm_storeu_pd(y, b);
  simde_mm_storeu_pd(z, c);
  x[0] = fma(x[0], y[0], z[0]);
  x[1] = fma(x[1], y[1], z[1]);
  return simde_mm_loadu_pd(x);
#endif
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_min(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_min_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_max(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_max_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_abs(rsimd_vf64 a) {
  return simde_mm_andnot_pd(simde_mm_set1_pd(-0.0), a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_neg(rsimd_vf64 a) {
  return simde_mm_xor_pd(simde_mm_set1_pd(-0.0), a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sqrt(rsimd_vf64 a) { return simde_mm_sqrt_pd(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_and(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_and_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_or(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_or_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_xor(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_xor_pd(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_andnot(rsimd_vf64 a, rsimd_vf64 b) {
  return simde_mm_andnot_pd(a, b);
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_eq(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmpeq_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ne(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmpneq_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_lt(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmplt_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_le(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmple_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_gt(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmpgt_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ge(rsimd_vf64 a, rsimd_vf64 b) { return simde_mm_cmpge_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_nan(rsimd_vf64 a) { return simde_mm_cmpunord_pd(a, a); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_na(rsimd_vf64 a) {
  /* Compare the 32-bit words with 1954, then copy each low-word result over
     its 64-bit lane: SSE2 has no 64-bit compare. */
  simde__m128i lo = simde_mm_cmpeq_epi32(simde_mm_castpd_si128(a),
                                         simde_mm_set1_epi64x(RSIMD_NA_LOW_WORD));
  lo = simde_mm_shuffle_epi32(lo, SIMDE_MM_SHUFFLE(2, 2, 0, 0));
  return simde_mm_and_pd(simde_mm_castsi128_pd(lo), simde_mm_cmpunord_pd(a, a));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_blend(rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m) {
  return simde_mm_or_pd(simde_mm_and_pd(m, b), simde_mm_andnot_pd(m, a));
}
RSIMD_INLINE rsimd_vi64 rsimd_vf64_as_vi64(rsimd_vf64 a) { return simde_mm_castpd_si128(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vi64_as_vf64(rsimd_vi64 a) { return simde_mm_castsi128_pd(a); }

RSIMD_INLINE rsimd_mf64 rsimd_mf64_and(rsimd_mf64 a, rsimd_mf64 b) { return simde_mm_and_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_mf64_or(rsimd_mf64 a, rsimd_mf64 b) { return simde_mm_or_pd(a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_mf64_not(rsimd_mf64 a) {
  return simde_mm_xor_pd(a, simde_mm_castsi128_pd(rsimd_s128_ones()));
}
RSIMD_INLINE rsimd_mf64 rsimd_mf64_andnot(rsimd_mf64 a, rsimd_mf64 b) {
  return simde_mm_andnot_pd(a, b);
}
RSIMD_INLINE int rsimd_mf64_any(rsimd_mf64 a) { return simde_mm_movemask_pd(a) != 0; }
RSIMD_INLINE int rsimd_mf64_all(rsimd_mf64 a) { return simde_mm_movemask_pd(a) == 3; }
RSIMD_INLINE int rsimd_mf64_count(rsimd_mf64 a) {
  return rsimd_popcount32((uint32_t) simde_mm_movemask_pd(a));
}
RSIMD_INLINE rsimd_mf64 rsimd_mi64_to_mf64(rsimd_mi64 m) { return simde_mm_castsi128_pd(m); }
RSIMD_INLINE rsimd_mi64 rsimd_mf64_to_mi64(rsimd_mf64 m) { return simde_mm_castpd_si128(m); }
#endif /* RSIMD_NO_F64_SIMD */

/* Bitwise operations on the integer register type. */
RSIMD_INLINE simde__m128i rsimd_s128_blend(simde__m128i a, simde__m128i b, simde__m128i m) {
  return simde_mm_or_si128(simde_mm_and_si128(m, b), simde_mm_andnot_si128(m, a));
}

#define RSIMD_S128_BITWISE(v)                                                     \
  RSIMD_INLINE rsimd_##v rsimd_##v##_and(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm_and_si128(a, b);                                             \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_or(rsimd_##v a, rsimd_##v b) {               \
    return simde_mm_or_si128(a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_xor(rsimd_##v a, rsimd_##v b) {              \
    return simde_mm_xor_si128(a, b);                                             \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_andnot(rsimd_##v a, rsimd_##v b) {           \
    return simde_mm_andnot_si128(a, b);                                          \
  }
#define RSIMD_S128_MASK(m)                                                        \
  RSIMD_INLINE rsimd_##m rsimd_##m##_and(rsimd_##m a, rsimd_##m b) {              \
    return simde_mm_and_si128(a, b);                                             \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_or(rsimd_##m a, rsimd_##m b) {               \
    return simde_mm_or_si128(a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_not(rsimd_##m a) {                           \
    return simde_mm_xor_si128(a, rsimd_s128_ones());                             \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_andnot(rsimd_##m a, rsimd_##m b) {           \
    return simde_mm_andnot_si128(a, b);                                          \
  }

/* 32-bit integers */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu(const int32_t *p) {
  return simde_mm_loadu_si128((const simde__m128i *) (const void *) p);
}
RSIMD_INLINE void rsimd_vi32_storeu(int32_t *p, rsimd_vi32 v) {
  simde_mm_storeu_si128((simde__m128i *) (void *) p, v);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_set1(int32_t x) { return simde_mm_set1_epi32(x); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_zero(void) { return simde_mm_setzero_si128(); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_add(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_add_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sub(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_sub_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mul(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_mullo_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_min(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_min_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_max(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_max_epi32(a, b); }
RSIMD_S128_BITWISE(vi32)
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_eq(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_cmpeq_epi32(a, b); }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_gt(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_cmpgt_epi32(a, b); }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_lt(rsimd_vi32 a, rsimd_vi32 b) { return simde_mm_cmplt_epi32(a, b); }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_blend(rsimd_vi32 a, rsimd_vi32 b, rsimd_mi32 m) {
  return rsimd_s128_blend(a, b, m);
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_is_na(rsimd_vi32 a) {
  return simde_mm_cmpeq_epi32(a, simde_mm_set1_epi32(INT32_MIN));
}
RSIMD_S128_MASK(mi32)
RSIMD_INLINE int rsimd_mi32_any(rsimd_mi32 a) { return simde_mm_movemask_epi8(a) != 0; }
RSIMD_INLINE int rsimd_mi32_all(rsimd_mi32 a) { return simde_mm_movemask_epi8(a) == 0xFFFF; }
RSIMD_INLINE int rsimd_mi32_count(rsimd_mi32 a) {
  return rsimd_popcount32((uint32_t) simde_mm_movemask_ps(simde_mm_castsi128_ps(a)));
}

/* 64-bit integers */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu(const int64_t *p) {
  return simde_mm_loadu_si128((const simde__m128i *) (const void *) p);
}
RSIMD_INLINE void rsimd_vi64_storeu(int64_t *p, rsimd_vi64 v) {
  simde_mm_storeu_si128((simde__m128i *) (void *) p, v);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_set1(int64_t x) { return simde_mm_set1_epi64x(x); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_zero(void) { return simde_mm_setzero_si128(); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_add(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm_add_epi64(a, b); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sub(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm_sub_epi64(a, b); }
RSIMD_S128_BITWISE(vi64)
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_eq(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm_cmpeq_epi64(a, b); }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_gt(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm_cmpgt_epi64(a, b); }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_lt(rsimd_vi64 a, rsimd_vi64 b) { return simde_mm_cmpgt_epi64(b, a); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_blend(rsimd_vi64 a, rsimd_vi64 b, rsimd_mi64 m) {
  return rsimd_s128_blend(a, b, m);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_min(rsimd_vi64 a, rsimd_vi64 b) {
  return rsimd_s128_blend(b, a, simde_mm_cmpgt_epi64(b, a));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_max(rsimd_vi64 a, rsimd_vi64 b) {
  return rsimd_s128_blend(b, a, simde_mm_cmpgt_epi64(a, b));
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_is_na(rsimd_vi64 a) {
  return simde_mm_cmpeq_epi64(a, simde_mm_set1_epi64x(INT64_MIN));
}
RSIMD_S128_MASK(mi64)
RSIMD_INLINE int rsimd_mi64_any(rsimd_mi64 a) { return simde_mm_movemask_epi8(a) != 0; }
RSIMD_INLINE int rsimd_mi64_all(rsimd_mi64 a) { return simde_mm_movemask_epi8(a) == 0xFFFF; }
RSIMD_INLINE int rsimd_mi64_count(rsimd_mi64 a) {
  return rsimd_popcount32((uint32_t) simde_mm_movemask_epi8(a)) / 8;
}

#undef RSIMD_S128_BITWISE
#undef RSIMD_S128_MASK

#include "vec_fixed.inc.h"
