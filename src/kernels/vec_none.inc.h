/* Vector layer of the none tier: one lane, plain C, no SIMDe. This makes
   the none tier compile the same kernel sources as every other tier while
   staying obviously correct. Included by kernels/common.inc.h only. */

#define RSIMD_WIDTH_F64 1
#define RSIMD_WIDTH_I64 1
#define RSIMD_WIDTH_I32 1
#define RSIMD_LANES_64 ((ptrdiff_t) 1)
#define RSIMD_LANES_32 ((ptrdiff_t) 1)

typedef double rsimd_vf64;
typedef int32_t rsimd_vi32;
typedef int64_t rsimd_vi64;
typedef int rsimd_mf64;
typedef int rsimd_mi32;
typedef int rsimd_mi64;
typedef int rsimd_p64; /* 1 if the lane is active */
typedef int rsimd_p32;

/* Predicates */
RSIMD_INLINE rsimd_p64 rsimd_p64_while(ptrdiff_t i, ptrdiff_t n) { return i < n; }
RSIMD_INLINE rsimd_p64 rsimd_p64_true(void) { return 1; }
RSIMD_INLINE int rsimd_p64_count(rsimd_p64 pg) { return pg; }
RSIMD_INLINE rsimd_p32 rsimd_p32_while(ptrdiff_t i, ptrdiff_t n) { return i < n; }
RSIMD_INLINE rsimd_p32 rsimd_p32_true(void) { return 1; }
RSIMD_INLINE int rsimd_p32_count(rsimd_p32 pg) { return pg; }

/* Masks */
#define RSIMD_NONE_MASK_OPS(m)                                                   \
  RSIMD_INLINE rsimd_##m rsimd_##m##_and(rsimd_##m a, rsimd_##m b) { return a && b; } \
  RSIMD_INLINE rsimd_##m rsimd_##m##_or(rsimd_##m a, rsimd_##m b) { return a || b; }  \
  RSIMD_INLINE rsimd_##m rsimd_##m##_not(rsimd_##m a) { return !a; }                  \
  RSIMD_INLINE rsimd_##m rsimd_##m##_andnot(rsimd_##m a, rsimd_##m b) { return !a && b; } \
  RSIMD_INLINE int rsimd_##m##_any(rsimd_##m a) { return a != 0; }                     \
  RSIMD_INLINE int rsimd_##m##_all(rsimd_##m a) { return a != 0; }                     \
  RSIMD_INLINE int rsimd_##m##_count(rsimd_##m a) { return a != 0; }
RSIMD_NONE_MASK_OPS(mf64)
RSIMD_NONE_MASK_OPS(mi32)
RSIMD_NONE_MASK_OPS(mi64)
#undef RSIMD_NONE_MASK_OPS

RSIMD_INLINE rsimd_mf64 rsimd_mi64_to_mf64(rsimd_mi64 m) { return m; }
RSIMD_INLINE rsimd_mi64 rsimd_mf64_to_mi64(rsimd_mf64 m) { return m; }

/* Doubles */
RSIMD_INLINE uint64_t rsimd_none_f64_bits(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return u;
}
RSIMD_INLINE double rsimd_none_bits_f64(uint64_t u) {
  double x;
  memcpy(&x, &u, sizeof x);
  return x;
}

RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu(const double *p) { return *p; }
RSIMD_INLINE void rsimd_vf64_storeu(double *p, rsimd_vf64 v) { *p = v; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_p(rsimd_p64 pg, const double *p, double fill) {
  return pg ? *p : fill;
}
RSIMD_INLINE void rsimd_vf64_storeu_p(rsimd_p64 pg, double *p, rsimd_vf64 v) {
  if (pg) *p = v;
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_set1(double x) { return x; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_zero(void) { return 0.0; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_add(rsimd_vf64 a, rsimd_vf64 b) { return a + b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sub(rsimd_vf64 a, rsimd_vf64 b) { return a - b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_mul(rsimd_vf64 a, rsimd_vf64 b) { return a * b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_div(rsimd_vf64 a, rsimd_vf64 b) { return a / b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_fma(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c) {
  return fma(a, b, c);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_min(rsimd_vf64 a, rsimd_vf64 b) { return a < b ? a : b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_max(rsimd_vf64 a, rsimd_vf64 b) { return a > b ? a : b; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_abs(rsimd_vf64 a) { return fabs(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_neg(rsimd_vf64 a) { return -a; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sqrt(rsimd_vf64 a) { return sqrt(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_and(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_none_bits_f64(rsimd_none_f64_bits(a) & rsimd_none_f64_bits(b));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_or(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_none_bits_f64(rsimd_none_f64_bits(a) | rsimd_none_f64_bits(b));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_xor(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_none_bits_f64(rsimd_none_f64_bits(a) ^ rsimd_none_f64_bits(b));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_andnot(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_none_bits_f64(~rsimd_none_f64_bits(a) & rsimd_none_f64_bits(b));
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_eq(rsimd_vf64 a, rsimd_vf64 b) { return a == b; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ne(rsimd_vf64 a, rsimd_vf64 b) { return !(a == b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_lt(rsimd_vf64 a, rsimd_vf64 b) { return a < b; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_le(rsimd_vf64 a, rsimd_vf64 b) { return a <= b; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_gt(rsimd_vf64 a, rsimd_vf64 b) { return a > b; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ge(rsimd_vf64 a, rsimd_vf64 b) { return a >= b; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_nan(rsimd_vf64 a) { return isnan(a) != 0; }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_na(rsimd_vf64 a) {
  return isnan(a) && (rsimd_none_f64_bits(a) & 0xFFFFFFFFu) == RSIMD_NA_LOW_WORD;
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_blend(rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m) {
  return m ? b : a;
}
RSIMD_INLINE double rsimd_vf64_reduce_add(rsimd_vf64 a) { return a; }
RSIMD_INLINE double rsimd_vf64_reduce_min(rsimd_vf64 a) { return a; }
RSIMD_INLINE double rsimd_vf64_reduce_max(rsimd_vf64 a) { return a; }
RSIMD_INLINE rsimd_vi64 rsimd_vf64_as_vi64(rsimd_vf64 a) {
  return (int64_t) rsimd_none_f64_bits(a);
}
RSIMD_INLINE rsimd_vf64 rsimd_vi64_as_vf64(rsimd_vi64 a) {
  return rsimd_none_bits_f64((uint64_t) a);
}

/* 32-bit integers; arithmetic goes through unsigned to wrap without
   undefined behaviour. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu(const int32_t *p) { return *p; }
RSIMD_INLINE void rsimd_vi32_storeu(int32_t *p, rsimd_vi32 v) { *p = v; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_p(rsimd_p32 pg, const int32_t *p, int32_t fill) {
  return pg ? *p : fill;
}
RSIMD_INLINE void rsimd_vi32_storeu_p(rsimd_p32 pg, int32_t *p, rsimd_vi32 v) {
  if (pg) *p = v;
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_set1(int32_t x) { return x; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_zero(void) { return 0; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_add(rsimd_vi32 a, rsimd_vi32 b) {
  return (int32_t) ((uint32_t) a + (uint32_t) b);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sub(rsimd_vi32 a, rsimd_vi32 b) {
  return (int32_t) ((uint32_t) a - (uint32_t) b);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mul(rsimd_vi32 a, rsimd_vi32 b) {
  return (int32_t) ((uint32_t) a * (uint32_t) b);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_min(rsimd_vi32 a, rsimd_vi32 b) { return a < b ? a : b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_max(rsimd_vi32 a, rsimd_vi32 b) { return a > b ? a : b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_and(rsimd_vi32 a, rsimd_vi32 b) { return a & b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_or(rsimd_vi32 a, rsimd_vi32 b) { return a | b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_xor(rsimd_vi32 a, rsimd_vi32 b) { return a ^ b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_andnot(rsimd_vi32 a, rsimd_vi32 b) { return ~a & b; }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_eq(rsimd_vi32 a, rsimd_vi32 b) { return a == b; }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_gt(rsimd_vi32 a, rsimd_vi32 b) { return a > b; }
RSIMD_INLINE rsimd_mi32 rsimd_vi32_cmp_lt(rsimd_vi32 a, rsimd_vi32 b) { return a < b; }
RSIMD_INLINE rsimd_vi32 rsimd_vi32_blend(rsimd_vi32 a, rsimd_vi32 b, rsimd_mi32 m) {
  return m ? b : a;
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_is_na(rsimd_vi32 a) { return a == INT32_MIN; }
RSIMD_INLINE int64_t rsimd_vi32_reduce_add(rsimd_vi32 a) { return a; }
RSIMD_INLINE int32_t rsimd_vi32_reduce_min(rsimd_vi32 a) { return a; }
RSIMD_INLINE int32_t rsimd_vi32_reduce_max(rsimd_vi32 a) { return a; }

/* 64-bit integers */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu(const int64_t *p) { return *p; }
RSIMD_INLINE void rsimd_vi64_storeu(int64_t *p, rsimd_vi64 v) { *p = v; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_p(rsimd_p64 pg, const int64_t *p, int64_t fill) {
  return pg ? *p : fill;
}
RSIMD_INLINE void rsimd_vi64_storeu_p(rsimd_p64 pg, int64_t *p, rsimd_vi64 v) {
  if (pg) *p = v;
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_set1(int64_t x) { return x; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_zero(void) { return 0; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_add(rsimd_vi64 a, rsimd_vi64 b) {
  return (int64_t) ((uint64_t) a + (uint64_t) b);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sub(rsimd_vi64 a, rsimd_vi64 b) {
  return (int64_t) ((uint64_t) a - (uint64_t) b);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_mul(rsimd_vi64 a, rsimd_vi64 b) {
  return (int64_t) ((uint64_t) a * (uint64_t) b);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_min(rsimd_vi64 a, rsimd_vi64 b) { return a < b ? a : b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_max(rsimd_vi64 a, rsimd_vi64 b) { return a > b ? a : b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_and(rsimd_vi64 a, rsimd_vi64 b) { return a & b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_or(rsimd_vi64 a, rsimd_vi64 b) { return a | b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_xor(rsimd_vi64 a, rsimd_vi64 b) { return a ^ b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_andnot(rsimd_vi64 a, rsimd_vi64 b) { return ~a & b; }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_eq(rsimd_vi64 a, rsimd_vi64 b) { return a == b; }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_gt(rsimd_vi64 a, rsimd_vi64 b) { return a > b; }
RSIMD_INLINE rsimd_mi64 rsimd_vi64_cmp_lt(rsimd_vi64 a, rsimd_vi64 b) { return a < b; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_blend(rsimd_vi64 a, rsimd_vi64 b, rsimd_mi64 m) {
  return m ? b : a;
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_is_na(rsimd_vi64 a) { return a == INT64_MIN; }
RSIMD_INLINE int64_t rsimd_vi64_reduce_add(rsimd_vi64 a) { return a; }
RSIMD_INLINE int64_t rsimd_vi64_reduce_min(rsimd_vi64 a) { return a; }
RSIMD_INLINE int64_t rsimd_vi64_reduce_max(rsimd_vi64 a) { return a; }

/* Width conversions and extras */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mulhi(rsimd_vi32 a, rsimd_vi32 b) {
  return (int32_t) (((int64_t) a * b) >> 32);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32(const int32_t *p) { return *p; }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return pg ? *p : fill;
}
RSIMD_INLINE void rsimd_vi64_storeu_i32(int32_t *p, rsimd_vi64 v) {
  *p = (int32_t) (uint32_t) (uint64_t) v;
}
RSIMD_INLINE void rsimd_vi64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vi64 v) {
  if (pg) rsimd_vi64_storeu_i32(p, v);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_floor(rsimd_vf64 a) { return floor(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32(const int32_t *p) { return (double) *p; }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return (double) (pg ? *p : fill);
}
RSIMD_INLINE void rsimd_vf64_storeu_i32(int32_t *p, rsimd_vf64 v) { *p = (int32_t) v; }
RSIMD_INLINE void rsimd_vf64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vf64 v) {
  if (pg) *p = (int32_t) v;
}
