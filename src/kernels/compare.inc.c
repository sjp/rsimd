/* Comparison kernels: x op y (op codes RSIMD_CMP_* in kernel_types.h) as
 * logical vectors with base R's missing values: NA where either operand is
 * missing (any NaN for doubles, INT32_MIN for int32 elements). Doubles
 * compare as IEEE numbers, so -0 == 0. Bytes (raw) compare as unsigned
 * values and are never missing.
 *
 * cmp_f64 reads each operand as doubles or int32 elements (flag
 * RSIMD_EW_I32(k)), converting int32 exactly with NA becoming NA_real_.
 * RSIMD_EW_SCALAR(k) broadcasts element 0 of operand k.
 */

#include "logical.inc.h"

static inline int rsimd_cmp_1(int op, double a, double b) {
  switch (op) {
  case RSIMD_CMP_EQ: return a == b;
  case RSIMD_CMP_NE: return a != b;
  case RSIMD_CMP_LT: return a < b;
  case RSIMD_CMP_LE: return a <= b;
  case RSIMD_CMP_GT: return a > b;
  default: return a >= b;
  }
}

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(cmp_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    double a = rsimd_ew_get(x, flags, 0, i, 1), b = rsimd_ew_get(y, flags, 1, i, 1);
    out[i] = isnan(a) || isnan(b) ? RSIMD_NA_I32 : rsimd_cmp_1(op, a, b);
  }
}

void RSIMD_KERNEL(cmp_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                           int *out) {
  const R_xlen_t sx = RSIMD_LGL_SCALAR(0) ? 0 : 1, sy = RSIMD_LGL_SCALAR(1) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int a = x[i * sx], b = y[i * sy];
    out[i] = a == RSIMD_NA_I32 || b == RSIMD_NA_I32 ? RSIMD_NA_I32 : rsimd_cmp_1(op, a, b);
  }
}

void RSIMD_KERNEL(cmp_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                          int *out);
void RSIMD_KERNEL(cmp_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                          int *out) {
  const R_xlen_t sx = RSIMD_LGL_SCALAR(0) ? 0 : 1, sy = RSIMD_LGL_SCALAR(1) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) out[i] = rsimd_cmp_1(op, x[i * sx], y[i * sy]);
}

#else /* vector tiers */

/* x op y for the int32 lanes a and b; ne, le and ge are the complements
   of eq, gt and lt. */
RSIMD_ALWAYS_INLINE rsimd_mi32 rsimd_cmp_vi32(const int op, rsimd_vi32 a, rsimd_vi32 b) {
  switch (op) {
  case RSIMD_CMP_EQ: return rsimd_vi32_cmp_eq(a, b);
  case RSIMD_CMP_NE: return rsimd_mi32_not(rsimd_vi32_cmp_eq(a, b));
  case RSIMD_CMP_LT: return rsimd_vi32_cmp_lt(a, b);
  case RSIMD_CMP_LE: return rsimd_mi32_not(rsimd_vi32_cmp_gt(a, b));
  case RSIMD_CMP_GT: return rsimd_vi32_cmp_gt(a, b);
  default: return rsimd_mi32_not(rsimd_vi32_cmp_lt(a, b));
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cmp_i32_)(const int op, const int *x, const int *y,
                                                R_xlen_t n, int flags, int *out) {
  RSIMD_LANE32_LOOP(int32_t, rsimd_vi32_loadu, rsimd_vi32_loadu_p, int32_t, rsimd_vi32_storeu,
                    rsimd_vi32_storeu_p,
                    r = rsimd_lgl_vi32(rsimd_cmp_vi32(op, a, b), rsimd_vi32_na2(a, b)));
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cmp_u8_)(const int op, const Rbyte *x, const Rbyte *y,
                                               R_xlen_t n, int flags, int *out) {
  const rsimd_mi32 none = rsimd_mi32_none();
  RSIMD_LANE32_LOOP(uint8_t, rsimd_vi32_loadu_u8, rsimd_vi32_loadu_u8_p, int32_t,
                    rsimd_vi32_storeu, rsimd_vi32_storeu_p,
                    r = rsimd_lgl_vi32(rsimd_cmp_vi32(op, a, b), none));
}

/* Calls KERNEL(op, ...) with op as a constant, so that the switch in the
   lane comparison compiles away. */
#define RSIMD_CMP_DISPATCH(KERNEL, ...)                                                    \
  switch (op) {                                                                            \
  case RSIMD_CMP_EQ: KERNEL(RSIMD_CMP_EQ, __VA_ARGS__); break;                             \
  case RSIMD_CMP_NE: KERNEL(RSIMD_CMP_NE, __VA_ARGS__); break;                             \
  case RSIMD_CMP_LT: KERNEL(RSIMD_CMP_LT, __VA_ARGS__); break;                             \
  case RSIMD_CMP_LE: KERNEL(RSIMD_CMP_LE, __VA_ARGS__); break;                             \
  case RSIMD_CMP_GT: KERNEL(RSIMD_CMP_GT, __VA_ARGS__); break;                             \
  default: KERNEL(RSIMD_CMP_GE, __VA_ARGS__); break;                                       \
  }

void RSIMD_KERNEL(cmp_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                           int *out) {
  RSIMD_CMP_DISPATCH(RSIMD_KERNEL(cmp_i32_), x, y, n, flags, out)
}

void RSIMD_KERNEL(cmp_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                          int *out);
void RSIMD_KERNEL(cmp_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                          int *out) {
  RSIMD_CMP_DISPATCH(RSIMD_KERNEL(cmp_u8_), x, y, n, flags, out)
}

#ifndef RSIMD_NO_F64_SIMD

/* x op y for the double lanes a and b: false where either is NaN (the
   caller puts NA there), as the ordered comparisons are. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cmp_vf64(const int op, rsimd_vf64 a, rsimd_vf64 b) {
  switch (op) {
  case RSIMD_CMP_EQ: return rsimd_vf64_cmp_eq(a, b);
  case RSIMD_CMP_NE: return rsimd_vf64_cmp_ne(a, b);
  case RSIMD_CMP_LT: return rsimd_vf64_cmp_lt(a, b);
  case RSIMD_CMP_LE: return rsimd_vf64_cmp_le(a, b);
  case RSIMD_CMP_GT: return rsimd_vf64_cmp_gt(a, b);
  default: return rsimd_vf64_cmp_ge(a, b);
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cmp_f64_)(const int op, const void *x, const void *y,
                                                R_xlen_t n, int flags, int *out) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_ew_bcast(y, flags, 1, 1);
  RSIMD_LANE64_LOOP_FLAGS({
    rsimd_mf64 na = rsimd_mf64_or(rsimd_vf64_is_nan(a), rsimd_vf64_is_nan(b));
    r = rsimd_lgl_vi64(rsimd_cmp_vf64(op, a, b), na);
  });
}

void RSIMD_KERNEL(cmp_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out) {
  RSIMD_CMP_DISPATCH(RSIMD_KERNEL(cmp_f64_), x, y, n, flags, out)
}

#else /* RSIMD_NO_F64_SIMD */
#define RSIMD_SKIP_cmp_f64 1
#endif

#undef RSIMD_CMP_DISPATCH

#endif /* vector tiers */
