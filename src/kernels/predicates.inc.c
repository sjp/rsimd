/* Predicate kernels: is.na, is.nan, is.finite, is.infinite, negative and
 * zero of double and int32 elements (op codes RSIMD_PRED_* in
 * kernel_types.h), as a logical vector (mode RSIMD_PRED_ELT) or folded to
 * any/all (RSIMD_PRED_ANY, RSIMD_PRED_ALL), which return as soon as a
 * vector decides the answer.
 *
 * A double is NA when it is any NaN (base R's is.na), NaN when it is a NaN
 * that is not NA_real_, negative when its sign bit is set and it is not a
 * NaN (so -0 and -Inf are negative), zero when it equals 0 (so is -0). An
 * int32 element is NA when it is INT32_MIN; it is never NaN or infinite.
 */

#include "logical.inc.h"

/* The predicate of one element: the none tier's kernels and the
   reference of the others. */
static inline int rsimd_pred_f64_1(int op, double x) {
  switch (op) {
  case RSIMD_PRED_NA: return isnan(x) != 0;
  case RSIMD_PRED_NAN: return rsimd_is_nan_f64(x);
  case RSIMD_PRED_FINITE: return isfinite(x) != 0;
  case RSIMD_PRED_INFINITE: return isinf(x) != 0;
  case RSIMD_PRED_NEGATIVE: return !isnan(x) && signbit(x) != 0;
  default: return x == 0;
  }
}
static inline int rsimd_pred_i32_1(int op, int32_t x) {
  switch (op) {
  case RSIMD_PRED_NA: return x == RSIMD_NA_I32;
  case RSIMD_PRED_FINITE: return x != RSIMD_NA_I32;
  case RSIMD_PRED_NEGATIVE: return x < 0 && x != RSIMD_NA_I32;
  default: return x == 0;
  }
}

#if RSIMD_TIER_IS(none)

int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int v = rsimd_pred_f64_1(op, x[i]);
    if (mode == RSIMD_PRED_ELT) out[i] = v;
    else if (mode == RSIMD_PRED_ANY && v) return 1;
    else if (mode == RSIMD_PRED_ALL && !v) return 0;
  }
  return mode == RSIMD_PRED_ALL;
}

int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int v = rsimd_pred_i32_1(op, x[i]);
    if (mode == RSIMD_PRED_ELT) out[i] = v;
    else if (mode == RSIMD_PRED_ANY && v) return 1;
    else if (mode == RSIMD_PRED_ALL && !v) return 0;
  }
  return mode == RSIMD_PRED_ALL;
}

#else /* vector tiers */

/* Runs the predicate over the chunk: m is the mask of the vector loaded
   from x at i (by the statement `load`), with the elementwise result
   stored by `store`. In the last vector the inactive lanes repeat
   element i, so they cannot change an any/all answer. */
#define RSIMD_PRED_LOOP(LANES, MASK_ANY, MASK_ALL, load, load_p, store, store_p)            \
  do {                                                                                     \
    ptrdiff_t i = 0;                                                                       \
    for (; i + (LANES) <= n; i += (LANES)) {                                               \
      load;                                                                                \
      if (mode == RSIMD_PRED_ELT) {                                                        \
        store;                                                                             \
      } else if (mode == RSIMD_PRED_ANY) {                                                 \
        if (MASK_ANY(m)) return 1;                                                         \
      } else if (!MASK_ALL(m)) {                                                           \
        return 0;                                                                          \
      }                                                                                    \
    }                                                                                      \
    if (i < n) {                                                                           \
      load_p;                                                                              \
      if (mode == RSIMD_PRED_ELT) {                                                        \
        store_p;                                                                           \
      } else if (mode == RSIMD_PRED_ANY) {                                                 \
        return MASK_ANY(m);                                                                \
      } else {                                                                             \
        return MASK_ALL(m);                                                                \
      }                                                                                    \
    }                                                                                      \
    return mode == RSIMD_PRED_ALL;                                                         \
  } while (0)

/* ---- int32 ---------------------------------------------------------------- */

RSIMD_ALWAYS_INLINE rsimd_mi32 rsimd_pred_vi32(const int op, rsimd_vi32 v) {
  switch (op) {
  case RSIMD_PRED_NA: return rsimd_vi32_is_na(v);
  case RSIMD_PRED_FINITE: return rsimd_mi32_not(rsimd_vi32_is_na(v));
  case RSIMD_PRED_NEGATIVE:
    return rsimd_mi32_andnot(rsimd_vi32_is_na(v), rsimd_vi32_cmp_lt(v, rsimd_vi32_zero()));
  default: return rsimd_vi32_cmp_eq(v, rsimd_vi32_zero());
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_i32_)(const int op, const int32_t *x, R_xlen_t n,
                                                int mode, int32_t *out) {
  const rsimd_mi32 none = rsimd_mi32_none();
  (void) none;
  RSIMD_PRED_LOOP(RSIMD_LANES_32, rsimd_mi32_any, rsimd_mi32_all,
                  rsimd_mi32 m = rsimd_pred_vi32(op, rsimd_vi32_loadu(x + i)),
                  rsimd_p32 pg = rsimd_p32_while(i, n);
                  rsimd_mi32 m = rsimd_pred_vi32(op, rsimd_vi32_loadu_p(pg, x + i, x[i])),
                  rsimd_vi32_storeu(out + i, rsimd_lgl_vi32(m, none)),
                  rsimd_vi32_storeu_p(pg, out + i, rsimd_lgl_vi32(m, none)));
}

int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_FINITE: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_FINITE, x, n, mode, out);
  case RSIMD_PRED_NEGATIVE: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NEGATIVE, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_ZERO, x, n, mode, out);
  }
}

/* ---- double --------------------------------------------------------------- */

#ifndef RSIMD_NO_F64_SIMD

RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_pred_vf64(const int op, rsimd_vf64 v) {
  switch (op) {
  case RSIMD_PRED_NA: return rsimd_vf64_is_nan(v);
  case RSIMD_PRED_NAN: return rsimd_mf64_andnot(rsimd_vf64_is_na(v), rsimd_vf64_is_nan(v));
  /* |x| < Inf and |x| == Inf are both false for NaN. */
  case RSIMD_PRED_FINITE: return rsimd_vf64_cmp_lt(rsimd_vf64_abs(v), rsimd_vf64_set1(HUGE_VAL));
  case RSIMD_PRED_INFINITE: return rsimd_vf64_cmp_eq(rsimd_vf64_abs(v), rsimd_vf64_set1(HUGE_VAL));
  /* The sign bit, as the sign of the 64-bit integer view. */
  case RSIMD_PRED_NEGATIVE:
    return rsimd_mf64_andnot(rsimd_vf64_is_nan(v),
                             rsimd_mi64_to_mf64(rsimd_vi64_cmp_lt(rsimd_vf64_as_vi64(v),
                                                                  rsimd_vi64_zero())));
  default: return rsimd_vf64_cmp_eq(v, rsimd_vf64_zero());
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_f64_)(const int op, const double *x, R_xlen_t n,
                                                int mode, int32_t *out) {
  const rsimd_mf64 none = rsimd_mf64_none();
  (void) none;
  RSIMD_PRED_LOOP(RSIMD_LANES_64, rsimd_mf64_any, rsimd_mf64_all,
                  rsimd_mf64 m = rsimd_pred_vf64(op, rsimd_vf64_loadu(x + i)),
                  rsimd_p64 pg = rsimd_p64_while(i, n);
                  rsimd_mf64 m = rsimd_pred_vf64(op, rsimd_vf64_loadu_p(pg, x + i, x[i])),
                  rsimd_vi64_storeu_i32(out + i, rsimd_lgl_vi64(m, none)),
                  rsimd_vi64_storeu_i32_p(pg, out + i, rsimd_lgl_vi64(m, none)));
}

int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_NAN: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NAN, x, n, mode, out);
  case RSIMD_PRED_FINITE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_FINITE, x, n, mode, out);
  case RSIMD_PRED_INFINITE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_INFINITE, x, n, mode, out);
  case RSIMD_PRED_NEGATIVE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NEGATIVE, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_ZERO, x, n, mode, out);
  }
}

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does doubles */
#define RSIMD_SKIP_pred_f64 1
#endif

#undef RSIMD_PRED_LOOP

#endif /* vector tiers */
