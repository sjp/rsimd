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
 *
 * The number classes (normal, subnormal, whole, even, odd, pow2) are
 * false for NA and NaN. The vector tiers test doubles through their bit
 * patterns (exponent field 0 for zero and subnormals, 0x7FF for infinities
 * and NaN) and trunc, which is exact on every tier: x is whole when
 * trunc(x) == x and x is finite, and a whole x is even when x / 2 is whole
 * (x * 0.5 is exact for whole x). A positive finite double is a power of
 * two when its mantissa field is 0 (normal) or has a single bit set
 * (subnormal, exponent field 0).
 */

#include "logical.inc.h"

/* The predicate of one element: the none tier's kernels and the
   reference of the others. */
static inline int rsimd_whole_f64_1(double x) { return isfinite(x) && trunc(x) == x; }
static inline int rsimd_pred_f64_1(int op, double x) {
  int e;
  switch (op) {
  case RSIMD_PRED_NA: return isnan(x) != 0;
  case RSIMD_PRED_NAN: return rsimd_is_nan_f64(x);
  case RSIMD_PRED_FINITE: return isfinite(x) != 0;
  case RSIMD_PRED_INFINITE: return isinf(x) != 0;
  case RSIMD_PRED_NEGATIVE: return !isnan(x) && signbit(x) != 0;
  case RSIMD_PRED_NORMAL: return isnormal(x) != 0;
  case RSIMD_PRED_SUBNORMAL: return fpclassify(x) == FP_SUBNORMAL;
  case RSIMD_PRED_WHOLE: return rsimd_whole_f64_1(x);
  case RSIMD_PRED_EVEN: return rsimd_whole_f64_1(x) && fmod(x, 2.0) == 0;
  case RSIMD_PRED_ODD: return rsimd_whole_f64_1(x) && fabs(fmod(x, 2.0)) == 1;
  case RSIMD_PRED_POW2: return x > 0 && isfinite(x) && frexp(x, &e) == 0.5;
  default: return x == 0;
  }
}
static inline int rsimd_pred_i32_1(int op, int32_t x) {
  switch (op) {
  case RSIMD_PRED_NA: return x == RSIMD_NA_I32;
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return x != RSIMD_NA_I32;
  case RSIMD_PRED_NEGATIVE: return x < 0 && x != RSIMD_NA_I32;
  case RSIMD_PRED_NORMAL: return x != 0 && x != RSIMD_NA_I32;
  case RSIMD_PRED_EVEN: return x % 2 == 0 && x != RSIMD_NA_I32;
  case RSIMD_PRED_ODD: return x % 2 != 0;
  case RSIMD_PRED_POW2: return x > 0 && (x & (x - 1)) == 0;
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return 0;
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

/* Runs the predicate PRED(op, v) over the chunk of x, read with LD / LD_P,
   giving masks of type MASK (rsimd_mi32 or rsimd_mf64) for vectors of
   W-bit lanes; the elementwise result LGL(m, none) is stored with ST /
   ST_P. The any and all modes combine the masks of four vectors before
   each horizontal test. In the last vector the inactive lanes repeat
   element i, so they cannot change an any/all answer. */
#define RSIMD_PRED_LOOP(W, MASK, LD, LD_P, PRED, ST, ST_P, LGL)                              \
  do {                                                                                     \
    const ptrdiff_t l_ = RSIMD_LANES_##W;                                                  \
    ptrdiff_t i = 0;                                                                       \
    if (mode != RSIMD_PRED_ELT) {                                                          \
      for (; i + 4 * l_ <= n; i += 4 * l_) {                                               \
        MASK m0_ = PRED(op, LD(x + i)), m1_ = PRED(op, LD(x + i + l_));                    \
        MASK m2_ = PRED(op, LD(x + i + 2 * l_)), m3_ = PRED(op, LD(x + i + 3 * l_));       \
        if (mode == RSIMD_PRED_ANY) {                                                      \
          if (MASK##_any(MASK##_or(MASK##_or(m0_, m1_), MASK##_or(m2_, m3_)))) return 1;   \
        } else if (!MASK##_all(MASK##_and(MASK##_and(m0_, m1_), MASK##_and(m2_, m3_)))) { \
          return 0;                                                                        \
        }                                                                                  \
      }                                                                                    \
    }                                                                                      \
    RSIMD_CHUNK_LOOP(W, i, n, {                                                            \
      MASK m = PRED(op, RSIMD_LDT(LD, LD_P, x, i));                                        \
      if (mode == RSIMD_PRED_ELT) {                                                        \
        RSIMD_STT(ST, ST_P, out + i, LGL(m, none));                                        \
      } else if (mode == RSIMD_PRED_ANY) {                                                 \
        if (MASK##_any(m)) return 1;                                                       \
      } else if (!MASK##_all(m)) {                                                         \
        return 0;                                                                          \
      }                                                                                    \
    });                                                                                    \
    return mode == RSIMD_PRED_ALL;                                                         \
  } while (0)

/* ---- int32 ---------------------------------------------------------------- */

RSIMD_ALWAYS_INLINE rsimd_mi32 rsimd_pred_vi32(const int op, rsimd_vi32 v) {
  const rsimd_vi32 zero = rsimd_vi32_zero(), one = rsimd_vi32_set1(1);
  switch (op) {
  case RSIMD_PRED_NA: return rsimd_vi32_is_na(v);
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return rsimd_mi32_not(rsimd_vi32_is_na(v));
  case RSIMD_PRED_NEGATIVE:
    return rsimd_mi32_andnot(rsimd_vi32_is_na(v), rsimd_vi32_cmp_lt(v, zero));
  case RSIMD_PRED_NORMAL:
    return rsimd_mi32_not(rsimd_mi32_or(rsimd_vi32_is_na(v), rsimd_vi32_cmp_eq(v, zero)));
  /* NA (INT32_MIN) has the low bit of an even number. */
  case RSIMD_PRED_EVEN:
    return rsimd_mi32_andnot(rsimd_vi32_is_na(v),
                             rsimd_vi32_cmp_eq(rsimd_vi32_and(v, one), zero));
  case RSIMD_PRED_ODD: return rsimd_vi32_cmp_eq(rsimd_vi32_and(v, one), one);
  case RSIMD_PRED_POW2:
    return rsimd_mi32_and(rsimd_vi32_cmp_gt(v, zero),
                          rsimd_vi32_cmp_eq(rsimd_vi32_and(v, rsimd_vi32_sub(v, one)), zero));
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return rsimd_mi32_none();
  default: return rsimd_vi32_cmp_eq(v, zero);
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_i32_)(const int op, const int32_t *x, R_xlen_t n,
                                                int mode, int32_t *out) {
  const rsimd_mi32 none = rsimd_mi32_none();
  (void) none;
  RSIMD_PRED_LOOP(32, rsimd_mi32, rsimd_vi32_loadu, rsimd_vi32_loadu_p, rsimd_pred_vi32,
                  rsimd_vi32_storeu, rsimd_vi32_storeu_p, rsimd_lgl_vi32);
}

int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_i32)(int op, const int *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_FINITE, x, n, mode, out);
  case RSIMD_PRED_NEGATIVE: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NEGATIVE, x, n, mode, out);
  case RSIMD_PRED_NORMAL: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NORMAL, x, n, mode, out);
  case RSIMD_PRED_EVEN: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_EVEN, x, n, mode, out);
  case RSIMD_PRED_ODD: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_ODD, x, n, mode, out);
  case RSIMD_PRED_POW2: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_POW2, x, n, mode, out);
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_NAN, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_i32_)(RSIMD_PRED_ZERO, x, n, mode, out);
  }
}

/* ---- double --------------------------------------------------------------- */

#ifndef RSIMD_NO_F64_SIMD

/* Fields of the bit pattern of a double. */
#define RSIMD_F64_EXP_BITS INT64_C(0x7FF0000000000000)
#define RSIMD_F64_MANT_BITS INT64_C(0x000FFFFFFFFFFFFF)

/* x is whole: trunc(x) == x, and x is finite (trunc(Inf) is Inf). */
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_whole(rsimd_vf64 v) {
  return rsimd_mf64_and(rsimd_vf64_cmp_eq(rsimd_vf64_trunc(v), v),
                        rsimd_vf64_cmp_lt(rsimd_vf64_abs(v), rsimd_vf64_set1(HUGE_VAL)));
}

/* A positive finite double (bit pattern in (0, 0x7FF0000000000000) as a
   signed integer) is 2^k when the mantissa field of a normal number is 0
   or the subnormal pattern (exponent field 0) has a single bit. */
RSIMD_INLINE rsimd_mi64 rsimd_vf64_is_pow2(rsimd_vf64 v) {
  const rsimd_vi64 b = rsimd_vf64_as_vi64(v), zero = rsimd_vi64_zero();
  const rsimd_vi64 ebits = rsimd_vi64_set1(RSIMD_F64_EXP_BITS);
  rsimd_mi64 sub = rsimd_vi64_cmp_eq(rsimd_vi64_and(b, ebits), zero);
  rsimd_vi64 t = rsimd_vi64_blend(rsimd_vi64_and(b, rsimd_vi64_set1(RSIMD_F64_MANT_BITS)),
                                  rsimd_vi64_and(b, rsimd_vi64_sub(b, rsimd_vi64_set1(1))), sub);
  return rsimd_mi64_and(rsimd_mi64_and(rsimd_vi64_cmp_gt(b, zero), rsimd_vi64_cmp_lt(b, ebits)),
                        rsimd_vi64_cmp_eq(t, zero));
}

RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_pred_vf64(const int op, rsimd_vf64 v) {
  const rsimd_vi64 ebits = rsimd_vi64_set1(RSIMD_F64_EXP_BITS);
  const rsimd_vi64 e = rsimd_vi64_and(rsimd_vf64_as_vi64(v), ebits);
  rsimd_vf64 h;
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
  case RSIMD_PRED_NORMAL:
    return rsimd_mi64_to_mf64(rsimd_mi64_not(
        rsimd_mi64_or(rsimd_vi64_cmp_eq(e, rsimd_vi64_zero()), rsimd_vi64_cmp_eq(e, ebits))));
  /* Exponent field 0 and not +-0. */
  case RSIMD_PRED_SUBNORMAL:
    return rsimd_mf64_andnot(rsimd_vf64_cmp_eq(v, rsimd_vf64_zero()),
                             rsimd_mi64_to_mf64(rsimd_vi64_cmp_eq(e, rsimd_vi64_zero())));
  case RSIMD_PRED_WHOLE: return rsimd_vf64_is_whole(v);
  case RSIMD_PRED_EVEN:
    h = rsimd_vf64_mul(v, rsimd_vf64_set1(0.5));
    return rsimd_mf64_and(rsimd_vf64_is_whole(v), rsimd_vf64_cmp_eq(rsimd_vf64_trunc(h), h));
  case RSIMD_PRED_ODD:
    h = rsimd_vf64_mul(v, rsimd_vf64_set1(0.5));
    return rsimd_mf64_andnot(rsimd_vf64_cmp_eq(rsimd_vf64_trunc(h), h), rsimd_vf64_is_whole(v));
  case RSIMD_PRED_POW2: return rsimd_mi64_to_mf64(rsimd_vf64_is_pow2(v));
  default: return rsimd_vf64_cmp_eq(v, rsimd_vf64_zero());
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_f64_)(const int op, const double *x, R_xlen_t n,
                                                int mode, int32_t *out) {
  const rsimd_mf64 none = rsimd_mf64_none();
  (void) none;
  RSIMD_PRED_LOOP(64, rsimd_mf64, rsimd_vf64_loadu, rsimd_vf64_loadu_p, rsimd_pred_vf64,
                  rsimd_vi64_storeu_i32, rsimd_vi64_storeu_i32_p, rsimd_lgl_vi64);
}

int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_f64)(int op, const double *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_NAN: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NAN, x, n, mode, out);
  case RSIMD_PRED_FINITE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_FINITE, x, n, mode, out);
  case RSIMD_PRED_INFINITE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_INFINITE, x, n, mode, out);
  case RSIMD_PRED_NEGATIVE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NEGATIVE, x, n, mode, out);
  case RSIMD_PRED_NORMAL: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_NORMAL, x, n, mode, out);
  case RSIMD_PRED_SUBNORMAL: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_SUBNORMAL, x, n, mode, out);
  case RSIMD_PRED_WHOLE: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_WHOLE, x, n, mode, out);
  case RSIMD_PRED_EVEN: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_EVEN, x, n, mode, out);
  case RSIMD_PRED_ODD: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_ODD, x, n, mode, out);
  case RSIMD_PRED_POW2: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_POW2, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_f64_)(RSIMD_PRED_ZERO, x, n, mode, out);
  }
}

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does doubles */
#define RSIMD_SKIP_pred_f64 1
#endif

#undef RSIMD_PRED_LOOP

#endif /* vector tiers */
