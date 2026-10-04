/* Type conversion kernels: one slot, `convert`, with the conversion and
 * mode codes and status bits of kernel_types.h.
 *
 * To integer (from double): truncate toward zero, then
 *   CHECKED     NA for NaN and, with RSIMD_CVT_WARN_INT, for x >= 2^31 or
 *               x <= -2^31 (base R's IntegerFromReal: -2^31 is NA, while
 *               -2147483647.5 truncates to -2147483647);
 *   SATURATING  clamp to [-INT32_MAX, INT32_MAX]; NaN gives NA;
 *   TRUNCATING  the truncated value modulo 2^32 in two's complement (a
 *               result of -2^31 is NA); NaN and infinities give NA.
 * To raw (from double or int32): truncate, then
 *   CHECKED     00 with RSIMD_CVT_WARN_RAW for a missing value or one
 *               outside 0..255; doubles go through CHECKED integer first,
 *               as base R does, so a double outside the integer range also
 *               sets RSIMD_CVT_WARN_INT;
 *   SATURATING  clamp to 0..255; missing values give 00;
 *   TRUNCATING  the low 8 bits (modulo 256); missing values and
 *               infinities give 00.
 * To logical: non-zero is TRUE, 0 FALSE, a missing value NA.
 * To double: exact; NA_integer_ becomes NA_real_.
 * Status bits are only set in CHECKED mode. The integer64 conversions are
 * in int64.inc.c (included before this file), which the slot calls.
 */

#include "logical.inc.h"

/* ---- Scalar forms: the none tier, and references ------------------------- */

static inline int32_t rsimd_cvt_f64_i32_1(double x, int mode, int *st) {
  double t, w;
  if (isnan(x)) return RSIMD_NA_I32;
  t = trunc(x);
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    if (x >= 2147483648.0 || x <= -2147483648.0) {
      *st |= RSIMD_CVT_WARN_INT;
      return RSIMD_NA_I32;
    }
    return (int32_t) t;
  case RSIMD_CVT_SATURATING:
    if (t >= 2147483647.0) return INT32_MAX;
    if (t <= -2147483647.0) return -INT32_MAX;
    return (int32_t) t;
  default:
    if (isinf(x)) return RSIMD_NA_I32;
    w = fmod(t, 4294967296.0);
    if (w < 0) w += 4294967296.0;
    /* w is in [0, 2^32): -2^31 (NA) for w == 2^31. */
    return w >= 2147483648.0 ? (int32_t) (w - 4294967296.0) : (int32_t) w;
  }
}

static inline uint8_t rsimd_cvt_i32_u8_1(int32_t x, int mode, int *st) {
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    if (x == RSIMD_NA_I32 || x < 0 || x > 255) {
      *st |= RSIMD_CVT_WARN_RAW;
      return 0;
    }
    return (uint8_t) x;
  case RSIMD_CVT_SATURATING:
    if (x == RSIMD_NA_I32 || x < 0) return 0;
    return x > 255 ? 255 : (uint8_t) x;
  default: return x == RSIMD_NA_I32 ? 0 : (uint8_t) ((uint32_t) x & 0xFFu);
  }
}

static inline uint8_t rsimd_cvt_f64_u8_1(double x, int mode, int *st) {
  double t, w;
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    return rsimd_cvt_i32_u8_1(rsimd_cvt_f64_i32_1(x, RSIMD_CVT_CHECKED, st), RSIMD_CVT_CHECKED,
                              st);
  case RSIMD_CVT_SATURATING:
    if (isnan(x)) return 0;
    t = trunc(x);
    return t <= 0 ? 0 : t >= 255 ? 255 : (uint8_t) t;
  default:
    if (!isfinite(x)) return 0;
    w = fmod(trunc(x), 256.0);
    if (w < 0) w += 256.0;
    return (uint8_t) w;
  }
}

static inline int rsimd_cvt_1(int op, int mode, const void *x, R_xlen_t i, void *out) {
  int st = 0;
  switch (op) {
  case RSIMD_CVT_F64_I32:
    ((int *) out)[i] = rsimd_cvt_f64_i32_1(((const double *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_F64_U8:
    ((Rbyte *) out)[i] = rsimd_cvt_f64_u8_1(((const double *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_F64_LGL: {
    double v = ((const double *) x)[i];
    ((int *) out)[i] = isnan(v) ? RSIMD_NA_I32 : v != 0;
    break;
  }
  case RSIMD_CVT_I32_F64: {
    int v = ((const int *) x)[i];
    ((double *) out)[i] = v == RSIMD_NA_I32 ? rsimd_na_real() : (double) v;
    break;
  }
  case RSIMD_CVT_I32_U8:
    ((Rbyte *) out)[i] = rsimd_cvt_i32_u8_1(((const int *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_I32_LGL: {
    int v = ((const int *) x)[i];
    ((int *) out)[i] = v == RSIMD_NA_I32 ? v : v != 0;
    break;
  }
  case RSIMD_CVT_U8_F64: ((double *) out)[i] = ((const Rbyte *) x)[i]; break;
  case RSIMD_CVT_U8_I32: ((int *) out)[i] = ((const Rbyte *) x)[i]; break;
  default: ((int *) out)[i] = ((const Rbyte *) x)[i] != 0; break;
  }
  return st;
}

#if RSIMD_TIER_IS(none)

int RSIMD_KERNEL(convert)(int op, int mode, const void *x, R_xlen_t n, void *out);
int RSIMD_KERNEL(convert)(int op, int mode, const void *x, R_xlen_t n, void *out) {
  int st = 0;
  R_xlen_t i;
  if (op >= RSIMD_CVT_I64_F64) return RSIMD_KERNEL(convert_i64_)(op, mode, x, n, out);
  for (i = 0; i < n; i++) st |= rsimd_cvt_1(op, mode, x, i, out);
  return st;
}

#elif !defined(RSIMD_NO_F64_SIMD) /* vector tiers with double vectors */

/* Two-step conversions go through a stack buffer of int32 elements. */
#define RSIMD_CVT_BLOCK 1024

/* x truncated and converted to an integer value in double lanes (stored
   with rsimd_vf64_storeu_i32, so every lane must be in the int32 range;
   -2^31 is NA). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_cvt_vf64_i32(rsimd_vf64 x, const int mode, int *st) {
  const rsimd_vf64 t = rsimd_vf64_trunc(x), na = rsimd_vf64_set1(RSIMD_NA_I32_AS_F64);
  const rsimd_mf64 nan = rsimd_vf64_is_nan(x);
  if (mode == RSIMD_CVT_CHECKED) {
    rsimd_mf64 oob = rsimd_mf64_or(rsimd_vf64_cmp_ge(x, rsimd_vf64_set1(2147483648.0)),
                                   rsimd_vf64_cmp_le(x, rsimd_vf64_set1(-2147483648.0)));
    if (rsimd_mf64_any(oob)) *st |= RSIMD_CVT_WARN_INT;
    return rsimd_vf64_blend(t, na, rsimd_mf64_or(nan, oob));
  }
  if (mode == RSIMD_CVT_SATURATING) {
    /* max() gives its second operand for NaN; NaN lanes are reset anyway. */
    rsimd_vf64 r = rsimd_vf64_min(rsimd_vf64_max(t, rsimd_vf64_set1(-2147483647.0)),
                                  rsimd_vf64_set1(2147483647.0));
    return rsimd_vf64_blend(r, na, nan);
  }
  {
    /* t - 2^32 floor(t / 2^32) is exact (t is an integer), in [0, 2^32). */
    const rsimd_vf64 two32 = rsimd_vf64_set1(4294967296.0);
    rsimd_vf64 w = rsimd_vf64_sub(
      t, rsimd_vf64_mul(two32, rsimd_vf64_floor(rsimd_vf64_mul(t, rsimd_vf64_set1(0x1p-32)))));
    w = rsimd_vf64_blend(w, rsimd_vf64_sub(w, two32),
                         rsimd_vf64_cmp_ge(w, rsimd_vf64_set1(2147483648.0)));
    return rsimd_vf64_blend(
      w, na, rsimd_mf64_not(rsimd_vf64_cmp_lt(rsimd_vf64_abs(x), rsimd_vf64_set1(HUGE_VAL))));
  }
}

/* x truncated to a raw value 0..255 in double lanes. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_cvt_vf64_u8(rsimd_vf64 x, const int mode, int *st) {
  const rsimd_vf64 t = rsimd_vf64_trunc(x), zero = rsimd_vf64_zero(),
                   top = rsimd_vf64_set1(255.0);
  if (mode == RSIMD_CVT_CHECKED) {
    /* -0.5 truncates to -0, which is not < 0. */
    rsimd_mf64 oob = rsimd_mf64_or(rsimd_vf64_cmp_ge(x, rsimd_vf64_set1(2147483648.0)),
                                   rsimd_vf64_cmp_le(x, rsimd_vf64_set1(-2147483648.0)));
    rsimd_mf64 bad = rsimd_mf64_or(rsimd_vf64_cmp_lt(t, zero), rsimd_vf64_cmp_gt(t, top));
    bad = rsimd_mf64_or(rsimd_vf64_is_nan(x), bad);
    if (rsimd_mf64_any(oob)) *st |= RSIMD_CVT_WARN_INT;
    if (rsimd_mf64_any(bad)) *st |= RSIMD_CVT_WARN_RAW;
    return rsimd_vf64_blend(t, zero, bad);
  }
  if (mode == RSIMD_CVT_SATURATING) {
    /* max(NaN, 0) is 0. */
    return rsimd_vf64_min(rsimd_vf64_max(t, zero), top);
  }
  {
    const rsimd_vf64 b = rsimd_vf64_set1(256.0);
    rsimd_vf64 w = rsimd_vf64_sub(
      t, rsimd_vf64_mul(b, rsimd_vf64_floor(rsimd_vf64_mul(t, rsimd_vf64_set1(0x1p-8)))));
    return rsimd_vf64_blend(
      w, zero, rsimd_mf64_not(rsimd_vf64_cmp_lt(rsimd_vf64_abs(x), rsimd_vf64_set1(HUGE_VAL))));
  }
}

/* x to int32 with fn (rsimd_cvt_vf64_i32 or rsimd_cvt_vf64_u8). */
#define RSIMD_CVT_F64_LOOP(fn)                                                             \
  do {                                                                                     \
    ptrdiff_t i = 0;                                                                       \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                                 \
      rsimd_vf64_storeu_i32(o32 + i, fn(rsimd_vf64_loadu(x64 + i), mode, &st));            \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                                \
      rsimd_vf64 v_ = rsimd_vf64_loadu_p(pg, x64 + i, x64[i]);                            \
      rsimd_vf64_storeu_i32_p(pg, o32 + i, fn(v_, mode, &st));                             \
    }                                                                                      \
  } while (0)

/* Elements of the int32 buffer b (values 0..255) to bytes. */
RSIMD_INLINE void rsimd_cvt_i32_to_u8(const int32_t *b, ptrdiff_t n, uint8_t *out) {
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {
    rsimd_vi32_storeu_u8(out + i, rsimd_vi32_loadu(b + i));
  }
  if (i < n) {
    rsimd_p32 pg = rsimd_p32_while(i, n);
    rsimd_vi32_storeu_u8_p(pg, out + i, rsimd_vi32_loadu_p(pg, b + i, 0));
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(f64_to_i32_)(const double *x64, R_xlen_t n, const int mode,
                                                  int32_t *o32) {
  int st = 0;
  RSIMD_CVT_F64_LOOP(rsimd_cvt_vf64_i32);
  return st;
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(f64_to_u8_)(const double *x, R_xlen_t len, const int mode,
                                                 uint8_t *out) {
  int32_t buf[RSIMD_CVT_BLOCK];
  int st = 0;
  R_xlen_t off;
  for (off = 0; off < len; off += RSIMD_CVT_BLOCK) {
    const R_xlen_t n = len - off < RSIMD_CVT_BLOCK ? len - off : RSIMD_CVT_BLOCK;
    const double *x64 = x + off;
    int32_t *o32 = buf;
    RSIMD_CVT_F64_LOOP(rsimd_cvt_vf64_u8);
    rsimd_cvt_i32_to_u8(buf, n, out + off);
  }
  return st;
}

/* Logical lanes from the TRUE and NA masks of the elements. */
static inline void RSIMD_KERNEL(f64_to_lgl_)(const double *x, R_xlen_t n, int32_t *out) {
  ptrdiff_t i = 0;
#define RSIMD_CVT_LGL64(v)                                                                 \
  rsimd_lgl_vi64(rsimd_mf64_not(rsimd_mf64_or(rsimd_vf64_is_nan(v),                         \
                                              rsimd_vf64_cmp_eq(v, rsimd_vf64_zero()))),    \
                 rsimd_vf64_is_nan(v))
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vf64 v = rsimd_vf64_loadu(x + i);
    rsimd_vi64_storeu_i32(out + i, RSIMD_CVT_LGL64(v));
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64 v = rsimd_vf64_loadu_p(pg, x + i, 0.0);
    rsimd_vi64_storeu_i32_p(pg, out + i, RSIMD_CVT_LGL64(v));
  }
#undef RSIMD_CVT_LGL64
}

static inline void RSIMD_KERNEL(i32_to_f64_)(const int32_t *x, R_xlen_t n, double *out) {
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vf64_storeu(out + i, rsimd_ew_from_i32(rsimd_vf64_loadu_i32(x + i), 1));
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64_storeu_p(pg, out + i, rsimd_ew_from_i32(rsimd_vf64_loadu_i32_p(pg, x + i, 0), 1));
  }
}

/* int32 (integer or logical) to raw. */
RSIMD_ALWAYS_INLINE rsimd_vi32 rsimd_cvt_vi32_u8(rsimd_vi32 v, const int mode, int *st) {
  if (mode == RSIMD_CVT_CHECKED) {
    /* NA (INT32_MIN) is < 0. */
    rsimd_mi32 bad = rsimd_mi32_or(rsimd_vi32_cmp_lt(v, rsimd_vi32_zero()),
                                   rsimd_vi32_cmp_gt(v, rsimd_vi32_set1(255)));
    if (rsimd_mi32_any(bad)) *st |= RSIMD_CVT_WARN_RAW;
    return rsimd_vi32_blend(v, rsimd_vi32_zero(), bad);
  }
  if (mode == RSIMD_CVT_SATURATING) {
    return rsimd_vi32_min(rsimd_vi32_max(v, rsimd_vi32_zero()), rsimd_vi32_set1(255));
  }
  /* The low byte of NA (0x80000000) is 0. */
  return rsimd_vi32_and(v, rsimd_vi32_set1(0xFF));
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(i32_to_u8_)(const int32_t *x, R_xlen_t n, const int mode,
                                                 uint8_t *out) {
  int st = 0;
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {
    rsimd_vi32_storeu_u8(out + i, rsimd_cvt_vi32_u8(rsimd_vi32_loadu(x + i), mode, &st));
  }
  if (i < n) {
    rsimd_p32 pg = rsimd_p32_while(i, n);
    rsimd_vi32 v = rsimd_vi32_loadu_p(pg, x + i, x[i]);
    rsimd_vi32_storeu_u8_p(pg, out + i, rsimd_cvt_vi32_u8(v, mode, &st));
  }
  return st;
}

/* Runs the statements in `...`, which set the rsimd_vi32 r from the lanes v (loaded with
   LD / LD_P from x of element type T), storing r as int32 elements. */
#define RSIMD_CVT_I32_LOOP(T, LD, LD_P, ...)                                               \
  do {                                                                                     \
    const T *xp_ = (const T *) x;                                                          \
    int32_t *op_ = (int32_t *) out;                                                        \
    ptrdiff_t i = 0;                                                                       \
    for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {                                 \
      rsimd_vi32 v = LD(xp_ + i), r;                                                       \
      __VA_ARGS__;                                                                          \
      rsimd_vi32_storeu(op_ + i, r);                                                       \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p32 pg = rsimd_p32_while(i, n);                                                \
      rsimd_vi32 v = LD_P(pg, xp_ + i, 0), r;                                              \
      __VA_ARGS__;                                                                          \
      rsimd_vi32_storeu_p(pg, op_ + i, r);                                                 \
    }                                                                                      \
  } while (0)

static inline void RSIMD_KERNEL(u8_to_f64_)(const uint8_t *x, R_xlen_t len, double *out) {
  int32_t buf[RSIMD_CVT_BLOCK];
  R_xlen_t off;
  for (off = 0; off < len; off += RSIMD_CVT_BLOCK) {
    const R_xlen_t n = len - off < RSIMD_CVT_BLOCK ? len - off : RSIMD_CVT_BLOCK;
    ptrdiff_t i = 0;
    for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {
      rsimd_vi32_storeu(buf + i, rsimd_vi32_loadu_u8(x + off + i));
    }
    if (i < n) {
      rsimd_p32 pg = rsimd_p32_while(i, n);
      rsimd_vi32_storeu_p(pg, buf + i, rsimd_vi32_loadu_u8_p(pg, x + off + i, 0));
    }
    for (i = 0; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
      rsimd_vf64_storeu(out + off + i, rsimd_vf64_loadu_i32(buf + i));
    }
    if (i < n) {
      rsimd_p64 pg = rsimd_p64_while(i, n);
      rsimd_vf64_storeu_p(pg, out + off + i, rsimd_vf64_loadu_i32_p(pg, buf + i, 0));
    }
  }
}

int RSIMD_KERNEL(convert)(int op, int mode, const void *x, R_xlen_t n, void *out);
int RSIMD_KERNEL(convert)(int op, int mode, const void *x, R_xlen_t n, void *out) {
  const rsimd_mi32 none = rsimd_mi32_none();
  (void) none;
  if (op >= RSIMD_CVT_I64_F64) return RSIMD_KERNEL(convert_i64_)(op, mode, x, n, out);
  switch (op) {
  case RSIMD_CVT_F64_I32:
    switch (mode) {
    case RSIMD_CVT_CHECKED: return RSIMD_KERNEL(f64_to_i32_)(x, n, RSIMD_CVT_CHECKED, out);
    case RSIMD_CVT_SATURATING: return RSIMD_KERNEL(f64_to_i32_)(x, n, RSIMD_CVT_SATURATING, out);
    default: return RSIMD_KERNEL(f64_to_i32_)(x, n, RSIMD_CVT_TRUNCATING, out);
    }
  case RSIMD_CVT_F64_U8:
    switch (mode) {
    case RSIMD_CVT_CHECKED: return RSIMD_KERNEL(f64_to_u8_)(x, n, RSIMD_CVT_CHECKED, out);
    case RSIMD_CVT_SATURATING: return RSIMD_KERNEL(f64_to_u8_)(x, n, RSIMD_CVT_SATURATING, out);
    default: return RSIMD_KERNEL(f64_to_u8_)(x, n, RSIMD_CVT_TRUNCATING, out);
    }
  case RSIMD_CVT_F64_LGL: RSIMD_KERNEL(f64_to_lgl_)(x, n, out); return 0;
  case RSIMD_CVT_I32_F64: RSIMD_KERNEL(i32_to_f64_)(x, n, out); return 0;
  case RSIMD_CVT_I32_U8:
    switch (mode) {
    case RSIMD_CVT_CHECKED: return RSIMD_KERNEL(i32_to_u8_)(x, n, RSIMD_CVT_CHECKED, out);
    case RSIMD_CVT_SATURATING: return RSIMD_KERNEL(i32_to_u8_)(x, n, RSIMD_CVT_SATURATING, out);
    default: return RSIMD_KERNEL(i32_to_u8_)(x, n, RSIMD_CVT_TRUNCATING, out);
    }
  case RSIMD_CVT_I32_LGL:
    RSIMD_CVT_I32_LOOP(int32_t, rsimd_vi32_loadu, rsimd_vi32_loadu_p, {
      rsimd_mi32 na = rsimd_vi32_is_na(v);
      r = rsimd_lgl_vi32(rsimd_mi32_not(rsimd_mi32_or(na, rsimd_vi32_cmp_eq(v, rsimd_vi32_zero()))),
                         na);
    });
    return 0;
  case RSIMD_CVT_U8_F64: RSIMD_KERNEL(u8_to_f64_)(x, n, out); return 0;
  case RSIMD_CVT_U8_I32:
    RSIMD_CVT_I32_LOOP(uint8_t, rsimd_vi32_loadu_u8, rsimd_vi32_loadu_u8_p, r = v);
    return 0;
  default:
    RSIMD_CVT_I32_LOOP(uint8_t, rsimd_vi32_loadu_u8, rsimd_vi32_loadu_u8_p,
                       r = rsimd_lgl_vi32(rsimd_mi32_not(rsimd_vi32_cmp_eq(v, rsimd_vi32_zero())),
                                          none));
    return 0;
  }
}

#undef RSIMD_CVT_F64_LOOP
#undef RSIMD_CVT_I32_LOOP
#undef RSIMD_CVT_BLOCK

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier converts */
#define RSIMD_SKIP_convert 1
#endif
