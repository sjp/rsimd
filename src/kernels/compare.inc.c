/* Comparison kernels: x op y (op codes RSIMD_CMP_* in kernel_types.h) as
 * logical vectors with base R's missing values: NA where either operand is
 * missing (any NaN for doubles, INT32_MIN for int32 elements). Doubles
 * compare as IEEE numbers, so -0 == 0. Bytes (raw) compare as unsigned
 * values and are never missing.
 *
 * cmp_f64 reads each operand as doubles or int32 elements (flag
 * RSIMD_EW_I32(k)), converting int32 exactly with NA becoming NA_real_.
 * RSIMD_EW_SCALAR(k) broadcasts element 0 of operand k.
 *
 * The Hamming kernels count the pairs that differ instead of storing the
 * comparisons. The vector tiers count in lane counters (rsimd_vi*_inc),
 * folded into r every RSIMD_HAMMING_BLOCK vectors, and do the last
 * partial vector with the scalar code. A double pair differs where
 * cmp_ne is true, which includes every pair with a NaN, so the count of
 * differing non-missing pairs is the count of cmp_ne lanes minus the
 * count of missing pairs.
 */

#include "logical.inc.h"

/* Element i of operand k of hamming_f64 as a double (int32 elements per
   RSIMD_EW_I32(k), NA_integer_ becoming NA_real_). */
static inline double rsimd_ham_get(const void *p, int flags, int k, R_xlen_t i) {
  R_xlen_t j = (flags & RSIMD_EW_SCALAR(k)) ? 0 : i;
  if (flags & RSIMD_EW_I32(k)) {
    int v = ((const int *) p)[j];
    return v == RSIMD_NA_I32 ? NAN : (double) v;
  }
  return ((const double *) p)[j];
}

/* The scalar Hamming loops over elements [i, n): the none tier's kernels
   and the vector tiers' tails. They return 1 when they stop at a missing
   pair (without na_rm). */
static inline int rsimd_ham_f64_from(const void *x, const void *y, R_xlen_t i, R_xlen_t n,
                                     int flags, rsimd_reduce_result *r, int na_rm) {
  for (; i < n; i++) {
    double a = rsimd_ham_get(x, flags, 0, i), b = rsimd_ham_get(y, flags, 1, i);
    if (isnan(a) || isnan(b)) {
      r->saw_na = 1;
      if (!na_rm) return 1;
      continue;
    }
    r->i64 += a != b;
  }
  return 0;
}
static inline int rsimd_ham_i32_from(const int *x, const int *y, R_xlen_t i, R_xlen_t n,
                                     int flags, rsimd_reduce_result *r, int na_rm) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  for (; i < n; i++) {
    int a = x[i * sx], b = y[i * sy];
    if (a == RSIMD_NA_I32 || b == RSIMD_NA_I32) {
      r->saw_na = 1;
      if (!na_rm) return 1;
      continue;
    }
    r->i64 += a != b;
  }
  return 0;
}
static inline void rsimd_ham_u8_from(const Rbyte *x, const Rbyte *y, R_xlen_t i, R_xlen_t n,
                                     int flags, rsimd_reduce_result *r) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  for (; i < n; i++) r->i64 += x[i * sx] != y[i * sy];
}

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

void RSIMD_KERNEL(hamming_f64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_f64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_ham_f64_from(x, y, 0, n, flags, r, o->na_rm);
}

void RSIMD_KERNEL(hamming_i32)(const int *x, const int *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_i32)(const int *x, const int *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_ham_i32_from(x, y, 0, n, flags, r, o->na_rm);
}

void RSIMD_KERNEL(hamming_u8)(const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                              rsimd_reduce_result *r);
void RSIMD_KERNEL(hamming_u8)(const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                              rsimd_reduce_result *r) {
  rsimd_ham_u8_from(x, y, 0, n, flags, r);
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

/* ---- Hamming distance ----------------------------------------------------- */

/* The int32 and byte kernels: d counts the differing non-missing pairs
   and m the missing ones (int32 only) in 32-bit lanes. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(hamming_i32_)(const int32_t *x, const int32_t *y, R_xlen_t n,
                                                    int flags, rsimd_reduce_result *r, int na_rm,
                                                    const int bytes) {
  const int sx = RSIMD_LGL_SCALAR(0), sy = RSIMD_LGL_SCALAR(1);
  const Rbyte *bx = (const Rbyte *) x, *by = (const Rbyte *) y;
  const rsimd_vi32 vx = rsimd_vi32_set1(bytes ? bx[0] : x[0]),
                   vy = rsimd_vi32_set1(bytes ? by[0] : y[0]);
  ptrdiff_t i = 0;
  while (i + RSIMD_LANES_32 <= n) {
    const ptrdiff_t end = (n - i) / RSIMD_LANES_32 > RSIMD_HAMMING_BLOCK
                              ? i + RSIMD_HAMMING_BLOCK * RSIMD_LANES_32
                              : n;
    rsimd_vi32 d = rsimd_vi32_zero(), m = rsimd_vi32_zero();
    int64_t missing = 0;
    for (; i + RSIMD_LANES_32 <= end; i += RSIMD_LANES_32) {
      rsimd_vi32 a = sx ? vx : bytes ? rsimd_vi32_loadu_u8(bx + i) : rsimd_vi32_loadu(x + i);
      rsimd_vi32 b = sy ? vy : bytes ? rsimd_vi32_loadu_u8(by + i) : rsimd_vi32_loadu(y + i);
      rsimd_mi32 eq = rsimd_vi32_cmp_eq(a, b);
      if (bytes) {
        d = rsimd_vi32_inc(d, rsimd_mi32_not(eq));
      } else {
        rsimd_mi32 na = rsimd_vi32_na2(a, b);
        d = rsimd_vi32_inc(d, rsimd_mi32_not(rsimd_mi32_or(eq, na)));
        m = rsimd_vi32_inc(m, na);
      }
    }
    r->i64 += rsimd_vi32_reduce_add(d);
    if (!bytes) missing = rsimd_vi32_reduce_add(m);
    if (missing > 0) {
      r->saw_na = 1;
      if (!na_rm) return;
    }
  }
  if (bytes) rsimd_ham_u8_from(bx, by, i, n, flags, r);
  else rsimd_ham_i32_from(x, y, i, n, flags, r, na_rm);
}

void RSIMD_KERNEL(hamming_i32)(const int *x, const int *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_i32)(const int *x, const int *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  if (flags == 0) RSIMD_KERNEL(hamming_i32_)(x, y, n, 0, r, o->na_rm, 0);
  else RSIMD_KERNEL(hamming_i32_)(x, y, n, flags, r, o->na_rm, 0);
}

void RSIMD_KERNEL(hamming_u8)(const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                              rsimd_reduce_result *r);
void RSIMD_KERNEL(hamming_u8)(const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags,
                              rsimd_reduce_result *r) {
  if (flags == 0) {
    RSIMD_KERNEL(hamming_i32_)((const int32_t *) x, (const int32_t *) y, n, 0, r, 1, 1);
  } else {
    RSIMD_KERNEL(hamming_i32_)((const int32_t *) x, (const int32_t *) y, n, flags, r, 1, 1);
  }
}

#ifndef RSIMD_NO_F64_SIMD

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(hamming_f64_)(const void *x, const void *y, R_xlen_t n,
                                                    int flags, rsimd_reduce_result *r,
                                                    int na_rm) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_ew_bcast(y, flags, 1, 1);
  ptrdiff_t i = 0;
  while (i + RSIMD_LANES_64 <= n) {
    const ptrdiff_t end = (n - i) / RSIMD_LANES_64 > RSIMD_HAMMING_BLOCK
                              ? i + RSIMD_HAMMING_BLOCK * RSIMD_LANES_64
                              : n;
    rsimd_vi64 ne = rsimd_vi64_zero(), m = rsimd_vi64_zero();
    int64_t missing;
    for (; i + RSIMD_LANES_64 <= end; i += RSIMD_LANES_64) {
      rsimd_vf64 a = rsimd_ew_ld(x, flags, 0, bc0, i, 1), b = rsimd_ew_ld(y, flags, 1, bc1, i, 1);
      ne = rsimd_vi64_inc(ne, rsimd_mf64_to_mi64(rsimd_vf64_cmp_ne(a, b)));
      m = rsimd_vi64_inc(
          m, rsimd_mf64_to_mi64(rsimd_mf64_or(rsimd_vf64_is_nan(a), rsimd_vf64_is_nan(b))));
    }
    missing = rsimd_vi64_reduce_add(m);
    r->i64 += rsimd_vi64_reduce_add(ne) - missing;
    if (missing > 0) {
      r->saw_na = 1;
      if (!na_rm) return;
    }
  }
  rsimd_ham_f64_from(x, y, i, n, flags, r, na_rm);
}

void RSIMD_KERNEL(hamming_f64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_f64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  if (flags == 0) RSIMD_KERNEL(hamming_f64_)(x, y, n, 0, r, o->na_rm);
  else RSIMD_KERNEL(hamming_f64_)(x, y, n, flags, r, o->na_rm);
}

#else /* RSIMD_NO_F64_SIMD */
#define RSIMD_SKIP_hamming_f64 1
#endif

#undef RSIMD_CMP_DISPATCH

#endif /* vector tiers */
