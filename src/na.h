#ifndef RSIMD_NA_H
#define RSIMD_NA_H

/* NA semantics, precision modes and integer overflow: the helpers every
 * kernel uses, so that each rule is written once.
 *
 * The first part is plain C and may be included anywhere (rvec.c uses it to
 * finish reductions). The second part, the vector helpers, is compiled only
 * inside tier translation units: kernels/common.inc.h includes this header
 * again after its vector layer. Each rule has a scalar form, used by the
 * none tier, and a vector form, used by every other tier; the none tier is
 * the reference the vector forms are tested against.
 *
 * Rules
 * -----
 * Missing values. R's NA_real_ is a NaN whose low 32 bits are 1954; any
 * other NaN is NaN. Integer and logical NA is INT32_MIN, integer64 NA is
 * INT64_MIN. SIMD min, max, blends and conversions do not keep NaN
 * payloads, so kernels never infer NA-ness from a result: they test the
 * inputs and record masks.
 *
 * Reductions record, per chunk, whether they saw an NA and whether they saw
 * any NaN (saw_na, saw_nan in rsimd_reduce_result) and how many elements
 * survived (count). rsimd_reduce_finish() then applies:
 *   - na.rm = FALSE: NA if any NA was seen, else NaN if any NaN was seen,
 *     whatever their order (base R's result depends on the order);
 *   - na.rm = TRUE: NA and NaN are both removed; with no survivors, sum is
 *     0, prod 1, mean NaN, min/max Inf/-Inf with base R's warning, var/sd
 *     NA (as with fewer than two survivors);
 *   - any: TRUE if any TRUE, else NA if any NA (na.rm = FALSE), else FALSE;
 *     all: FALSE if any FALSE, else NA if any NA, else TRUE. Kernels for
 *     any, all and any_na may stop early once the answer is known.
 * na.rm = TRUE always detects missing values. With na.rm = FALSE and
 * na_check = 0 (the caller promises there are none) kernels skip the masks;
 * a missing value then gives an unspecified but memory-safe result (an
 * integer NA counts as INT32_MIN; a double NA may come back as NaN).
 *
 * Elementwise doubles: arithmetic propagates NaN by itself, but which
 * payload survives when both operands are NaN differs between CPUs. Kernels
 * whose result is missing when an input is (add, sub, mul, div, fma, pmin,
 * pmax, clamp, lerp, math) put NA_real_ wherever an input is NA
 * (rsimd_na_merge_f64 / rsimd_vf64_na_merge), unless na_check = 0.
 * Elementwise integers: NA in gives NA out, for the checked and the
 * wrapping (_wrap) variants alike, unless na_check = 0.
 *
 * Integer overflow. The valid int32 range is [-INT32_MAX, INT32_MAX]:
 * INT32_MIN is NA. Checked add, sub and mul give NA for a result outside
 * it, including a wrapped result of exactly INT32_MIN (-INT32_MAX - 1L has
 * no two's complement overflow but is still NA in R), and report overflow
 * so the entry point can warn once, as base R does. neg and abs cannot
 * overflow except from INT32_MIN, which is NA, and wrapping maps INT32_MIN
 * to itself, so NA in gives NA out without a mask. The _wrap variants wrap
 * and never report overflow; a wrapped result of INT32_MIN reads as NA in
 * R. Integer %/% and %% have no SIMD instruction: they divide in double,
 * which is exact for int32 operands (x / y is never within 2^-53 relative
 * of an integer it is not equal to), take floor, and give NA for a zero
 * divisor or an NA operand (base R).
 *
 * Integer sums accumulate in 64-bit lanes; rsimd_merge_i64_sum() moves to
 * double accumulation in the (only theoretically reachable, beyond 2^32
 * elements) case that the 64-bit total would overflow.
 *
 * Precision modes (rsimd_opts.precision) for floating-point sums:
 *   - fast: four vector accumulators of W lanes each, filled in blocks of
 *     4W elements; the remaining elements go into accumulator 0 a vector
 *     at a time and then as one predicated vector. The chunk's value is
 *     (acc0 + acc1) + (acc2 + acc3), then its lanes are added pairwise,
 *     ((l0 + l1) + (l2 + l3)) + ..., and chunks are added in order. The
 *     result depends on W, so it differs between tiers, and on chunking.
 *   - pairwise: the input is cut into leaves of RSIMD_PAIRWISE_LEAF (128)
 *     elements at absolute positions; each leaf is summed as in fast mode
 *     and the leaf sums are combined by a binary counter
 *     (rsimd_pairwise_push), which is a balanced tree independent of
 *     chunking as long as chunks start at multiples of 128 (the chunk loop
 *     guarantees that). Differs between tiers only through the leaf sums.
 *   - compensated: Neumaier's improved Kahan summation in every lane of
 *     four accumulators; a chunk's lanes are combined with the same scheme
 *     in a fixed order and the chunk is added to the running (f64, comp)
 *     pair with it too. Its result is s + c (just s when s is not finite,
 *     where c would be NaN). Needs no long double, and gives
 *     sum(c(1e16, 1, -1e16)) == 1.
 * The none tier uses the same fast scheme with W = 1 and a plain sequential
 * Neumaier loop.
 */

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "kernel_types.h"

#define RSIMD_NA_I32 INT32_MIN
#define RSIMD_NA_I64 INT64_MIN
/* Leaf size of pairwise summation; reduction chunks start at multiples of
   it. */
#define RSIMD_PAIRWISE_LEAF 128
/* Bits of NA_real_: an exponent of all ones and low word 1954. */
#define RSIMD_NA_REAL_BITS UINT64_C(0x7FF00000000007A2)
/* -2147483648.0, NA_integer_ as a double, for kernels that compute integer
   results in double lanes. */
#define RSIMD_NA_I32_AS_F64 (-2147483648.0)

/* ---- Scalar helpers ------------------------------------------------------ */

/* R_IsNA(): a NaN whose low 32 bits are 1954. */
static inline int rsimd_is_na_f64(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return isnan(x) && (uint32_t) u == 1954u;
}
/* R_IsNaN(): a NaN that is not NA. */
static inline int rsimd_is_nan_f64(double x) { return isnan(x) && !rsimd_is_na_f64(x); }
static inline double rsimd_na_real(void) {
  uint64_t u = RSIMD_NA_REAL_BITS;
  double x;
  memcpy(&x, &u, sizeof x);
  return x;
}

/* r, or NA_real_ if x or y is NA (elementwise doubles). */
static inline double rsimd_na_merge_f64(double r, double x, double y) {
  return rsimd_is_na_f64(x) || rsimd_is_na_f64(y) ? rsimd_na_real() : r;
}

/* Integer arithmetic. `check` is na_check: NA operands give NA. The
   checked forms set *ovf when a result is out of range (NA operands never
   count as overflow). Wrapping goes through uint32_t, so there is no
   signed-overflow undefined behaviour. */
static inline int32_t rsimd_wrap_add_i32(int32_t x, int32_t y) {
  return (int32_t) ((uint32_t) x + (uint32_t) y);
}
static inline int32_t rsimd_wrap_sub_i32(int32_t x, int32_t y) {
  return (int32_t) ((uint32_t) x - (uint32_t) y);
}
static inline int32_t rsimd_wrap_mul_i32(int32_t x, int32_t y) {
  return (int32_t) ((uint32_t) x * (uint32_t) y);
}
static inline int rsimd_na2_i32(int32_t x, int32_t y) {
  return x == RSIMD_NA_I32 || y == RSIMD_NA_I32;
}

static inline int32_t rsimd_add_i32(int32_t x, int32_t y, int check, int *ovf) {
  int32_t r = rsimd_wrap_add_i32(x, y);
  if (check && rsimd_na2_i32(x, y)) return RSIMD_NA_I32;
  if (((x ^ r) & (y ^ r)) < 0 || r == RSIMD_NA_I32) {
    *ovf = 1;
    return RSIMD_NA_I32;
  }
  return r;
}
static inline int32_t rsimd_sub_i32(int32_t x, int32_t y, int check, int *ovf) {
  int32_t r = rsimd_wrap_sub_i32(x, y);
  if (check && rsimd_na2_i32(x, y)) return RSIMD_NA_I32;
  if (((x ^ y) & (x ^ r)) < 0 || r == RSIMD_NA_I32) {
    *ovf = 1;
    return RSIMD_NA_I32;
  }
  return r;
}
static inline int32_t rsimd_mul_i32(int32_t x, int32_t y, int check, int *ovf) {
  int64_t p = (int64_t) x * y;
  if (check && rsimd_na2_i32(x, y)) return RSIMD_NA_I32;
  if (p > INT32_MAX || p < -INT32_MAX) {
    *ovf = 1;
    return RSIMD_NA_I32;
  }
  return (int32_t) p;
}
/* Wrapping: NA (INT32_MIN) maps to itself, so these need no check. */
static inline int32_t rsimd_neg_i32(int32_t x) { return (int32_t) (0u - (uint32_t) x); }
static inline int32_t rsimd_abs_i32(int32_t x) { return x < 0 ? rsimd_neg_i32(x) : x; }

static inline int32_t rsimd_add_wrap_i32(int32_t x, int32_t y, int check) {
  return check && rsimd_na2_i32(x, y) ? RSIMD_NA_I32 : rsimd_wrap_add_i32(x, y);
}
static inline int32_t rsimd_sub_wrap_i32(int32_t x, int32_t y, int check) {
  return check && rsimd_na2_i32(x, y) ? RSIMD_NA_I32 : rsimd_wrap_sub_i32(x, y);
}
static inline int32_t rsimd_mul_wrap_i32(int32_t x, int32_t y, int check) {
  return check && rsimd_na2_i32(x, y) ? RSIMD_NA_I32 : rsimd_wrap_mul_i32(x, y);
}

/* x %/% y (mod = 0) or x %% y (mod = 1) as in base R: floor division, the
   remainder takes the sign of y; NA for y == 0 or an NA operand (or, with
   check = 0, a quotient outside the int32 range). */
static inline int32_t rsimd_intdiv_i32(int32_t x, int32_t y, int mod, int check) {
  double xd = x, yd = y, q;
  if (y == 0 || (check && rsimd_na2_i32(x, y))) return RSIMD_NA_I32;
  q = floor(xd / yd);
  if (!(q >= -INT32_MAX && q <= INT32_MAX)) return RSIMD_NA_I32;
  return (int32_t) (mod ? xd - q * yd : q);
}

/* Floating-point accumulation. */

/* Adds x to the Neumaier pair (*s, *c). */
static inline void rsimd_neumaier_add(double *s, double *c, double x) {
  double t = *s + x;
  *c += fabs(*s) >= fabs(x) ? (*s - t) + x : (x - t) + *s;
  *s = t;
}
/* The value of a Neumaier pair: s + c, or s when s is Inf or NaN (c is NaN
   then, and the sum is s). */
static inline double rsimd_neumaier_value(double s, double c) {
  return isfinite(s) ? s + c : s;
}

/* Pairwise mode: leaf sums are pushed in input order; s[k] holds the sum
   of 2^k leaves for every set bit k of `leaves`. */
static inline void rsimd_pairwise_push(rsimd_pairwise *pw, double v) {
  uint64_t k = pw->leaves;
  int level = 0;
  while (k & 1u) {
    v = pw->s[level] + v;
    k >>= 1;
    level++;
  }
  pw->s[level] = v;
  pw->leaves++;
}
/* The sum of all leaves pushed: the partial sums from the most recent
   (lowest level) up, each added to the right of the older ones. 0 when no
   leaf was pushed. */
static inline double rsimd_pairwise_total(const rsimd_pairwise *pw) {
  uint64_t k = pw->leaves;
  double acc = 0.0;
  int level, any = 0;
  for (level = 0; k != 0; level++, k >>= 1) {
    if (k & 1u) {
      acc = any ? pw->s[level] + acc : pw->s[level];
      any = 1;
    }
  }
  return acc;
}

/* The floating-point value a reduction has accumulated in `precision`
   mode: f64 (fast), the Neumaier pair (compensated), the leaf tree
   (pairwise, or f64 if no leaf was pushed). */
static inline double rsimd_reduce_value(const rsimd_reduce_result *r, int precision) {
  switch (precision) {
  case RSIMD_PREC_COMPENSATED: return rsimd_neumaier_value(r->f64, r->comp);
  case RSIMD_PREC_PAIRWISE: return r->pw.leaves ? rsimd_pairwise_total(&r->pw) : r->f64;
  default: return r->f64;
  }
}
/* Stores rsimd_reduce_value() in f64 and clears comp and the leaf tree, so
   an entry point can post-process the value (mean, var ...) in f64. */
static inline void rsimd_reduce_settle(rsimd_reduce_result *r, int precision) {
  r->f64 = rsimd_reduce_value(r, precision);
  r->comp = 0.0;
  r->pw.leaves = 0;
}

/* Adds a chunk's integer sum to r->i64; on 64-bit overflow the running
   total moves to r->f64 and r->overflow is set (the result is then
   f64 + i64, a double). */
static inline void rsimd_merge_i64_sum(rsimd_reduce_result *r, int64_t chunk) {
  if ((chunk > 0 && r->i64 > INT64_MAX - chunk) || (chunk < 0 && r->i64 < INT64_MIN - chunk)) {
    r->f64 += (double) r->i64;
    r->i64 = chunk;
    r->overflow = 1;
  } else {
    r->i64 += chunk;
  }
}

/* ---- Scalar folds (the none tier) ---------------------------------------- */

/* What a floating-point accumulation adds per element: x, x^2, |x| or x*y
   (sum, sum_sq, sum_abs, dot). */
enum { RSIMD_TERM_X = 0, RSIMD_TERM_SQ = 1, RSIMD_TERM_ABS = 2, RSIMD_TERM_XY = 3 };

/* Early exit of the logical fold: none, at the first TRUE (any), at the
   first FALSE (all). */
enum { RSIMD_STOP_NONE = 0, RSIMD_STOP_TRUE = 1, RSIMD_STOP_FALSE = 2 };

/* Elements of the input that a fold missed when it is told na_check: an
   element is missing when x (or, for x*y, y) is NaN, and saw_na is set when
   one of those is NA. */
typedef struct {
  int check, narm;
  int nan, na;
  ptrdiff_t removed;
} rsimd_fold_state;

static inline double rsimd_fold_term(rsimd_fold_state *st, int term, double x, double y) {
  double t = term == RSIMD_TERM_SQ ? x * x
             : term == RSIMD_TERM_ABS ? fabs(x)
             : term == RSIMD_TERM_XY ? x * y
             : x;
  if (st->check) {
    int missing = isnan(x) || (term == RSIMD_TERM_XY && isnan(y));
    if (missing) {
      st->nan = 1;
      if (rsimd_is_na_f64(x) || (term == RSIMD_TERM_XY && rsimd_is_na_f64(y))) st->na = 1;
      if (st->narm) {
        t = 0.0;
        st->removed++;
      }
    }
  }
  return t;
}

static inline void rsimd_fold_state_done(const rsimd_fold_state *st, ptrdiff_t n,
                                         rsimd_reduce_result *r) {
  if (st->nan) r->saw_nan = 1;
  if (st->na) r->saw_na = 1;
  r->count += n - st->removed;
}

/* Fast mode with W = 1: four accumulators, blocks of four, the rest into
   accumulator 0. Returns the sum of the n elements. */
static inline double rsimd_fold_fast_f64_(rsimd_fold_state *st, const double *x, const double *y,
                                          ptrdiff_t n, int term) {
  double a[4] = {0.0, 0.0, 0.0, 0.0};
  ptrdiff_t i = 0;
  int j;
  for (; i + 4 <= n; i += 4) {
    for (j = 0; j < 4; j++) a[j] += rsimd_fold_term(st, term, x[i + j], y ? y[i + j] : 0.0);
  }
  for (; i < n; i++) a[0] += rsimd_fold_term(st, term, x[i], y ? y[i] : 0.0);
  return (a[0] + a[1]) + (a[2] + a[3]);
}

/* Folds the n elements of one chunk into r in mode o->precision. y is used
   only for RSIMD_TERM_XY. In pairwise mode the chunk must start at a
   multiple of RSIMD_PAIRWISE_LEAF in the whole input (the chunk loops
   ensure that). */
static inline void rsimd_fold_f64(const double *x, const double *y, ptrdiff_t n, int term,
                                  rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_state st = {o->na_check || o->na_rm, o->na_rm, 0, 0, 0};
  ptrdiff_t i;
  if (term != RSIMD_TERM_XY) y = NULL;
  switch (o->precision) {
  case RSIMD_PREC_PAIRWISE:
    for (i = 0; i < n; i += RSIMD_PAIRWISE_LEAF) {
      ptrdiff_t len = n - i < RSIMD_PAIRWISE_LEAF ? n - i : RSIMD_PAIRWISE_LEAF;
      rsimd_pairwise_push(&r->pw, rsimd_fold_fast_f64_(&st, x + i, y ? y + i : NULL, len, term));
    }
    break;
  case RSIMD_PREC_COMPENSATED:
    for (i = 0; i < n; i++) {
      rsimd_neumaier_add(&r->f64, &r->comp, rsimd_fold_term(&st, term, x[i], y ? y[i] : 0.0));
    }
    break;
  default: r->f64 += rsimd_fold_fast_f64_(&st, x, y, n, term); break;
  }
  rsimd_fold_state_done(&st, n, r);
}

/* Integer and logical sum: exact in 64 bits, NA excluded under na.rm. */
static inline void rsimd_fold_sum_i32(const int32_t *x, ptrdiff_t n, rsimd_reduce_result *r,
                                      const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  int64_t s = 0;
  ptrdiff_t i, removed = 0;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I32) {
      r->saw_na = 1;
      if (o->na_rm) {
        removed++;
        continue;
      }
    }
    s += x[i];
  }
  rsimd_merge_i64_sum(r, s);
  r->count += n - removed;
}

/* Logical flags for any/all: any_true, any_false, saw_na. A value other
   than 0 and NA is TRUE. With na_check = 0 and na.rm = FALSE an NA counts
   as TRUE. Stops after the first TRUE or FALSE as `stop` says. */
static inline void rsimd_fold_lgl(const int32_t *x, ptrdiff_t n, int stop, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  ptrdiff_t i;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I32) {
      r->saw_na = 1;
    } else if (x[i] == 0) {
      r->any_false = 1;
      if (stop == RSIMD_STOP_FALSE) return;
    } else {
      r->any_true = 1;
      if (stop == RSIMD_STOP_TRUE) return;
    }
  }
}

#endif /* RSIMD_NA_H */

/* ---- Vector helpers (tier translation units only) ------------------------ */

#if defined(RSIMD_KERNELS_COMMON_INC_H) && !defined(RSIMD_NA_VECTOR_H)
#define RSIMD_NA_VECTOR_H

/* Most lanes a 64-bit vector can have (2048-bit SVE). */
#define RSIMD_MAX_LANES_64 32

#if defined(__GNUC__) || defined(__clang__)
#define RSIMD_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#define RSIMD_ALWAYS_INLINE static inline
#endif

/* All-false masks. */
RSIMD_INLINE rsimd_mi32 rsimd_mi32_none(void) {
  return rsimd_vi32_cmp_gt(rsimd_vi32_zero(), rsimd_vi32_zero());
}
RSIMD_INLINE rsimd_mi64 rsimd_mi64_none(void) {
  return rsimd_vi64_cmp_gt(rsimd_vi64_zero(), rsimd_vi64_zero());
}

/* Integer arithmetic (32-bit lanes); r is the wrapped result of the op.
   The overflow masks follow the scalar rules above. */
RSIMD_INLINE rsimd_mi32 rsimd_vi32_na2(rsimd_vi32 x, rsimd_vi32 y) {
  return rsimd_mi32_or(rsimd_vi32_is_na(x), rsimd_vi32_is_na(y));
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_add_ovf(rsimd_vi32 x, rsimd_vi32 y, rsimd_vi32 r) {
  rsimd_vi32 t = rsimd_vi32_and(rsimd_vi32_xor(x, r), rsimd_vi32_xor(y, r));
  return rsimd_mi32_or(rsimd_vi32_cmp_lt(t, rsimd_vi32_zero()), rsimd_vi32_is_na(r));
}
RSIMD_INLINE rsimd_mi32 rsimd_vi32_sub_ovf(rsimd_vi32 x, rsimd_vi32 y, rsimd_vi32 r) {
  rsimd_vi32 t = rsimd_vi32_and(rsimd_vi32_xor(x, y), rsimd_vi32_xor(x, r));
  return rsimd_mi32_or(rsimd_vi32_cmp_lt(t, rsimd_vi32_zero()), rsimd_vi32_is_na(r));
}
/* r = rsimd_vi32_mul(x, y): the product fits in int32 when its high word
   is the sign extension of r (all ones for negative r, else zero). */
RSIMD_INLINE rsimd_mi32 rsimd_vi32_mul_ovf(rsimd_vi32 x, rsimd_vi32 y, rsimd_vi32 r) {
  rsimd_vi32 sign = rsimd_vi32_blend(rsimd_vi32_zero(), rsimd_vi32_set1(-1),
                                     rsimd_vi32_cmp_lt(r, rsimd_vi32_zero()));
  rsimd_mi32 fits = rsimd_vi32_cmp_eq(rsimd_vi32_mulhi(x, y), sign);
  return rsimd_mi32_or(rsimd_mi32_not(fits), rsimd_vi32_is_na(r));
}
/* r with NA_integer_ in the lanes of m. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_set_na(rsimd_vi32 r, rsimd_mi32 m) {
  return rsimd_vi32_blend(r, rsimd_vi32_set1(RSIMD_NA_I32), m);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_neg_wrap(rsimd_vi32 x) {
  return rsimd_vi32_sub(rsimd_vi32_zero(), x);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_abs_wrap(rsimd_vi32 x) {
  return rsimd_vi32_max(x, rsimd_vi32_neg_wrap(x));
}

#ifndef RSIMD_NO_F64_SIMD
RSIMD_INLINE rsimd_mf64 rsimd_mf64_none(void) {
  return rsimd_vf64_cmp_lt(rsimd_vf64_zero(), rsimd_vf64_zero());
}

/* r with NA_real_ wherever x or y is NA (elementwise doubles). */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_na_merge(rsimd_vf64 r, rsimd_vf64 x, rsimd_vf64 y) {
  rsimd_mf64 m = rsimd_mf64_or(rsimd_vf64_is_na(x), rsimd_vf64_is_na(y));
  return rsimd_vf64_blend(r, rsimd_vf64_set1(rsimd_na_real()), m);
}

/* %/% (mod = 0) or %% (mod = 1) of int32 values held in double lanes (from
   rsimd_vf64_loadu_i32); NA lanes come back as RSIMD_NA_I32_AS_F64, so the
   result can be stored with rsimd_vf64_storeu_i32. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_intdiv(rsimd_vf64 x, rsimd_vf64 y, int mod, int check) {
  rsimd_vf64 q = rsimd_vf64_floor(rsimd_vf64_div(x, y));
  rsimd_mf64 ok = rsimd_mf64_and(rsimd_vf64_cmp_ge(q, rsimd_vf64_set1(-(double) INT32_MAX)),
                                 rsimd_vf64_cmp_le(q, rsimd_vf64_set1((double) INT32_MAX)));
  rsimd_mf64 bad = rsimd_mf64_or(rsimd_vf64_cmp_eq(y, rsimd_vf64_zero()), rsimd_mf64_not(ok));
  rsimd_vf64 r = mod ? rsimd_vf64_sub(x, rsimd_vf64_mul(q, y)) : q;
  if (check) {
    rsimd_vf64 na = rsimd_vf64_set1(RSIMD_NA_I32_AS_F64);
    bad = rsimd_mf64_or(bad, rsimd_mf64_or(rsimd_vf64_cmp_eq(x, na), rsimd_vf64_cmp_eq(y, na)));
  }
  return rsimd_vf64_blend(r, rsimd_vf64_set1(RSIMD_NA_I32_AS_F64), bad);
}

/* The lanes of v added pairwise: ((l0 + l1) + (l2 + l3)) + ...; with an
   odd count the last lane moves up a level unchanged. */
RSIMD_INLINE double rsimd_vf64_sum_tree(rsimd_vf64 v) {
  double b[RSIMD_MAX_LANES_64];
  ptrdiff_t w = RSIMD_LANES_64, j;
  rsimd_vf64_storeu(b, v);
  while (w > 1) {
    ptrdiff_t h = w / 2;
    for (j = 0; j < h; j++) b[j] = b[2 * j] + b[2 * j + 1];
    if (w & 1) b[h++] = b[w - 1];
    w = h;
  }
  return b[0];
}

/* One Neumaier step in every lane of (s, c). */
#define RSIMD_VF64_NEUMAIER(s, c, v)                                                       \
  do {                                                                                     \
    rsimd_vf64 t_ = rsimd_vf64_add((s), (v));                                              \
    rsimd_mf64 big_ = rsimd_vf64_cmp_ge(rsimd_vf64_abs(s), rsimd_vf64_abs(v));             \
    rsimd_vf64 d_ = rsimd_vf64_blend(rsimd_vf64_add(rsimd_vf64_sub((v), t_), (s)),         \
                                     rsimd_vf64_add(rsimd_vf64_sub((s), t_), (v)), big_);  \
    (c) = rsimd_vf64_add((c), d_);                                                         \
    (s) = t_;                                                                              \
  } while (0)
/* 1 if any of the n doubles is NA, stopping at the first. */
RSIMD_INLINE int rsimd_vf64_any_na(const double *x, ptrdiff_t n) {
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    if (rsimd_mf64_any(rsimd_vf64_is_na(rsimd_vf64_loadu(x + i)))) return 1;
  }
  return i < n && rsimd_mf64_any(rsimd_vf64_is_na(
                    rsimd_vf64_loadu_p(rsimd_p64_while(i, n), x + i, 0.0)));
}

/* The vector folds. They mirror the scalar folds, with the term and the
   na_check/na.rm choice as compile-time constants (each combination is a
   separate copy of the loop). Missing values are tracked as one running
   NaN mask; the chunk is rescanned for NA only if it had a NaN. */

/* t = term of vectors vx, vy; updates mnan and removed. */
#define RSIMD_VFOLD_TERM_(t, vx, vy)                                                     \
  do {                                                                                   \
    (t) = (vx);                                                                          \
    if (term == RSIMD_TERM_SQ) (t) = rsimd_vf64_mul((vx), (vx));                         \
    else if (term == RSIMD_TERM_ABS) (t) = rsimd_vf64_abs(vx);                           \
    else if (term == RSIMD_TERM_XY) (t) = rsimd_vf64_mul((vx), (vy));                    \
    if (check) {                                                                         \
      rsimd_mf64 m_ = rsimd_vf64_is_nan(vx);                                             \
      if (term == RSIMD_TERM_XY) m_ = rsimd_mf64_or(m_, rsimd_vf64_is_nan(vy));          \
      mnan = rsimd_mf64_or(mnan, m_);                                                    \
      if (narm) {                                                                        \
        (t) = rsimd_vf64_blend((t), rsimd_vf64_zero(), m_);                              \
        removed += rsimd_mf64_count(m_);                                                 \
      }                                                                                  \
    }                                                                                    \
  } while (0)

/* Loads of x and y at i: full vectors, or the predicated tail (fill 0). */
#define RSIMD_VFOLD_LOAD_(vx, vy, i)                                                     \
  rsimd_vf64 vx = rsimd_vf64_loadu(x + (i));                                             \
  rsimd_vf64 vy = term == RSIMD_TERM_XY ? rsimd_vf64_loadu(y + (i)) : vx
#define RSIMD_VFOLD_LOAD_P_(vx, vy, pg, i)                                               \
  rsimd_vf64 vx = rsimd_vf64_loadu_p((pg), x + (i), 0.0);                                \
  rsimd_vf64 vy = term == RSIMD_TERM_XY ? rsimd_vf64_loadu_p((pg), y + (i), 0.0) : vx

RSIMD_ALWAYS_INLINE void rsimd_vfold_done_(rsimd_mf64 mnan, ptrdiff_t removed, const double *x,
                                           const double *y, ptrdiff_t n, int term,
                                           rsimd_reduce_result *r) {
  if (rsimd_mf64_any(mnan)) {
    r->saw_nan = 1;
    if (!r->saw_na &&
        (rsimd_vf64_any_na(x, n) || (term == RSIMD_TERM_XY && rsimd_vf64_any_na(y, n)))) {
      r->saw_na = 1;
    }
  }
  r->count += n - removed;
}

/* Fast mode: returns the sum of the n elements (see the rules above). */
RSIMD_ALWAYS_INLINE double rsimd_vfold_fast_(const double *x, const double *y, ptrdiff_t n,
                                             const int term, const int check, const int narm,
                                             rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  rsimd_vf64 a0 = rsimd_vf64_zero(), a1 = a0, a2 = a0, a3 = a0, t;
  rsimd_mf64 mnan = rsimd_mf64_none();
  ptrdiff_t i = 0, removed = 0;
  for (; i + 4 * W <= n; i += 4 * W) {
    RSIMD_VFOLD_LOAD_(x0, y0, i);
    RSIMD_VFOLD_LOAD_(x1, y1, i + W);
    RSIMD_VFOLD_LOAD_(x2, y2, i + 2 * W);
    RSIMD_VFOLD_LOAD_(x3, y3, i + 3 * W);
    RSIMD_VFOLD_TERM_(t, x0, y0);
    a0 = rsimd_vf64_add(a0, t);
    RSIMD_VFOLD_TERM_(t, x1, y1);
    a1 = rsimd_vf64_add(a1, t);
    RSIMD_VFOLD_TERM_(t, x2, y2);
    a2 = rsimd_vf64_add(a2, t);
    RSIMD_VFOLD_TERM_(t, x3, y3);
    a3 = rsimd_vf64_add(a3, t);
  }
  for (; i + W <= n; i += W) {
    RSIMD_VFOLD_LOAD_(vx, vy, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    a0 = rsimd_vf64_add(a0, t);
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    RSIMD_VFOLD_LOAD_P_(vx, vy, pg, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    a0 = rsimd_vf64_add(a0, t);
  }
  if (check) rsimd_vfold_done_(mnan, removed, x, y, n, term, r);
  else r->count += n;
  return rsimd_vf64_sum_tree(rsimd_vf64_add(rsimd_vf64_add(a0, a1), rsimd_vf64_add(a2, a3)));
}

/* Compensated mode: Neumaier per lane of four accumulators, then the lanes
   (accumulator 0 lanes 0..W-1, then accumulator 1 ...) with the scalar
   scheme, then into r's running pair. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_comp_(const double *x, const double *y, ptrdiff_t n,
                                           const int term, const int check, const int narm,
                                           rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  rsimd_vf64 s0 = rsimd_vf64_zero(), s1 = s0, s2 = s0, s3 = s0;
  rsimd_vf64 c0 = s0, c1 = s0, c2 = s0, c3 = s0, t;
  rsimd_mf64 mnan = rsimd_mf64_none();
  ptrdiff_t i = 0, removed = 0, j;
  double sb[4][RSIMD_MAX_LANES_64], cb[4][RSIMD_MAX_LANES_64], S = 0.0, C = 0.0;
  int k;
  for (; i + 4 * W <= n; i += 4 * W) {
    RSIMD_VFOLD_LOAD_(x0, y0, i);
    RSIMD_VFOLD_LOAD_(x1, y1, i + W);
    RSIMD_VFOLD_LOAD_(x2, y2, i + 2 * W);
    RSIMD_VFOLD_LOAD_(x3, y3, i + 3 * W);
    RSIMD_VFOLD_TERM_(t, x0, y0);
    RSIMD_VF64_NEUMAIER(s0, c0, t);
    RSIMD_VFOLD_TERM_(t, x1, y1);
    RSIMD_VF64_NEUMAIER(s1, c1, t);
    RSIMD_VFOLD_TERM_(t, x2, y2);
    RSIMD_VF64_NEUMAIER(s2, c2, t);
    RSIMD_VFOLD_TERM_(t, x3, y3);
    RSIMD_VF64_NEUMAIER(s3, c3, t);
  }
  for (; i + W <= n; i += W) {
    RSIMD_VFOLD_LOAD_(vx, vy, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    RSIMD_VF64_NEUMAIER(s0, c0, t);
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    RSIMD_VFOLD_LOAD_P_(vx, vy, pg, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    RSIMD_VF64_NEUMAIER(s0, c0, t);
  }
  if (check) rsimd_vfold_done_(mnan, removed, x, y, n, term, r);
  else r->count += n;
  rsimd_vf64_storeu(sb[0], s0);
  rsimd_vf64_storeu(sb[1], s1);
  rsimd_vf64_storeu(sb[2], s2);
  rsimd_vf64_storeu(sb[3], s3);
  rsimd_vf64_storeu(cb[0], c0);
  rsimd_vf64_storeu(cb[1], c1);
  rsimd_vf64_storeu(cb[2], c2);
  rsimd_vf64_storeu(cb[3], c3);
  for (k = 0; k < 4; k++) {
    for (j = 0; j < W; j++) {
      rsimd_neumaier_add(&S, &C, sb[k][j]);
      C += cb[k][j];
    }
  }
  rsimd_neumaier_add(&r->f64, &r->comp, S);
  r->comp += C;
}

RSIMD_ALWAYS_INLINE void rsimd_vfold_mode_(const double *x, const double *y, ptrdiff_t n,
                                           const int term, const int check, const int narm,
                                           rsimd_reduce_result *r, int precision) {
  ptrdiff_t i;
  switch (precision) {
  case RSIMD_PREC_PAIRWISE:
    for (i = 0; i < n; i += RSIMD_PAIRWISE_LEAF) {
      ptrdiff_t len = n - i < RSIMD_PAIRWISE_LEAF ? n - i : RSIMD_PAIRWISE_LEAF;
      rsimd_pairwise_push(&r->pw, rsimd_vfold_fast_(x + i, y + i, len, term, check, narm, r));
    }
    break;
  case RSIMD_PREC_COMPENSATED: rsimd_vfold_comp_(x, y, n, term, check, narm, r); break;
  default: r->f64 += rsimd_vfold_fast_(x, y, n, term, check, narm, r); break;
  }
}

/* The vector form of rsimd_fold_f64(); `term` should be a constant. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_f64(const double *x, const double *y, ptrdiff_t n,
                                         const int term, rsimd_reduce_result *r,
                                         const rsimd_opts *o) {
  if (term != RSIMD_TERM_XY) y = x; /* never read */
  if (o->na_rm) rsimd_vfold_mode_(x, y, n, term, 1, 1, r, o->precision);
  else if (o->na_check) rsimd_vfold_mode_(x, y, n, term, 1, 0, r, o->precision);
  else rsimd_vfold_mode_(x, y, n, term, 0, 0, r, o->precision);
}
#endif /* RSIMD_NO_F64_SIMD */

/* The vector form of rsimd_fold_sum_i32(): 64-bit lanes. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_sum_i32_(const int32_t *x, ptrdiff_t n, const int check,
                                              const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vi64 na = rsimd_vi64_set1(RSIMD_NA_I32), zero = rsimd_vi64_zero();
  rsimd_vi64 a0 = zero, a1 = zero;
  rsimd_mi64 mna = rsimd_mi64_none();
  ptrdiff_t i = 0, removed = 0;
#define RSIMD_VFOLD_I32_(acc, v)                                                         \
  do {                                                                                   \
    if (check) {                                                                         \
      rsimd_mi64 m_ = rsimd_vi64_cmp_eq((v), na);                                        \
      mna = rsimd_mi64_or(mna, m_);                                                      \
      if (narm) {                                                                        \
        (v) = rsimd_vi64_blend((v), zero, m_);                                           \
        removed += rsimd_mi64_count(m_);                                                 \
      }                                                                                  \
    }                                                                                    \
    (acc) = rsimd_vi64_add((acc), (v));                                                  \
  } while (0)
  for (; i + 2 * W <= n; i += 2 * W) {
    rsimd_vi64 v0 = rsimd_vi64_loadu_i32(x + i), v1 = rsimd_vi64_loadu_i32(x + i + W);
    RSIMD_VFOLD_I32_(a0, v0);
    RSIMD_VFOLD_I32_(a1, v1);
  }
  for (; i < n; i += W) {
    rsimd_vi64 v = rsimd_vi64_loadu_i32_p(rsimd_p64_while(i, n), x + i, 0);
    RSIMD_VFOLD_I32_(a0, v);
  }
#undef RSIMD_VFOLD_I32_
  if (rsimd_mi64_any(mna)) r->saw_na = 1;
  rsimd_merge_i64_sum(r, rsimd_vi64_reduce_add(rsimd_vi64_add(a0, a1)));
  r->count += n - removed;
}
RSIMD_ALWAYS_INLINE void rsimd_vfold_sum_i32(const int32_t *x, ptrdiff_t n, rsimd_reduce_result *r,
                                             const rsimd_opts *o) {
  if (o->na_rm) rsimd_vfold_sum_i32_(x, n, 1, 1, r);
  else if (o->na_check) rsimd_vfold_sum_i32_(x, n, 1, 0, r);
  else rsimd_vfold_sum_i32_(x, n, 0, 0, r);
}

/* The vector form of rsimd_fold_lgl(). The tail is loaded twice, with fill
   0 for the TRUE and NA tests and fill 1 for the FALSE test, so filled
   lanes count as nothing. Early exit is checked every 4 vectors. */
RSIMD_INLINE void rsimd_vfold_lgl(const int32_t *x, ptrdiff_t n, int stop, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const int check = o->na_check || o->na_rm;
  const rsimd_vi32 zero = rsimd_vi32_zero();
  rsimd_mi32 mt = rsimd_mi32_none(), mf = mt, mna = mt;
  ptrdiff_t i = 0;
#define RSIMD_VFOLD_LGL_(vt, vf)                                                         \
  do {                                                                                   \
    rsimd_mi32 na_ = check ? rsimd_vi32_is_na(vt) : rsimd_mi32_none();                   \
    mna = rsimd_mi32_or(mna, na_);                                                       \
    mt = rsimd_mi32_or(mt, rsimd_mi32_not(rsimd_mi32_or(rsimd_vi32_cmp_eq((vt), zero), na_))); \
    mf = rsimd_mi32_or(mf, rsimd_vi32_cmp_eq((vf), zero));                               \
  } while (0)
  while (i < n) {
    ptrdiff_t end = n - i > 4 * W ? i + 4 * W : n;
    for (; i + W <= end; i += W) {
      rsimd_vi32 v = rsimd_vi32_loadu(x + i);
      RSIMD_VFOLD_LGL_(v, v);
    }
    if (i < end) {
      rsimd_p32 pg = rsimd_p32_while(i, end);
      rsimd_vi32 vt = rsimd_vi32_loadu_p(pg, x + i, 0), vf = rsimd_vi32_loadu_p(pg, x + i, 1);
      RSIMD_VFOLD_LGL_(vt, vf);
      i = end;
    }
    if ((stop == RSIMD_STOP_TRUE && rsimd_mi32_any(mt)) ||
        (stop == RSIMD_STOP_FALSE && rsimd_mi32_any(mf))) {
      break;
    }
  }
#undef RSIMD_VFOLD_LGL_
  if (rsimd_mi32_any(mt)) r->any_true = 1;
  if (rsimd_mi32_any(mf)) r->any_false = 1;
  if (rsimd_mi32_any(mna)) r->saw_na = 1;
}

#endif /* vector helpers */
