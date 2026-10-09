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
 *     whatever their order (base R's result depends on the order); var and
 *     sd are NA for either, as in base R;
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
 * whose result is missing when an input is (add, sub, mul, div, %/%, %%,
 * copysign, fma, mul_add, add_mul, lerp, math) put NA_real_ wherever an
 * input is NA (rsimd_na_merge_f64 / rsimd_vf64_na_merge), unless
 * na_check = 0. pmin, pmax and clamp instead return one of their operands,
 * payload and all, as base R does: a missing operand propagates (NA or
 * NaN as it is), the second one when both are missing; the _num forms
 * (na.rm = TRUE) skip a missing operand. Unary ops (neg, abs, sqrt,
 * floor ...) keep the payload of a NaN operand by themselves.
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
 *     at a time and then as one predicated vector. Product terms (x^2,
 *     x*y, (x - y)^2, (x - c)^2) are accumulated with a fused multiply-add
 *     on the tiers with a native one (RSIMD_NATIVE_FMA: avx2, avx512,
 *     arm64 neon, sve, sve2), so those differ from sse2 and none in the
 *     last bits; pairwise leaves and compensated mode add the rounded
 *     term. The chunk's value is
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
#include "soft_fma.h"

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
   remainder takes the sign of y; NA for y == 0, an NA operand, or a
   quotient outside [-INT32_MAX, INT32_MAX], which only an INT32_MIN operand
   can produce (INT32_MIN is NA, so only with check = 0). */
static inline int32_t rsimd_intdiv_i32(int32_t x, int32_t y, int mod, int check) {
  double xd = x, yd = y, q;
  if (y == 0 || (check && rsimd_na2_i32(x, y))) return RSIMD_NA_I32;
  q = floor(xd / yd);
  if (!(q >= -INT32_MAX && q <= INT32_MAX)) return RSIMD_NA_I32;
  return (int32_t) (mod ? xd - q * yd : q);
}

/* Division of int32 values by a constant d other than 0 and INT32_MIN
   without dividing (Hacker's Delight, 2nd ed., section 10-1): the
   truncated quotient of x is ((mulhi(m, x) + add * x) >> s) plus one when
   that is negative, with m == 0 standing for |d| == 1. */
typedef struct {
  int32_t d, m;
  int s, add;
} rsimd_divmagic_i32;
static inline rsimd_divmagic_i32 rsimd_divmagic_i32_make(int32_t d) {
  const uint32_t two31 = UINT32_C(0x80000000);
  rsimd_divmagic_i32 g;
  uint32_t ad = d < 0 ? 0u - (uint32_t) d : (uint32_t) d, t, anc, q1, r1, q2, r2, delta, m;
  int p = 31;
  g.d = d;
  g.m = 0;
  g.s = 0;
  g.add = 0;
  if (ad == 1) return g;
  t = two31 + ((uint32_t) d >> 31);
  anc = t - 1 - t % ad;
  q1 = two31 / anc;
  r1 = two31 - q1 * anc;
  q2 = two31 / ad;
  r2 = two31 - q2 * ad;
  do {
    p++;
    q1 *= 2;
    r1 *= 2;
    if (r1 >= anc) {
      q1++;
      r1 -= anc;
    }
    q2 *= 2;
    r2 *= 2;
    if (r2 >= ad) {
      q2++;
      r2 -= ad;
    }
    delta = ad - r2;
  } while (q1 < delta || (q1 == delta && r1 == 0));
  m = q2 + 1;
  if (d < 0) m = 0u - m;
  memcpy(&g.m, &m, sizeof m);
  g.s = p - 32;
  g.add = d > 0 && g.m < 0 ? 1 : (d < 0 && g.m > 0 ? -1 : 0);
  return g;
}

/* r, or NA_real_ if x, y or z is NA (elementwise doubles). */
static inline double rsimd_na_merge3_f64(double r, double x, double y, double z) {
  return rsimd_is_na_f64(z) ? rsimd_na_real() : rsimd_na_merge_f64(r, x, y);
}

/* x %/% y for a finite quotient q = x / y with |q| >= 2^52, where q is
   an integer and the exact floor F of the real quotient X can lie
   between two doubles: returns F rounded once. The residual
   x - q * y is exact (q is the rounded quotient), so its sign says
   whether X < q. If it is not, F rounds to q. Otherwise X lies in
   [m, q), with m the midpoint between q and the next double qd below
   it. F is qd when that is q - 1. If not, F rounds to q unless F = m,
   that is m <= X < m + 1, a tie that goes to the even one of q and qd.
   With h = (q - qd) / 2, (X - m) y = res + h y, which is exact
   (Sterbenz) whenever it can be smaller than y. A tiny y is scaled up
   first so that neither residual underflows. */
static inline double rsimd_idiv_big_f64(double x, double y) {
  double q, qd, res, h, d;
  uint64_t b;
  if (fabs(y) < 0x1p-900) {
    x *= 0x1p200;
    y *= 0x1p200;
  }
  q = x / y;
  res = rsimd_fma(-q, y, x);
  if (res == 0 || (res < 0) == (y < 0)) return q;
  qd = nextafter(q, -HUGE_VAL);
  if (q - qd == 1.0) return qd;
  h = (q - qd) * 0.5;
  d = rsimd_fma(h, y, res);
  if (y < 0) d = -d;
  if (d < 0 || d >= fabs(y)) return q;
  memcpy(&b, &q, sizeof b);
  return (b & 1) ? qd : q;
}

/* Base R's x %/% y for doubles (myfloor in R's arithmetic.c): the exact
   floor of the real quotient, rounded once. Below 2^52 that is
   floor(x / y), or one less when x / y rounded up to an integer; the
   residual x - k * y is computed with a single rounding (fma), so its
   sign is exact. Larger finite quotients go to rsimd_idiv_big_f64().
   -2 %/% Inf is -1 and 2 %/% Inf is 0; x %/% 0 is x / 0. A zero
   quotient is +0. */
static inline double rsimd_idiv_f64(double x, double y) {
  double q = x / y, k, res;
  if (fabs(q) >= 0x1p52 && fabs(q) < HUGE_VAL) return rsimd_idiv_big_f64(x, y);
  k = floor(q);
  res = isinf(y) ? x : rsimd_fma(-k, y, x);
  if ((res < 0 && y > 0) || (res > 0 && y < 0)) k -= 1;
  return k + 0.0;
}
/* Base R's x %% y for doubles (myfmod), exactly: x - floor(x / y) * y
   with the real quotient, rounded once, which has the sign of y. fmod is
   exact, and adding y to a remainder of the wrong sign rounds once. Zero
   is +0, except that an infinite y returns x (or y when x has the other
   sign), as base R does; x %% 0 and Inf %% y are NaN. */
static inline double rsimd_mod_f64(double x, double y) {
  double r;
  if (isnan(x) || isnan(y)) return x + y;
  if (isinf(y) && isfinite(x)) return (x < 0 && y > 0) || (x > 0 && y < 0) ? y : x;
  r = fmod(x, y);
  if ((r < 0 && y > 0) || (r > 0 && y < 0)) r += y;
  return r + 0.0;
}

/* Base R's pmin(x, y) and pmax(x, y) for one pair: a missing operand is
   returned as it is (y when both are); ties return x, so pmin(0, -0) is 0.
   The _num forms (na.rm = TRUE) return the other operand when one is
   missing, and y when both are. */
static inline double rsimd_pmin_f64(double x, double y) {
  if (isnan(y)) return y;
  if (isnan(x)) return x;
  return y < x ? y : x;
}
static inline double rsimd_pmax_f64(double x, double y) {
  if (isnan(y)) return y;
  if (isnan(x)) return x;
  return y > x ? y : x;
}
static inline double rsimd_pmin_num_f64(double x, double y) {
  if (isnan(x)) return y;
  if (isnan(y)) return x;
  return y < x ? y : x;
}
static inline double rsimd_pmax_num_f64(double x, double y) {
  if (isnan(x)) return y;
  if (isnan(y)) return x;
  return y > x ? y : x;
}
static inline int32_t rsimd_pmin_i32(int32_t x, int32_t y) {
  return rsimd_na2_i32(x, y) ? RSIMD_NA_I32 : (y < x ? y : x);
}
static inline int32_t rsimd_pmax_i32(int32_t x, int32_t y) {
  return rsimd_na2_i32(x, y) ? RSIMD_NA_I32 : (y > x ? y : x);
}
static inline int32_t rsimd_pmin_num_i32(int32_t x, int32_t y) {
  if (x == RSIMD_NA_I32) return y;
  if (y == RSIMD_NA_I32) return x;
  return y < x ? y : x;
}
static inline int32_t rsimd_pmax_num_i32(int32_t x, int32_t y) {
  if (x == RSIMD_NA_I32) return y;
  if (y == RSIMD_NA_I32) return x;
  return y > x ? y : x;
}

/* |x| with the sign of s; a NaN s (NA included) is returned as it is. A
   NaN x keeps its payload. */
static inline double rsimd_copysign_f64(double x, double s) {
  return isnan(s) ? s : copysign(x, s);
}
/* Base R's sign(): -1, 0 or 1, +0 for -0, NaN payloads kept. */
static inline double rsimd_sign_f64(double x) {
  if (x > 0) return 1.0;
  if (x < 0) return -1.0;
  return x == 0 ? 0.0 : x;
}
/* Base R's round(x, digits) (fround() in R's nmath) for an integer
   d = floor(digits + 0.5) in [-308, 308], given p10 = R_pow_di(10, d) and
   big, the smallest power of 2 at which fround() returns x as it is
   (|x| * 10^d has more than 15 significant digits). x is rounded down
   and up to d places, xd = floor(|x| * p10) / p10 and
   xu = ceil(|x| * p10) / p10, and the closer one wins, xu on a tie when
   floor(|x| * p10) is odd. Zero, infinities and NaN are returned as they
   are. This is a copy of fround() as R has had it since 4.0.0, not a call
   to it, so it has to be updated if R changes that algorithm; the tests
   compare it with round() bit for bit. */
static inline double rsimd_round_digits_f64(double x, double p10, double big) {
  double a = fabs(x), x10, i10, xd, xu, du, dd, r;
  if (!(a > 0 && a < big)) return x;
  x10 = a * p10;
  i10 = floor(x10);
  xd = i10 / p10;
  xu = ceil(x10) / p10;
  du = xu - a;
  dd = a - xd;
  r = du < dd || (fmod(i10, 2.0) == 1 && du == dd) ? xu : xd;
  return x < 0 ? -r : r;
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

/* What a floating-point accumulation adds per element: x, x^2, |x|, x*y,
   x - c for a constant c, (x - y)^2 or (x - c)^2 (sum, sum_sq, sum_abs, dot,
   the refinement pass of mean, dist, the second pass of var). */
enum {
  RSIMD_TERM_X = 0,
  RSIMD_TERM_SQ = 1,
  RSIMD_TERM_ABS = 2,
  RSIMD_TERM_XY = 3,
  RSIMD_TERM_DEV = 4,
  RSIMD_TERM_SQDIFF = 5,
  RSIMD_TERM_DEVSQ = 6
};
/* Terms that read the second operand y. */
#define RSIMD_TERM_PAIR(term) ((term) == RSIMD_TERM_XY || (term) == RSIMD_TERM_SQDIFF)
/* Terms that are a product, which fast mode accumulates with a fused
   multiply-add on tiers that have one. */
#define RSIMD_TERM_PRODUCT(term)                                                         \
  ((term) == RSIMD_TERM_SQ || (term) == RSIMD_TERM_XY || (term) == RSIMD_TERM_SQDIFF ||  \
   (term) == RSIMD_TERM_DEVSQ)

/* Element i of an operand stored as double (i32 = 0) or as int32 (i32 =
   1, integer or logical); an int32 NA reads as NA_real_. */
static inline double rsimd_elt_f64(const void *p, int i32, ptrdiff_t i) {
  if (i32) {
    int32_t v = ((const int32_t *) p)[i];
    return v == RSIMD_NA_I32 ? rsimd_na_real() : (double) v;
  }
  return ((const double *) p)[i];
}

/* Early exit of the logical fold: none, at the first TRUE (any), at the
   first FALSE (all). */
enum { RSIMD_STOP_NONE = 0, RSIMD_STOP_TRUE = 1, RSIMD_STOP_FALSE = 2 };

/* Elements of the input that a fold missed when it is told na_check: an
   element is missing when x is NaN, or for a pair term when the term is
   NaN (so also Inf * 0 and Inf - Inf, as in sum(x * y, na.rm = TRUE)), and
   saw_na is set when x (or y) is NA. */
typedef struct {
  int check, narm;
  int nan, na;
  ptrdiff_t removed;
} rsimd_fold_state;

static inline double rsimd_fold_term(rsimd_fold_state *st, int term, double x, double y,
                                     double c) {
  double d = term == RSIMD_TERM_SQDIFF ? x - y : x - c;
  double t = term == RSIMD_TERM_SQ ? x * x
             : term == RSIMD_TERM_ABS ? fabs(x)
             : term == RSIMD_TERM_XY ? x * y
             : term == RSIMD_TERM_DEV ? d
             : term == RSIMD_TERM_SQDIFF || term == RSIMD_TERM_DEVSQ ? d * d
             : x;
  if (st->check) {
    int missing = RSIMD_TERM_PAIR(term) ? isnan(t) : isnan(x);
    if (missing) {
      st->nan = 1;
      if (rsimd_is_na_f64(x) || (RSIMD_TERM_PAIR(term) && rsimd_is_na_f64(y))) st->na = 1;
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
   accumulator 0. Returns the sum of the n elements, which start at element
   i0 of x and y (each double or int32 as xi32 and yi32 say; y may be NULL
   for a term that does not read it). */
static inline double rsimd_fold_fast_f64_(rsimd_fold_state *st, const void *x, int xi32,
                                          const void *y, int yi32, ptrdiff_t i0, ptrdiff_t n,
                                          int term, double c) {
  double a[4] = {0.0, 0.0, 0.0, 0.0};
  ptrdiff_t i = 0;
  int j;
#define RSIMD_FOLD_TERM_AT_(k)                                                           \
  rsimd_fold_term(st, term, rsimd_elt_f64(x, xi32, i0 + (k)),                           \
                  y ? rsimd_elt_f64(y, yi32, i0 + (k)) : 0.0, c)
  while (i + 4 <= n) {
    /* Blocks of 4096 (a multiple of 4), after each of which the fold
       stops at an NA without na.rm: the sum is then NA whatever follows. */
    const ptrdiff_t end = n - i > 4096 ? i + 4096 : n;
    for (; i + 4 <= end; i += 4) {
      for (j = 0; j < 4; j++) a[j] += RSIMD_FOLD_TERM_AT_(i + j);
    }
    if (st->na && !st->narm) return 0.0;
  }
  for (; i < n; i++) a[0] += RSIMD_FOLD_TERM_AT_(i);
#undef RSIMD_FOLD_TERM_AT_
  return (a[0] + a[1]) + (a[2] + a[3]);
}

/* Folds the n elements of one chunk into r in mode o->precision. Without
   na.rm it stops at an NA (the count is still n), as the result is then
   NA. x and y are double (xi32, yi32 = 0) or int32 (= 1) elements; y is
   used only by the pair terms (RSIMD_TERM_PAIR) and c only by
   RSIMD_TERM_DEV and RSIMD_TERM_DEVSQ. In pairwise mode the chunk must start at a multiple of
   RSIMD_PAIRWISE_LEAF in the whole input (the chunk loops ensure that). */
static inline void rsimd_fold_gen(const void *x, int xi32, const void *y, int yi32, ptrdiff_t n,
                                  int term, double c, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  rsimd_fold_state st = {o->na_check || o->na_rm, o->na_rm, 0, 0, 0};
  ptrdiff_t i;
  if (!RSIMD_TERM_PAIR(term)) y = NULL;
  switch (o->precision) {
  case RSIMD_PREC_PAIRWISE:
    for (i = 0; i < n; i += RSIMD_PAIRWISE_LEAF) {
      ptrdiff_t len = n - i < RSIMD_PAIRWISE_LEAF ? n - i : RSIMD_PAIRWISE_LEAF;
      rsimd_pairwise_push(&r->pw, rsimd_fold_fast_f64_(&st, x, xi32, y, yi32, i, len, term, c));
      if (st.na && !st.narm) break;
    }
    break;
  case RSIMD_PREC_COMPENSATED:
    for (i = 0; i < n; i++) {
      double t = rsimd_fold_term(&st, term, rsimd_elt_f64(x, xi32, i),
                                 y ? rsimd_elt_f64(y, yi32, i) : 0.0, c);
      rsimd_neumaier_add(&r->f64, &r->comp, t);
      if (st.na && !st.narm) break;
    }
    break;
  default: r->f64 += rsimd_fold_fast_f64_(&st, x, xi32, y, yi32, 0, n, term, c); break;
  }
  rsimd_fold_state_done(&st, n, r);
}
static inline void rsimd_fold_f64_c(const double *x, const double *y, ptrdiff_t n, int term,
                                    double c, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_gen(x, 0, y, 0, n, term, c, r, o);
}
static inline void rsimd_fold_f64(const double *x, const double *y, ptrdiff_t n, int term,
                                  rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_gen(x, 0, y, 0, n, term, 0.0, r, o);
}

/* Integer and logical sum (or, with absval, sum of absolute values):
   exact in 64 bits, NA excluded under na.rm. Without na.rm the fold stops
   at the first NA (r->saw_na set, the count still n). */
static inline void rsimd_fold_sum_i32(const int32_t *x, ptrdiff_t n, int absval,
                                      rsimd_reduce_result *r, const rsimd_opts *o) {
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
      /* The sum is NA whatever follows. */
      r->count += n;
      return;
    }
    s += absval && x[i] < 0 ? -(int64_t) x[i] : x[i];
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

/* Elements per block of the folds and min/max kernels that stop early at
   an NA: a multiple of 4W on every tier, so that the blocks do not change
   the order of the additions. */
#define RSIMD_FOLD_BLOCK 4096

#if defined(__GNUC__) || defined(__clang__)
#define RSIMD_ALWAYS_INLINE static inline __attribute__((always_inline))
#define RSIMD_NOINLINE __attribute__((noinline))
#else
#define RSIMD_ALWAYS_INLINE static inline
#define RSIMD_NOINLINE
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
/* x %/% d (mod = 0) or x %% d (mod = 1) for the constant d of g, as
   rsimd_intdiv_i32(): the truncated quotient from rsimd_divmagic_i32 is
   stepped down by one where the remainder is non-zero and of the other
   sign than d. NA lanes of x give NA when check is set, and so does
   INT32_MIN for |d| == 1, whose quotient is out of range. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_intdiv_const(rsimd_vi32 x, const rsimd_divmagic_i32 *g, int mod,
                                                int check) {
  const rsimd_vi32 zero = rsimd_vi32_zero(), d = rsimd_vi32_set1(g->d);
  rsimd_vi32 q, r;
  rsimd_mi32 low;
  if (g->m == 0) {
    q = g->d > 0 ? x : rsimd_vi32_neg_wrap(x);
  } else {
    q = rsimd_vi32_mulhi(x, rsimd_vi32_set1(g->m));
    if (g->add > 0) q = rsimd_vi32_add(q, x);
    else if (g->add < 0) q = rsimd_vi32_sub(q, x);
    q = rsimd_vi32_sra(q, g->s);
    q = rsimd_vi32_add(q, rsimd_vi32_srl(q, 31));
  }
  r = rsimd_vi32_sub(x, rsimd_vi32_mul(q, d));
  low = rsimd_mi32_andnot(rsimd_vi32_cmp_eq(r, zero),
                          rsimd_vi32_cmp_lt(rsimd_vi32_xor(r, d), zero));
  r = mod ? rsimd_vi32_blend(r, rsimd_vi32_add(r, d), low)
          : rsimd_vi32_blend(q, rsimd_vi32_sub(q, rsimd_vi32_set1(1)), low);
  if (check || g->m == 0) r = rsimd_vi32_set_na(r, rsimd_vi32_is_na(x));
  return r;
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
  const rsimd_vf64 na = rsimd_vf64_set1(RSIMD_NA_I32_AS_F64);
  rsimd_vf64 q = rsimd_vf64_floor(rsimd_vf64_div(x, y));
  rsimd_vf64 r = mod ? rsimd_vf64_sub(x, rsimd_vf64_mul(q, y)) : q;
  rsimd_mf64 bad = rsimd_vf64_cmp_eq(y, rsimd_vf64_zero());
  /* Only an INT32_MIN operand can put a quotient outside [-INT32_MAX,
     INT32_MAX] (INT32_MIN %/% -1, which is 2^31, or INT32_MIN %/% 1, which
     is -2^31, the NA pattern), so with the NA check that range test is not
     needed. */
  if (check) {
    bad = rsimd_mf64_or(bad, rsimd_mf64_or(rsimd_vf64_cmp_eq(x, na), rsimd_vf64_cmp_eq(y, na)));
  } else {
    bad = rsimd_mf64_or(bad, rsimd_vf64_cmp_gt(rsimd_vf64_abs(q), rsimd_vf64_set1((double) INT32_MAX)));
  }
  return rsimd_vf64_blend(r, na, bad);
}

/* r with NA_real_ wherever x, y or z is NA. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_na_merge3(rsimd_vf64 r, rsimd_vf64 x, rsimd_vf64 y,
                                             rsimd_vf64 z) {
  rsimd_mf64 m = rsimd_mf64_or(rsimd_mf64_or(rsimd_vf64_is_na(x), rsimd_vf64_is_na(y)),
                               rsimd_vf64_is_na(z));
  return rsimd_vf64_blend(r, rsimd_vf64_set1(rsimd_na_real()), m);
}

/* x %/% y for doubles as rsimd_idiv_f64(), for the lanes where
   |x / y| < 2^52. When slow is not NULL, *slow is set to the lanes with
   a larger quotient (infinities included, which are rare and which the
   scalar form handles too), which the caller recomputes with
   rsimd_idiv_f64(). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_vf64_idiv(rsimd_vf64 x, rsimd_vf64 y, rsimd_mf64 *slow) {
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_vf64 q = rsimd_vf64_div(x, y), k = rsimd_vf64_floor(q);
  rsimd_mf64 yinf = rsimd_vf64_cmp_eq(rsimd_vf64_abs(y), rsimd_vf64_set1(HUGE_VAL));
  rsimd_vf64 res = rsimd_vf64_blend(rsimd_vf64_fma(rsimd_vf64_neg(k), y, x), x, yinf);
  rsimd_mf64 low =
    rsimd_mf64_or(rsimd_mf64_and(rsimd_vf64_cmp_lt(res, zero), rsimd_vf64_cmp_gt(y, zero)),
                  rsimd_mf64_and(rsimd_vf64_cmp_gt(res, zero), rsimd_vf64_cmp_lt(y, zero)));
  k = rsimd_vf64_blend(k, rsimd_vf64_sub(k, rsimd_vf64_set1(1.0)), low);
  if (slow) *slow = rsimd_vf64_cmp_ge(rsimd_vf64_abs(q), rsimd_vf64_set1(0x1p52));
  return rsimd_vf64_add(k, zero);
}

/* x %% y for doubles as rsimd_mod_f64(), for the lanes where |x / y| <
   2^52, so that the corrected quotient k and k - 1 are exact and
   x - k * y rounds once. *slow is set to the finite lanes with a larger
   (or overflowing) quotient, which the caller recomputes with
   rsimd_mod_f64(). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_vf64_mod(rsimd_vf64 x, rsimd_vf64 y, rsimd_mf64 *slow) {
  const rsimd_vf64 zero = rsimd_vf64_zero(), inf = rsimd_vf64_set1(HUGE_VAL);
  rsimd_vf64 k = rsimd_vf64_idiv(x, y, NULL);
  rsimd_vf64 r = rsimd_vf64_add(rsimd_vf64_fma(rsimd_vf64_neg(k), y, x), zero);
  rsimd_vf64 ay = rsimd_vf64_abs(y);
  rsimd_mf64 yinf = rsimd_mf64_and(rsimd_vf64_cmp_eq(ay, inf),
                                   rsimd_vf64_cmp_lt(rsimd_vf64_abs(x), inf));
  rsimd_mf64 differ =
    rsimd_mf64_or(rsimd_mf64_and(rsimd_vf64_cmp_lt(x, zero), rsimd_vf64_cmp_gt(y, zero)),
                  rsimd_mf64_and(rsimd_vf64_cmp_gt(x, zero), rsimd_vf64_cmp_lt(y, zero)));
  rsimd_mf64 finite = rsimd_mf64_and(rsimd_vf64_cmp_lt(rsimd_vf64_abs(x), inf),
                                     rsimd_vf64_cmp_lt(ay, inf));
  rsimd_mf64 small = rsimd_vf64_cmp_lt(rsimd_vf64_abs(rsimd_vf64_div(x, y)),
                                       rsimd_vf64_set1(4503599627370496.0));
  *slow = rsimd_mf64_andnot(small,
                            rsimd_mf64_and(finite, rsimd_vf64_cmp_ne(y, zero)));
  return rsimd_vf64_blend(r, rsimd_vf64_blend(x, y, differ), yinf);
}

/* Base R's pmin/pmax rules (rsimd_pmin_f64 ...), lane by lane: the
   layer's min(y, x) is y < x ? y : x, which already returns x when x is
   NaN or on a tie. num = 1 gives the _num forms. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_pminmax(rsimd_vf64 x, rsimd_vf64 y, int max, int num) {
  rsimd_vf64 r = max ? rsimd_vf64_max(y, x) : rsimd_vf64_min(y, x);
  return rsimd_vf64_blend(r, y, rsimd_vf64_is_nan(num ? x : y));
}

/* rsimd_copysign_f64() lane by lane. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_copysign(rsimd_vf64 x, rsimd_vf64 s) {
  const rsimd_vf64 sign = rsimd_vf64_set1(-0.0);
  rsimd_vf64 r = rsimd_vf64_or(rsimd_vf64_andnot(sign, x), rsimd_vf64_and(sign, s));
  return rsimd_vf64_blend(r, s, rsimd_vf64_is_nan(s));
}

/* rsimd_sign_f64() lane by lane. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sign(rsimd_vf64 x) {
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_vf64 r = rsimd_vf64_blend(x, rsimd_vf64_set1(1.0), rsimd_vf64_cmp_gt(x, zero));
  r = rsimd_vf64_blend(r, rsimd_vf64_set1(-1.0), rsimd_vf64_cmp_lt(x, zero));
  return rsimd_vf64_blend(r, zero, rsimd_vf64_cmp_eq(x, zero));
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
/* As rsimd_vf64_any_na() for an operand of doubles (i32 = 0) or int32
   elements (i32 = 1). */
RSIMD_ALWAYS_INLINE int rsimd_any_na_elt_(const void *p, const int i32, ptrdiff_t n) {
  ptrdiff_t i = 0;
  const int32_t *q = (const int32_t *) p;
  if (!i32) return rsimd_vf64_any_na((const double *) p, n);
  for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {
    if (rsimd_mi32_any(rsimd_vi32_is_na(rsimd_vi32_loadu(q + i)))) return 1;
  }
  for (; i < n; i++) {
    if (q[i] == RSIMD_NA_I32) return 1;
  }
  return 0;
}

/* Elements i, i + 1 ... of an operand of doubles (i32 = 0) or int32
   elements (i32 = 1, converted exactly; NA becomes NA_real_), as a full
   vector or as the tail of n - i < W elements with `fill` in the inactive
   lanes. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_vf64_load_elt_(const void *p, const int i32, ptrdiff_t i) {
  rsimd_vf64 v;
  if (!i32) return rsimd_vf64_loadu((const double *) p + i);
  v = rsimd_vf64_loadu_i32((const int32_t *) p + i);
  return rsimd_vf64_blend(v, rsimd_vf64_set1(rsimd_na_real()),
                          rsimd_vf64_cmp_eq(v, rsimd_vf64_set1(RSIMD_NA_I32_AS_F64)));
}
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_vf64_load_elt_tail_(const void *p, const int i32, ptrdiff_t i,
                                                        ptrdiff_t n, double fill) {
  double b[RSIMD_MAX_LANES_64] = {0.0}; /* every lane used is set; GCC cannot tell on SVE */
  ptrdiff_t j;
  if (!i32) return rsimd_vf64_loadu_p(rsimd_p64_while(i, n), (const double *) p + i, fill);
  for (j = 0; j < RSIMD_LANES_64; j++) b[j] = i + j < n ? rsimd_elt_f64(p, 1, i + j) : fill;
  return rsimd_vf64_loadu(b);
}

/* The vector folds. They mirror the scalar folds, with the term, the
   operand types and the na_check/na.rm choice as compile-time constants
   (each combination is a separate copy of the loop). Missing values are
   tracked as one running NaN mask. Without na.rm each block with a NaN is
   rescanned for NA, and at the end only the vectors after the last block
   are; with na.rm the chunk is rescanned for NA if it had a NaN. */

/* m = the missing lanes: those of vx, or for a pair term those where its
   term vt is NaN; recorded in mnan, and under narm counted as removed, in
   the lanes of rcnt. */
#define RSIMD_VFOLD_MISSING_(m, vx, vt)                                                  \
  do {                                                                                   \
    (m) = rsimd_vf64_is_nan(RSIMD_TERM_PAIR(term) ? (vt) : (vx));                        \
    mnan = rsimd_mf64_or(mnan, (m));                                                     \
    if (narm) rcnt = rsimd_vi64_inc(rcnt, rsimd_mf64_to_mi64(m));                        \
  } while (0)

/* t = term of vectors vx, vy; updates mnan and rcnt. */
#define RSIMD_VFOLD_TERM_(t, vx, vy)                                                     \
  do {                                                                                   \
    (t) = (vx);                                                                          \
    if (term == RSIMD_TERM_SQ) (t) = rsimd_vf64_mul((vx), (vx));                         \
    else if (term == RSIMD_TERM_ABS) (t) = rsimd_vf64_abs(vx);                           \
    else if (term == RSIMD_TERM_XY) (t) = rsimd_vf64_mul((vx), (vy));                    \
    else if (term == RSIMD_TERM_DEV) (t) = rsimd_vf64_sub((vx), vc);                     \
    else if (term == RSIMD_TERM_SQDIFF || term == RSIMD_TERM_DEVSQ) {                    \
      rsimd_vf64 d_ = rsimd_vf64_sub((vx), term == RSIMD_TERM_SQDIFF ? (vy) : vc);       \
      (t) = rsimd_vf64_mul(d_, d_);                                                      \
    }                                                                                    \
    if (check) {                                                                         \
      rsimd_mf64 m_;                                                                     \
      RSIMD_VFOLD_MISSING_(m_, vx, (t));                                                 \
      if (narm) (t) = rsimd_vf64_blend((t), rsimd_vf64_zero(), m_);                      \
    }                                                                                    \
  } while (0)

/* acc += term of vx, vy: with a fused multiply-add for a product term when
   use_fma is set (both factors of a removed lane are set to 0), else as
   the rounded term. */
#define RSIMD_VFOLD_ACC_(acc, vx, vy)                                                    \
  do {                                                                                   \
    if (use_fma && RSIMD_TERM_PRODUCT(term)) {                                           \
      rsimd_vf64 f_ = term == RSIMD_TERM_SQDIFF ? rsimd_vf64_sub((vx), (vy))             \
                      : term == RSIMD_TERM_DEVSQ ? rsimd_vf64_sub((vx), vc)              \
                                                 : (vx);                                 \
      rsimd_vf64 g_ = term == RSIMD_TERM_XY ? (vy) : f_;                                 \
      if (check) {                                                                       \
        rsimd_mf64 m_;                                                                   \
        /* The term of dist is NaN where its difference f_ is. */                       \
        RSIMD_VFOLD_MISSING_(m_, vx,                                                     \
                             term == RSIMD_TERM_XY ? rsimd_vf64_mul((vx), (vy)) : f_);   \
        if (narm) {                                                                      \
          f_ = rsimd_vf64_blend(f_, rsimd_vf64_zero(), m_);                              \
          g_ = rsimd_vf64_blend(g_, rsimd_vf64_zero(), m_);                              \
        }                                                                                \
      }                                                                                  \
      (acc) = rsimd_vf64_fma(f_, g_, (acc));                                             \
    } else {                                                                             \
      rsimd_vf64 t_;                                                                     \
      RSIMD_VFOLD_TERM_(t_, vx, vy);                                                     \
      (acc) = rsimd_vf64_add((acc), t_);                                                 \
    }                                                                                    \
  } while (0)

/* Loads of x and y at i: full vectors, or the tail, filled so that the
   inactive lanes' term is 0 (x = c for the terms with a constant, else
   0). */
#define RSIMD_VFOLD_LOAD_(vx, vy, i)                                                     \
  rsimd_vf64 vx = rsimd_vf64_load_elt_(x, xi32, (i));                                    \
  rsimd_vf64 vy = RSIMD_TERM_PAIR(term) ? rsimd_vf64_load_elt_(y, yi32, (i)) : vx
#define RSIMD_VFOLD_LOAD_TAIL_(vx, vy, i)                                                \
  rsimd_vf64 vx = rsimd_vf64_load_elt_tail_(                                             \
    x, xi32, (i), n, term == RSIMD_TERM_DEV || term == RSIMD_TERM_DEVSQ ? c : 0.0);      \
  rsimd_vf64 vy =                                                                        \
    RSIMD_TERM_PAIR(term) ? rsimd_vf64_load_elt_tail_(y, yi32, (i), n, 0.0) : vx

/* Element i of an operand of doubles or int32 elements, as a pointer. */
RSIMD_ALWAYS_INLINE const void *rsimd_elt_ptr_(const void *p, const int i32, ptrdiff_t i) {
  return i32 ? (const void *) ((const int32_t *) p + i) : (const void *) ((const double *) p + i);
}

/* For a fold without na.rm whose block [from, to) had a missing lane: 1,
   with r->saw_na and r->saw_nan set and the whole chunk of n counted, if
   the block has an NA, as the result is then NA whatever follows; else 0.
   Kept out of line: it runs at most once per block with a NaN. */
static RSIMD_NOINLINE int rsimd_vfold_block_na_(const void *x, int xi32, const void *y, int yi32,
                                                int pair, ptrdiff_t from, ptrdiff_t to,
                                                ptrdiff_t n, rsimd_reduce_result *r) {
  if (rsimd_any_na_elt_(rsimd_elt_ptr_(x, xi32, from), xi32, to - from) ||
      (pair && rsimd_any_na_elt_(rsimd_elt_ptr_(y, yi32, from), yi32, to - from))) {
    r->saw_na = 1;
    r->saw_nan = 1;
    r->count += n;
    return 1;
  }
  return 0;
}

/* The end of a fold of n elements with check: nan is set when a block
   rescan found only NaN, mnan holds the missing lanes since; the elements
   before from are known to hold no NA (the blocks rescanned without one
   or with no missing lane), so only [from, n) is scanned for an NA. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_done_(rsimd_mf64 mnan, int nan, ptrdiff_t from,
                                           ptrdiff_t removed, const void *x, const int xi32,
                                           const void *y, const int yi32, ptrdiff_t n, int term,
                                           rsimd_reduce_result *r) {
  if (nan) r->saw_nan = 1;
  if (rsimd_mf64_any(mnan)) {
    r->saw_nan = 1;
    if (!r->saw_na &&
        (rsimd_any_na_elt_(rsimd_elt_ptr_(x, xi32, from), xi32, n - from) ||
         (RSIMD_TERM_PAIR(term) &&
          rsimd_any_na_elt_(rsimd_elt_ptr_(y, yi32, from), yi32, n - from)))) {
      r->saw_na = 1;
    }
  }
  r->count += n - removed;
}

/* Fast mode: returns the sum of the n elements (see the rules above);
   product terms use fused multiply-adds when use_fma is set. */
RSIMD_ALWAYS_INLINE double rsimd_vfold_fast_(const void *x, const int xi32, const void *y,
                                             const int yi32, ptrdiff_t n, const int term,
                                             const int check, const int narm, const int use_fma,
                                             double c, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int stop = check && !narm;
  const rsimd_vf64 vc = rsimd_vf64_set1(c);
  rsimd_vf64 a0 = rsimd_vf64_zero(), a1 = a0, a2 = a0, a3 = a0;
  rsimd_mf64 mnan = rsimd_mf64_none();
  rsimd_vi64 rcnt = rsimd_vi64_zero();
  ptrdiff_t i = 0, from, tail;
  int nan = 0;
  while (i + 4 * W <= n) {
    const ptrdiff_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (from = i; i + 4 * W <= end; i += 4 * W) {
      RSIMD_VFOLD_LOAD_(x0, y0, i);
      RSIMD_VFOLD_LOAD_(x1, y1, i + W);
      RSIMD_VFOLD_LOAD_(x2, y2, i + 2 * W);
      RSIMD_VFOLD_LOAD_(x3, y3, i + 3 * W);
      RSIMD_VFOLD_ACC_(a0, x0, y0);
      RSIMD_VFOLD_ACC_(a1, x1, y1);
      RSIMD_VFOLD_ACC_(a2, x2, y2);
      RSIMD_VFOLD_ACC_(a3, x3, y3);
    }
    if (stop && rsimd_mf64_any(mnan)) {
      if (rsimd_vfold_block_na_(x, xi32, y, yi32, RSIMD_TERM_PAIR(term), from, i, n, r)) {
        return 0.0;
      }
      mnan = rsimd_mf64_none();
      nan = 1;
    }
  }
  /* Without stop, the only block is the whole chunk, and not rescanned. */
  tail = stop ? i : 0;
  for (; i + W <= n; i += W) {
    RSIMD_VFOLD_LOAD_(vx, vy, i);
    RSIMD_VFOLD_ACC_(a0, vx, vy);
  }
  if (i < n) {
    RSIMD_VFOLD_LOAD_TAIL_(vx, vy, i);
    RSIMD_VFOLD_ACC_(a0, vx, vy);
  }
  if (check) {
    rsimd_vfold_done_(mnan, nan, tail, rsimd_vi64_reduce_add(rcnt), x, xi32, y, yi32, n, term, r);
  } else {
    r->count += n;
  }
  return rsimd_vf64_sum_tree(rsimd_vf64_add(rsimd_vf64_add(a0, a1), rsimd_vf64_add(a2, a3)));
}

/* Compensated mode: Neumaier per lane of four accumulators, then the lanes
   (accumulator 0 lanes 0..W-1, then accumulator 1 ...) with the scalar
   scheme, then into r's running pair. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_comp_(const void *x, const int xi32, const void *y,
                                           const int yi32, ptrdiff_t n, const int term,
                                           const int check, const int narm, double c,
                                           rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 vc = rsimd_vf64_set1(c);
  rsimd_vf64 s0 = rsimd_vf64_zero(), s1 = s0, s2 = s0, s3 = s0;
  rsimd_vf64 c0 = s0, c1 = s0, c2 = s0, c3 = s0, t;
  const int stop = check && !narm;
  rsimd_mf64 mnan = rsimd_mf64_none();
  rsimd_vi64 rcnt = rsimd_vi64_zero();
  ptrdiff_t i = 0, j, from, tail;
  double sb[4][RSIMD_MAX_LANES_64], cb[4][RSIMD_MAX_LANES_64], S = 0.0, C = 0.0;
  int k, nan = 0;
  while (i + 4 * W <= n) {
    const ptrdiff_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (from = i; i + 4 * W <= end; i += 4 * W) {
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
    if (stop && rsimd_mf64_any(mnan)) {
      if (rsimd_vfold_block_na_(x, xi32, y, yi32, RSIMD_TERM_PAIR(term), from, i, n, r)) return;
      mnan = rsimd_mf64_none();
      nan = 1;
    }
  }
  tail = stop ? i : 0;
  for (; i + W <= n; i += W) {
    RSIMD_VFOLD_LOAD_(vx, vy, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    RSIMD_VF64_NEUMAIER(s0, c0, t);
  }
  if (i < n) {
    RSIMD_VFOLD_LOAD_TAIL_(vx, vy, i);
    RSIMD_VFOLD_TERM_(t, vx, vy);
    RSIMD_VF64_NEUMAIER(s0, c0, t);
  }
  if (check) {
    rsimd_vfold_done_(mnan, nan, tail, rsimd_vi64_reduce_add(rcnt), x, xi32, y, yi32, n, term, r);
  } else {
    r->count += n;
  }
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

RSIMD_ALWAYS_INLINE void rsimd_vfold_mode_(const void *x, const int xi32, const void *y,
                                           const int yi32, ptrdiff_t n, const int term,
                                           const int check, const int narm, double c,
                                           rsimd_reduce_result *r, int precision) {
  ptrdiff_t i;
  switch (precision) {
  case RSIMD_PREC_PAIRWISE:
    /* Leaves add the rounded terms (no fused multiply-add). */
    for (i = 0; i < n; i += RSIMD_PAIRWISE_LEAF) {
      ptrdiff_t len = n - i < RSIMD_PAIRWISE_LEAF ? n - i : RSIMD_PAIRWISE_LEAF;
      rsimd_pairwise_push(&r->pw, rsimd_vfold_fast_(rsimd_elt_ptr_(x, xi32, i), xi32,
                                                    rsimd_elt_ptr_(y, yi32, i), yi32, len, term,
                                                    check, narm, 0, c, r));
      /* Without na.rm, the sum is NA whatever follows. */
      if (check && !narm && r->saw_na) {
        r->count += n - i - len;
        break;
      }
    }
    break;
  case RSIMD_PREC_COMPENSATED:
    rsimd_vfold_comp_(x, xi32, y, yi32, n, term, check, narm, c, r);
    break;
  default:
    r->f64 += rsimd_vfold_fast_(x, xi32, y, yi32, n, term, check, narm, RSIMD_NATIVE_FMA, c, r);
    break;
  }
}

/* The vector form of rsimd_fold_gen(); term, xi32 and yi32 should be
   constants. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_gen(const void *x, const int xi32, const void *y,
                                         const int yi32, ptrdiff_t n, const int term, double c,
                                         rsimd_reduce_result *r, const rsimd_opts *o) {
  if (!RSIMD_TERM_PAIR(term)) y = x; /* never read */
  if (o->na_rm) rsimd_vfold_mode_(x, xi32, y, yi32, n, term, 1, 1, c, r, o->precision);
  else if (o->na_check) rsimd_vfold_mode_(x, xi32, y, yi32, n, term, 1, 0, c, r, o->precision);
  else rsimd_vfold_mode_(x, xi32, y, yi32, n, term, 0, 0, c, r, o->precision);
}
/* The vector form of rsimd_fold_f64_c(); `term` should be a constant. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_f64_c(const double *x, const double *y, ptrdiff_t n,
                                           const int term, double c, rsimd_reduce_result *r,
                                           const rsimd_opts *o) {
  rsimd_vfold_gen(x, 0, y, 0, n, term, c, r, o);
}
/* The vector form of rsimd_fold_f64(). */
RSIMD_ALWAYS_INLINE void rsimd_vfold_f64(const double *x, const double *y, ptrdiff_t n,
                                         const int term, rsimd_reduce_result *r,
                                         const rsimd_opts *o) {
  rsimd_vfold_gen(x, 0, y, 0, n, term, 0.0, r, o);
}
#endif /* RSIMD_NO_F64_SIMD */

/* Vectors per block of rsimd_vfold_sum_i32_(): each 32-bit lane of its
   accumulators then sums at most RSIMD_SUM_I32_BLOCK halves of at most
   2^16 in magnitude, far from overflow. */
#define RSIMD_SUM_I32_BLOCK 1024

/* The vector form of rsimd_fold_sum_i32(): 32-bit lanes. Each element is
   split into its high half (x >> 16, signed; logical for absval, whose |x|
   can be 2^31) and its low half (x & 0xffff), summed in separate 32-bit
   accumulators over blocks of RSIMD_SUM_I32_BLOCK vectors and combined
   exactly in 64 bits as hi * 2^16 + lo. Without na.rm a block with an NA
   ends the fold, as the sum is then NA whatever follows. The halves of two
   vectors are added together before they reach the accumulators: gcc 16
   folds an accumulator update per vector into ssra on neon and then
   copies the accumulators with a mov around each one in the loops without
   the NA compare, which made na_check = FALSE slower than the check. */
RSIMD_ALWAYS_INLINE void rsimd_vfold_sum_i32_(const int32_t *x, ptrdiff_t n, const int absval,
                                              const int check, const int narm,
                                              rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const rsimd_vi32 na = rsimd_vi32_set1(RSIMD_NA_I32), zero = rsimd_vi32_zero();
  const rsimd_vi32 low = rsimd_vi32_set1(0xffff);
  int64_t s = 0;
  ptrdiff_t i = 0, removed = 0;
/* Sets h and l to the halves of v, after the NA handling and absval. Under
   na.rm the NAs are only counted: a nonzero count is how they are seen. */
#define RSIMD_VFOLD_I32_(h, l, v)                                                        \
  do {                                                                                   \
    if (check) {                                                                         \
      rsimd_mi32 m_ = rsimd_vi32_cmp_eq((v), na);                                        \
      if (narm) {                                                                        \
        (v) = rsimd_vi32_blend((v), zero, m_);                                           \
        cnt = rsimd_vi32_inc(cnt, m_);                                                   \
      } else {                                                                           \
        mna = rsimd_mi32_or(mna, m_);                                                    \
      }                                                                                  \
    }                                                                                    \
    if (absval) {                                                                        \
      (v) = rsimd_vi32_blend((v), rsimd_vi32_sub(zero, (v)), rsimd_vi32_cmp_lt((v), zero)); \
      (h) = rsimd_vi32_srl((v), 16);                                                     \
    } else {                                                                             \
      (h) = rsimd_vi32_sra((v), 16);                                                     \
    }                                                                                    \
    (l) = rsimd_vi32_and((v), low);                                                      \
  } while (0)
  while (i < n) {
    const ptrdiff_t end = n - i > RSIMD_SUM_I32_BLOCK * W ? i + RSIMD_SUM_I32_BLOCK * W : n;
    rsimd_vi32 hi = zero, lo = zero, cnt = zero, h0, h1, l0, l1;
    rsimd_mi32 mna = rsimd_mi32_none();
    for (; i + 2 * W <= end; i += 2 * W) {
      rsimd_vi32 v0 = rsimd_vi32_loadu(x + i), v1 = rsimd_vi32_loadu(x + i + W);
      RSIMD_VFOLD_I32_(h0, l0, v0);
      RSIMD_VFOLD_I32_(h1, l1, v1);
      hi = rsimd_vi32_add(hi, rsimd_vi32_add(h0, h1));
      lo = rsimd_vi32_add(lo, rsimd_vi32_add(l0, l1));
    }
    for (; i < end; i += W) {
      rsimd_vi32 v = rsimd_vi32_loadu_p(rsimd_p32_while(i, end), x + i, 0);
      RSIMD_VFOLD_I32_(h0, l0, v);
      hi = rsimd_vi32_add(hi, h0);
      lo = rsimd_vi32_add(lo, l0);
    }
    i = end;
    if (narm) {
      const ptrdiff_t k = (ptrdiff_t) rsimd_vi32_reduce_add(cnt);
      if (k) r->saw_na = 1;
      removed += k;
    } else if (check && rsimd_mi32_any(mna)) {
      r->saw_na = 1;
      r->count += n;
      return;
    }
    s += rsimd_vi32_reduce_add(hi) * 65536 + rsimd_vi32_reduce_add(lo);
  }
#undef RSIMD_VFOLD_I32_
  rsimd_merge_i64_sum(r, s);
  r->count += n - removed;
}
RSIMD_ALWAYS_INLINE void rsimd_vfold_sum_i32(const int32_t *x, ptrdiff_t n, const int absval,
                                             rsimd_reduce_result *r, const rsimd_opts *o) {
  if (o->na_rm) rsimd_vfold_sum_i32_(x, n, absval, 1, 1, r);
  else if (o->na_check) rsimd_vfold_sum_i32_(x, n, absval, 1, 0, r);
  else rsimd_vfold_sum_i32_(x, n, absval, 0, 0, r);
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
