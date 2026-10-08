/* .Call entry points of the reductions: classify the input, fold it chunk
   by chunk through the active implementation, and finish with base R's
   result types (rvec.h). Empty inputs run no kernel: the identity of the
   reduction is finished directly. */

#include <math.h>
#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "api_complex.h"

static void bad_type(rsimd_etype type) {
  Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[type]);
}

static int is_numeric(rsimd_etype t) {
  return t == RSIMD_F64 || t == RSIMD_I32 || t == RSIMD_LGL;
}

/* Complex numbers per block of Mod values (8 KiB of doubles on the
   stack). */
#define RSIMD_MOD_BLOCK 1024

/* In a chunk loop of a reduction whose result is NA once an NA is seen
   without na.rm: stop reading. */
#define RSIMD_STOP_AT_NA(r) \
  if ((r).saw_na && !o.na_rm) break

/* RSIMD_MATH_FAST for accuracy code 1 (option rsimd.math_accuracy =
   "fast"), 0 for 0; NULL reads the option. */
static int accuracy_bit(SEXP accuracy) {
  return rsimd_arg_accuracy(accuracy) ? RSIMD_MATH_FAST : 0;
}

/* sum(x, na.rm): double for double x; integer for integer and logical x,
   or double when the total does not fit in an integer; integer64 for
   integer64 x, exact, or NA with bit64's warning when the total does not
   fit; complex for complex x, each part summed as a double vector would be
   (with na.rm an element goes when either part is NA or NaN). precision is
   the integer code RSIMD_PREC_*. */
static SEXP simd_sum_impl(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->sum_f64(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->sum_i32(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        {
                          rsimd_active->sum_i64((const int64_t *) px, len, &r, &o);
                          RSIMD_STOP_AT_NA(r);
                        });
    break;
  case RSIMD_C128: return rsimd_c128_sum(&in, &o);
  default: bad_type(in.type);
  }
  return rsimd_reduce_finish(RSIMD_RED_SUM, in.type, in.n, &r, &o);
}

SEXP C_simd_sum(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_sum_impl(x, na_rm, na_check, precision));
}

/* ---- The careful product ----------------------------------------------
   The kernels keep several partial products in double, so a partial
   product can overflow to Inf or underflow to 0 although the whole product
   is finite, and Inf times a zero factor in another lane is NaN: the
   result then depends on the tier and on where each element lands. When
   the fast product is 0, +-Inf or NaN with no missing value seen, the
   entry points redo it here in order, the same on every tier: zeros and
   infinities are counted apart (a zero gives a zero and an infinity an
   infinity, with the sign of the product, both NaN), and the other
   factors are multiplied with their binary exponents kept apart in e, so
   no partial product leaves the double range. The magnitudes multiplied
   stay in [2^-511, 2^511], so every product is a normal double and rounds
   as the unbounded product would. */

/* The range of the magnitudes multiplied. */
#define CPROD_HI 0x1p511
#define CPROD_LO 0x1p-511

typedef struct {
  double m;     /* the magnitude is m * 2^e */
  int64_t e;
  int neg;      /* the sign bits of the factors, xored */
  int zero, inf, nan;
} cprod_state;

static void cprod_init(cprod_state *s) {
  s->m = 1.0;
  s->e = 0;
  s->neg = s->zero = s->inf = s->nan = 0;
}

/* Multiplies s by v; NaN is skipped under narm. */
static inline void cprod_add(cprod_state *s, double v, int narm) {
  double a;
  int k;
  if (isnan(v)) {
    if (!narm) s->nan = 1;
    return;
  }
  s->neg ^= signbit(v) != 0;
  a = fabs(v);
  if (a == 0.0) s->zero = 1;
  else if (isinf(a)) s->inf = 1;
  else {
    if (a > CPROD_HI || a < CPROD_LO) {
      a = frexp(a, &k);
      s->e += k;
    }
    s->m *= a;
    if (s->m > CPROD_HI || s->m < CPROD_LO) {
      s->m = frexp(s->m, &k);
      s->e += k;
    }
  }
}

/* The product of s, or `fast` when s met a NaN (with na_check = FALSE and
   no na.rm, where the fast product stands). */
static double cprod_value(const cprod_state *s, double fast) {
  double v;
  if (s->nan) return fast;
  if (s->zero && s->inf) return R_NaN;
  if (s->zero) v = 0.0;
  else if (s->inf || s->e > 4096) v = R_PosInf;
  else if (s->e < -4096) v = 0.0;
  else v = ldexp(s->m, (int) s->e);
  return s->neg ? -v : v;
}

/* Whether the fast product in r must be redone: 0, +-Inf or NaN, with no
   missing value that decides the result. */
static int cprod_needed(const rsimd_reduce_result *r, const rsimd_opts *o) {
  if (!o->na_rm && (r->saw_na || r->saw_nan)) return 0;
  return r->f64 == 0.0 || !R_FINITE(r->f64);
}

/* prod(x, na.rm): double, or complex for complex x. The precision mode
   applies to complex x only. */
static SEXP simd_prod_impl(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  rsimd_reduce_result_init(&r, RSIMD_RED_PROD);
  switch (in.type) {
  case RSIMD_C128: return rsimd_c128_prod(&in, &o);
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->prod_f64(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->prod_i32(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  default: bad_type(in.type);
  }
  if (cprod_needed(&r, &o)) {
    cprod_state s;
    cprod_init(&s);
    if (in.type == RSIMD_F64) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
        for (R_xlen_t i = 0; i < len; i++) cprod_add(&s, px[i], o.na_rm);
      });
    } else {
      /* Without na.rm an NA here was promised away (na_check = FALSE), and
         the kernels multiplied it as RSIMD_NA_I32_AS_F64. */
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        for (R_xlen_t i = 0; i < len; i++) {
          if (px[i] != NA_INTEGER || !o.na_rm) cprod_add(&s, (double) px[i], 0);
        }
      });
    }
    r.f64 = cprod_value(&s, r.f64);
  }
  return rsimd_reduce_finish(RSIMD_RED_PROD, in.type, in.n, &r, &o);
}

SEXP C_simd_prod(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_prod_impl(x, na_rm, na_check, precision));
}

/* prod(x + y) (op 0) or prod(x - y) (op 1) without the vector of sums,
   with the length-1 broadcast rule: double for double, integer and logical
   operands, whose sums are computed in double (so integer sums do not
   overflow), and complex when both are complex (the R side converts the
   other), as simd_prod of the sums. A pair whose sum is NaN is missing
   (removed under na.rm), so the result is prod(x + y, na.rm) of those
   double sums. */
static SEXP simd_prod2_impl(SEXP x, SEXP y, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_bin b;
  int ew = rsimd_arg_int1(op, "op") ? RSIMD_EW_SUB : RSIMD_EW_ADD, flags;

  rsimd_bin_init(&b, x, y);
  rsimd_opts_init(&o, na_rm, na_check, b.x.no_na_hint && b.y.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  if (b.x.type == RSIMD_C128 && b.y.type == RSIMD_C128) return rsimd_c128_prod2(ew, &b, &o);
  if (!is_numeric(b.x.type)) bad_type(b.x.type);
  if (!is_numeric(b.y.type)) bad_type(b.y.type);
  flags = (b.x_scalar ? RSIMD_EW_SCALAR(0) : 0) | (b.y_scalar ? RSIMD_EW_SCALAR(1) : 0) |
          (b.x.type != RSIMD_F64 ? RSIMD_EW_I32(0) : 0) |
          (b.y.type != RSIMD_F64 ? RSIMD_EW_I32(1) : 0);
  rsimd_reduce_result_init(&r, RSIMD_RED_PROD);
#define RSIMD_PROD2_LOOP_(TX, TY, ...)                                                   \
  RSIMD_FOREACH_CHUNK2T(&b, TX, TY, px, py, len, off, __VA_ARGS__)
#define RSIMD_PROD2_TYPES_(...)                                                          \
  if (b.x.type == RSIMD_F64) {                                                           \
    if (b.y.type == RSIMD_F64) RSIMD_PROD2_LOOP_(double, double, __VA_ARGS__);           \
    else RSIMD_PROD2_LOOP_(double, int, __VA_ARGS__);                                    \
  } else {                                                                               \
    if (b.y.type == RSIMD_F64) RSIMD_PROD2_LOOP_(int, double, __VA_ARGS__);              \
    else RSIMD_PROD2_LOOP_(int, int, __VA_ARGS__);                                       \
  }
  RSIMD_PROD2_TYPES_({ rsimd_active->prod2_f64(ew, px, py, len, flags, &r, &o); })
  if (cprod_needed(&r, &o)) {
    /* The pairs as the kernels form them: x - y as x + (-y), and an NA
       operand (double or integer) giving a NaN sum. */
    const int xi = b.x.type != RSIMD_F64, yi = b.y.type != RSIMD_F64;
    const ptrdiff_t xs = !b.x_scalar, ys = !b.y_scalar;
    cprod_state s;
    cprod_init(&s);
    RSIMD_PROD2_TYPES_({
      for (R_xlen_t i = 0; i < len; i++) {
        double u = rsimd_elt_f64(px, xi, i * xs), v = rsimd_elt_f64(py, yi, i * ys);
        cprod_add(&s, ew == RSIMD_EW_SUB ? u - v : u + v, o.na_rm);
      }
    })
    r.f64 = cprod_value(&s, r.f64);
  }
#undef RSIMD_PROD2_TYPES_
#undef RSIMD_PROD2_LOOP_
  return rsimd_reduce_finish(RSIMD_RED_PROD, RSIMD_F64, b.n, &r, &o);
}

SEXP C_simd_prod2(SEXP x, SEXP y, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_prod2_impl(x, y, op, na_rm, na_check, precision));
}

/* The mean of the double x (r->count elements left) when its sum overflowed
   to +-Inf or NaN: that of its infinite elements when it has some (Inf,
   -Inf, or NaN for both), else the compensated sum (s, c) of x scaled by
   2^-k, with 2^k >= 2 r->count so that no partial sum can overflow,
   divided by the count as a pair (so rep(y, n) gives y) and scaled back.
   Rare, so scalar. */
static double mean_unscaled_f64(const rsimd_in *in, const rsimd_opts *o,
                                const rsimd_reduce_result *r) {
  double s = 0.0, c = 0.0, n = (double) r->count, q;
  int k = 1, pinf = 0, ninf = 0;
  while (k < 62 && ldexp(1.0, k - 1) < n) k++;
  RSIMD_FOREACH_CHUNK(in, double, px, len, off, {
    for (R_xlen_t i = 0; i < len; i++) {
      double x = px[i];
      if (isnan(x)) {
        if (o->na_rm) continue;
        return x;
      }
      if (isinf(x)) {
        if (x > 0) pinf = 1;
        else ninf = 1;
      } else {
        rsimd_neumaier_add(&s, &c, ldexp(x, -k));
      }
    }
  });
  if (pinf || ninf) return pinf && ninf ? R_NaN : pinf ? R_PosInf : R_NegInf;
  q = s / n;
  return ldexp(q + (fma(-q, n, s) + c) / n, k);
}

/* The mean of x, from the sum fold r (sum_f64 or sum_i32 over every chunk,
   with r->count > 0 elements left): the exact integer sum divided in long
   double, or the double sum divided by the count refined by the mean of
   the deviations from it, in pairwise and compensated modes or in every
   mode when refine is set; mean_unscaled_f64() when the sum overflowed. */
static double mean_value(const rsimd_in *in, const rsimd_opts *o, const rsimd_reduce_result *r,
                         int refine) {
  double m;
  if (in->type != RSIMD_F64) {
    long double s = (long double) r->i64;
    if (r->overflow) s += (long double) r->f64;
    return (double) (s / (long double) r->count);
  }
  m = rsimd_reduce_value(r, o->precision) / (double) r->count;
  /* Rescaled, the mean needs no refinement, whose deviations from m near
     the edge of the range would lose m's low bits or overflow. */
  if (!R_FINITE(m)) return mean_unscaled_f64(in, o, r);
  if (refine || o->precision != RSIMD_PREC_FAST) {
    rsimd_reduce_result d;
    rsimd_reduce_result_init(&d, RSIMD_RED_SUM);
    RSIMD_FOREACH_CHUNK(in, double, px, len, off,
                        { rsimd_active->sum_dev_f64(px, len, m, &d, o); });
    m += rsimd_reduce_value(&d, o->precision) / (double) r->count;
  }
  return m;
}

/* mean(x, na.rm): double, or complex for complex x (each part's sum, as
   sum gives it, divided by the count). Doubles: the sum divided by the count; in
   pairwise and compensated modes refined, as base R does, by adding the
   mean of the deviations from it (a second pass) when it is finite.
   Integers and logicals: the exact sum divided in long double. */
static SEXP simd_mean_impl(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  double m;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  rsimd_reduce_result_init(&r, RSIMD_RED_MEAN);
  switch (in.type) {
  case RSIMD_C128: return rsimd_c128_mean(&in, &o);
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->sum_f64(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->sum_i32(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  default: bad_type(in.type);
  }
  /* Missing values and empty input. */
  if (r.count == 0 || (!o.na_rm && (r.saw_na || r.saw_nan))) {
    return rsimd_reduce_finish(RSIMD_RED_MEAN, in.type, in.n, &r, &o);
  }
  m = mean_value(&in, &o, &r, 0);
  return Rf_ScalarReal(m);
}

SEXP C_simd_mean(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_mean_impl(x, na_rm, na_check, precision));
}

/* Folds Mod(x) of a complex input into the min/max result r: Mod of each
   block (math1_c128, as simd_abs computes it, in fast mode with `fast`)
   into a buffer that the double kernel reads. */
static void minmax_mod(const rsimd_in *in, int fast, rsimd_reduce_result *r, const rsimd_opts *o) {
  double buf[RSIMD_MOD_BLOCK];
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i += RSIMD_MOD_BLOCK) {
      R_xlen_t l = len - i < RSIMD_MOD_BLOCK ? len - i : RSIMD_MOD_BLOCK;
      rsimd_active->math1_c128(RSIMD_CMATH_MOD | fast, px + i, l, buf);
      rsimd_active->minmax_f64(buf, l, 0, r, o);
    }
  });
}

/* min (op 0), max (op 1) or range (op 2) of x, or with absval of abs(x)
   (Mod(x) for complex x, a double result, in the accuracy mode
   `accuracy`): both extrema come from one pass. Integer and logical
   results are integer, except for empty input (after na.rm), which gives
   Inf and -Inf with base R's warnings. integer64 results are integer64;
   empty input gives +INT64_MAX and -INT64_MAX with bit64's warnings (one
   for range). */
static SEXP simd_minmax_impl(SEXP x, SEXP op, SEXP na_rm, SEXP na_check, SEXP absval,
                             SEXP accuracy) {
  rsimd_reduce_result r, rmax;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op"), ab = rsimd_arg_lgl1(absval, "absval");
  SEXP lo, hi, out;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.extrema = which == 0 ? RSIMD_EXT_MIN : which == 1 ? RSIMD_EXT_MAX : RSIMD_EXT_BOTH;
  rsimd_reduce_result_init(&r, RSIMD_RED_MIN);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        {
                          rsimd_active->minmax_f64(px, len, ab, &r, &o);
                          RSIMD_STOP_AT_NA(r);
                        });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->minmax_i32(px, len, ab, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        {
                          rsimd_active->minmax_i64((const int64_t *) px, len, ab, &r, &o);
                          RSIMD_STOP_AT_NA(r);
                        });
    break;
  case RSIMD_C128:
    if (!ab) bad_type(in.type);
    minmax_mod(&in, accuracy_bit(accuracy), &r, &o);
    in.type = RSIMD_F64;
    break;
  default: bad_type(in.type);
  }
  rmax = r;
  rmax.f64 = r.f64_hi;
  rmax.i64 = r.i64_hi;
  if (which == 0) return rsimd_reduce_finish(RSIMD_RED_MIN, in.type, in.n, &r, &o);
  if (which == 1) return rsimd_reduce_finish(RSIMD_RED_MAX, in.type, in.n, &rmax, &o);
  if (in.type == RSIMD_I64) {
    /* Both extrema (NA, or the empty-input pair with one warning) as
       integer64. */
    int64_t pair[2];
    int missing = !o.na_rm && r.saw_na;
    if (!missing && r.count == 0) {
      rsimd_warn("no non-NA value, returning c(+9223372036854775807, -9223372036854775807)");
      pair[0] = INT64_MAX;
      pair[1] = -INT64_MAX;
    } else {
      pair[0] = missing ? RSIMD_NA_I64 : r.i64;
      pair[1] = missing ? RSIMD_NA_I64 : r.i64_hi;
    }
    out = PROTECT(rsimd_alloc_like(RSIMD_I64, 2));
    memcpy(rsimd_out_ptr(out), pair, sizeof pair);
    rsimd_set_i64_class(out);
    UNPROTECT(1);
    return out;
  }
  lo = PROTECT(rsimd_reduce_finish(RSIMD_RED_MIN, in.type, in.n, &r, &o));
  hi = PROTECT(rsimd_reduce_finish(RSIMD_RED_MAX, in.type, in.n, &rmax, &o));
  out = rsimd_range_pair(lo, hi);
  UNPROTECT(2);
  return out;
}

SEXP C_simd_minmax(SEXP x, SEXP op, SEXP na_rm, SEXP na_check, SEXP absval, SEXP accuracy) {
  rsimd_entry();
  return rsimd_exit(simd_minmax_impl(x, op, na_rm, na_check, absval, accuracy));
}

/* which.min (dir = 0) or which.max (dir = 1) of Mod(x) for complex x, into
   r->idx: Mod of each block into a buffer, its extremum from the double
   kernel, and the block's first element equal to it when it beats the
   extremum of the earlier blocks (which win ties). */
static void which_mod(const rsimd_in *in, int dir, int fast, rsimd_reduce_result *r,
                      const rsimd_opts *o) {
  double buf[RSIMD_MOD_BLOCK], best = 0;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i += RSIMD_MOD_BLOCK) {
      R_xlen_t l = len - i < RSIMD_MOD_BLOCK ? len - i : RSIMD_MOD_BLOCK;
      rsimd_reduce_result rb;
      double v;
      rsimd_active->math1_c128(RSIMD_CMATH_MOD | fast, px + i, l, buf);
      rsimd_reduce_result_init(&rb, RSIMD_RED_MIN);
      rsimd_active->minmax_f64(buf, l, 0, &rb, o);
      if (rb.count == 0) continue;
      v = dir ? rb.f64_hi : rb.f64;
      if (r->idx < 0 || (dir ? v > best : v < best)) {
        best = v;
        r->idx = off + i + rsimd_active->find_f64(buf, l, 0, v);
      }
    }
  });
}

/* which.min (max = FALSE) or which.max (max = TRUE), or with absval of
   abs(x) (Mod(x) for complex x, in the accuracy mode `accuracy`): missing
   values are ignored. The extremum comes from the min/max kernels, then
   the first element equal to it is looked up (so 0 and -0 tie, as in base
   R). Raw vectors use a one-pass kernel. */
static SEXP simd_which_impl(SEXP x, SEXP max, SEXP absval, SEXP accuracy) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int dir = rsimd_arg_lgl1(max, "max"), ab = rsimd_arg_lgl1(absval, "absval");
  int op = dir ? RSIMD_RED_WHICH_MAX : RSIMD_RED_WHICH_MIN;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, Rf_ScalarLogical(TRUE), in.no_na_hint);
  o.na_rm = 1;
  o.extrema = dir ? RSIMD_EXT_MAX : RSIMD_EXT_MIN;
  rsimd_reduce_result_init(&r, op);
  switch (in.type) {
  case RSIMD_F64: {
    double v;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        { rsimd_active->minmax_f64(px, len, ab, &r, &o); });
    if (r.count == 0) break;
    v = dir ? r.f64_hi : r.f64;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      R_xlen_t k = rsimd_active->find_f64(px, len, ab, v);
      if (k >= 0) {
        r.idx = off + k;
        break;
      }
    });
    break;
  }
  case RSIMD_I32:
  case RSIMD_LGL: {
    int v;
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->minmax_i32(px, len, ab, &r, &o); });
    if (r.count == 0) break;
    v = (int) (dir ? r.i64_hi : r.i64);
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      R_xlen_t k = rsimd_active->find_i32(px, len, ab, v);
      if (k >= 0) {
        r.idx = off + k;
        break;
      }
    });
    break;
  }
  case RSIMD_I64: {
    int64_t v;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        { rsimd_active->minmax_i64((const int64_t *) px, len, ab, &r, &o); });
    if (r.count == 0) break;
    v = dir ? r.i64_hi : r.i64;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      R_xlen_t k = rsimd_active->find_i64((const int64_t *) px, len, ab, v);
      if (k >= 0) {
        r.idx = off + k;
        break;
      }
    });
    break;
  }
  case RSIMD_U8:
    if (ab) bad_type(in.type);
    RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off,
                        { rsimd_active->which_u8(px, len, off, dir, &r); });
    break;
  case RSIMD_C128:
    if (!ab) bad_type(in.type);
    which_mod(&in, dir, accuracy_bit(accuracy), &r, &o);
    in.type = RSIMD_F64;
    break;
  default: bad_type(in.type);
  }
  return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
}

SEXP C_simd_which(SEXP x, SEXP max, SEXP absval, SEXP accuracy) {
  rsimd_entry();
  return rsimd_exit(simd_which_impl(x, max, absval, accuracy));
}

/* any (all = FALSE) or all (all = TRUE) with three-valued logic. Double and
   raw inputs are read as logical values without being converted, with
   base R's warning; integer64 inputs (non-zero TRUE) without one, as in
   bit64; reading stops once the answer is known. */
static SEXP simd_anyall_impl(SEXP x, SEXP all, SEXP na_rm) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int is_all = rsimd_arg_lgl1(all, "all");
  int op = is_all ? RSIMD_RED_ALL : RSIMD_RED_ANY;
  int stop = is_all ? RSIMD_STOP_FALSE : RSIMD_STOP_TRUE;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, Rf_ScalarLogical(TRUE), in.no_na_hint);
  rsimd_reduce_result_init(&r, op);
#define RSIMD_DONE_ (is_all ? r.any_false : r.any_true)
  switch (in.type) {
  case RSIMD_F64:
  case RSIMD_U8:
    if (in.n > 0) {
      rsimd_warn("coercing argument of type '%s' to logical",
                 in.type == RSIMD_F64 ? "double" : "raw");
    }
    if (in.type == RSIMD_U8) {
      RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off, {
        rsimd_active->anyall_u8(px, len, stop, &r);
        if (RSIMD_DONE_) break;
      });
    } else {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
        rsimd_active->anyall_f64(px, len, stop, &r, &o);
        if (RSIMD_DONE_) break;
      });
    }
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->anyall_i32(px, len, stop, &r, &o);
      if (RSIMD_DONE_) break;
    });
    break;
  case RSIMD_I64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->anyall_i64((const int64_t *) px, len, stop, &r, &o);
      if (RSIMD_DONE_) break;
    });
    break;
  default: bad_type(in.type);
  }
#undef RSIMD_DONE_
  return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
}

SEXP C_simd_anyall(SEXP x, SEXP all, SEXP na_rm) {
  rsimd_entry();
  return rsimd_exit(simd_anyall_impl(x, all, na_rm));
}

/* Runs the missing-value kernel of x's type over every chunk in `mode`
   (RSIMD_NAMODE_*), stopping after the first missing value for ANY. */
static void na_scan(rsimd_in *in, int mode, void *out, rsimd_reduce_result *r) {
  switch (in->type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(in, double, px, len, off, {
      rsimd_active->na_f64(px, len, mode, off, out, r);
      if (r->saw_na) break;
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(in, int, px, len, off, {
      rsimd_active->na_i32(px, len, mode, off, out, r);
      if (r->saw_na) break;
    });
    break;
  case RSIMD_C128:
    RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
      rsimd_active->na_c128(px, len, mode, off, out, r);
      if (r->saw_na) break;
    });
    break;
  case RSIMD_I64:
    RSIMD_FOREACH_CHUNK(in, double, px, len, off, {
      rsimd_active->na_i64((const int64_t *) px, len, mode, off, out, r);
      if (r->saw_na) break;
    });
    break;
  default: bad_type(in->type);
  }
}

/* The 1-based indices of the `count` elements that the missing-value
   kernels find in mode `flag` (0, or RSIMD_NAMODE_TRUE): integer, or
   double when one of them exceeds INT_MAX. They are written as doubles
   for a long vector, and converted to integer if the last fits after
   all. */
static SEXP which_indices(rsimd_in *in, int flag, R_xlen_t count) {
  rsimd_reduce_result r;
  rsimd_etype t = in->n > INT_MAX && count > 0 ? RSIMD_F64 : RSIMD_I32;
  R_xlen_t i;
  SEXP out = PROTECT(rsimd_alloc_like(t, count));
  rsimd_reduce_result_init(&r, RSIMD_RED_COUNT_NA);
  if (count > 0) {
    na_scan(in, (t == RSIMD_F64 ? RSIMD_NAMODE_WHICH_F64 : RSIMD_NAMODE_WHICH_I32) | flag,
            rsimd_out_ptr(out), &r);
  }
  if (t == RSIMD_F64 && ((const double *) rsimd_out_ptr(out))[count - 1] <= INT_MAX) {
    const double *src = (const double *) rsimd_out_ptr(out);
    SEXP small = PROTECT(rsimd_alloc_like(RSIMD_I32, count));
    int *dst = (int *) rsimd_out_ptr(small);
    for (i = 0; i < count; i++) dst[i] = (int) src[i];
    UNPROTECT(2);
    return small;
  }
  UNPROTECT(1);
  return out;
}

/* any_na (mode 0, logical), count_na (mode 1, double) or which_na (mode 2,
   the 1-based indices of the missing elements: integer, or double when one
   of them exceeds INT_MAX). NaN counts as missing, and so does an integer64
   NA. Raw vectors and inputs R knows to be NA-free are not read. */
static SEXP simd_na_impl(SEXP x, SEXP mode) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int m = rsimd_arg_int1(mode, "mode");

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, Rf_ScalarLogical(TRUE), in.no_na_hint);
  rsimd_reduce_result_init(&r, m == 0 ? RSIMD_RED_ANY_NA : RSIMD_RED_COUNT_NA);
  if (!in.no_na_hint) na_scan(&in, m == 0 ? RSIMD_NAMODE_ANY : RSIMD_NAMODE_COUNT, NULL, &r);
  if (m == 0) return rsimd_reduce_finish(RSIMD_RED_ANY_NA, in.type, in.n, &r, &o);
  if (m == 1) return rsimd_reduce_finish(RSIMD_RED_COUNT_NA, in.type, in.n, &r, &o);
  return which_indices(&in, 0, (R_xlen_t) r.i64);
}

SEXP C_simd_na(SEXP x, SEXP mode) {
  rsimd_entry();
  return rsimd_exit(simd_na_impl(x, mode));
}

/* sum(x) (count, mode 0) or which(x) (mode 1) of a logical x: the number
   of TRUE elements, integer or double beyond INT_MAX, and NA_integer_ when
   an element is NA and na.rm is FALSE; or their 1-based indices, NA
   ignored, as which_na() gives indices. An input known to be NA-free is
   not scanned for NA. */
static SEXP simd_true_impl(SEXP x, SEXP mode, SEXP na_rm) {
  rsimd_reduce_result r;
  rsimd_in in;
  int m = rsimd_arg_int1(mode, "mode"), narm = rsimd_arg_lgl1(na_rm, "na.rm");
  R_xlen_t count;

  rsimd_in_init(&in, x, "x");
  if (in.type != RSIMD_LGL) bad_type(in.type);
  rsimd_reduce_result_init(&r, RSIMD_RED_COUNT_NA);
  if (m == 0 && !narm && !in.no_na_hint) {
    na_scan(&in, RSIMD_NAMODE_ANY, NULL, &r);
    if (r.saw_na) return Rf_ScalarInteger(NA_INTEGER);
  }
  na_scan(&in, RSIMD_NAMODE_COUNT | RSIMD_NAMODE_TRUE, NULL, &r);
  count = (R_xlen_t) r.i64;
  if (m == 1) return which_indices(&in, RSIMD_NAMODE_TRUE, count);
  return count > INT_MAX ? Rf_ScalarReal((double) count) : Rf_ScalarInteger((int) count);
}

SEXP C_simd_true(SEXP x, SEXP mode, SEXP na_rm) {
  rsimd_entry();
  return rsimd_exit(simd_true_impl(x, mode, na_rm));
}

/* ---- Sums of squares, norms, distances, variance -------------------------- */

/* sum_sq (op 0), norm (op 1) or sum_abs (op 2) of x. sum_abs is a sum of
   |x| with simd_sum's result types (integer for integer and logical x,
   double when the total does not fit); the others are double. */
static SEXP simd_sum_sq_impl(SEXP x, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op");

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    if (which == 2) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                          {
                            rsimd_active->sumabs_f64(px, len, &r, &o);
                            RSIMD_STOP_AT_NA(r);
                          });
    } else {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
        rsimd_active->sumsq_f64(px, len, &r, &o);
        RSIMD_STOP_AT_NA(r);
      });
    }
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    if (which == 2) {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        rsimd_active->sumabs_i32(px, len, &r, &o);
        RSIMD_STOP_AT_NA(r);
      });
    } else {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        rsimd_active->sumsq_i32(px, len, &r, &o);
        RSIMD_STOP_AT_NA(r);
      });
    }
    break;
  default: bad_type(in.type);
  }
  if (which == 2) return rsimd_reduce_finish(RSIMD_RED_SUM, in.type, in.n, &r, &o);
  if (which == 1) {
    rsimd_reduce_settle(&r, o.precision);
    r.f64 = sqrt(r.f64);
    return rsimd_reduce_finish(RSIMD_RED_NORM, in.type, in.n, &r, &o);
  }
  return rsimd_reduce_finish(RSIMD_RED_SUM_SQ, in.type, in.n, &r, &o);
}

SEXP C_simd_sum_sq(SEXP x, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_sum_sq_impl(x, op, na_rm, na_check, precision));
}

/* dot (op 0), dist (op 1) or cosine (op 2) of x and y, which must have the
   same length (no broadcast). Double, integer and logical operands mix
   without conversion: the kernels read int32 elements as doubles. A pair
   whose term is NaN (a missing element, Inf * 0, Inf - Inf) is missing
   (removed under na.rm, which always scans for them). cosine is
   dot / (norm(x) * norm(y)), NaN for a zero vector or any Inf. */
static SEXP simd_dot_impl(SEXP x, SEXP y, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r[3];
  rsimd_opts o;
  rsimd_bin b;
  int which = rsimd_arg_int1(op, "op"), types, k;

  rsimd_in_init(&b.x, x, "x");
  rsimd_in_init(&b.y, y, "y");
  if (!is_numeric(b.x.type)) bad_type(b.x.type);
  if (!is_numeric(b.y.type)) bad_type(b.y.type);
  if (b.x.n != b.y.n) {
    Rf_error("lengths of 'x' (%lld) and 'y' (%lld) must be equal", (long long) b.x.n,
             (long long) b.y.n);
  }
  b.n = b.x.n;
  b.x_scalar = b.y_scalar = 0;
  rsimd_opts_init(&o, which == 2 ? R_NilValue : na_rm, na_check,
                  b.x.no_na_hint && b.y.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  /* The ops are symmetric, so a double operand goes first. */
  if (b.x.type != RSIMD_F64 && b.y.type == RSIMD_F64) {
    rsimd_in t = b.x;
    b.x = b.y;
    b.y = t;
  }
  types = b.x.type != RSIMD_F64   ? RSIMD_PAIR_I32_I32
          : b.y.type != RSIMD_F64 ? RSIMD_PAIR_F64_I32
                                  : RSIMD_PAIR_F64_F64;
  for (k = 0; k < 3; k++) rsimd_reduce_result_init(&r[k], RSIMD_RED_SUM);

#define RSIMD_PAIR_LOOP_(TX, TY)                                                         \
  do {                                                                                   \
    if (which == 0) {                                                                    \
      RSIMD_FOREACH_CHUNK2T(&b, TX, TY, px, py, len, off,                                \
                            { rsimd_active->dot_f64(px, py, len, types, r, &o); });      \
    } else if (which == 1) {                                                             \
      RSIMD_FOREACH_CHUNK2T(&b, TX, TY, px, py, len, off,                                \
                            { rsimd_active->dist_f64(px, py, len, types, r, &o); });     \
    } else {                                                                             \
      RSIMD_FOREACH_CHUNK2T(&b, TX, TY, px, py, len, off,                                \
                            { rsimd_active->cosine_f64(px, py, len, types, r, &o); });   \
    }                                                                                    \
  } while (0)
  switch (types) {
  case RSIMD_PAIR_F64_F64: RSIMD_PAIR_LOOP_(double, double); break;
  case RSIMD_PAIR_F64_I32: RSIMD_PAIR_LOOP_(double, int); break;
  default: RSIMD_PAIR_LOOP_(int, int); break;
  }
#undef RSIMD_PAIR_LOOP_

  /* Missing values are recorded as for doubles (an int32 NA reads as NA). */
  if (which == 0) return rsimd_reduce_finish(RSIMD_RED_DOT, RSIMD_F64, b.n, &r[0], &o);
  rsimd_reduce_settle(&r[0], o.precision);
  if (which == 1) {
    r[0].f64 = sqrt(r[0].f64);
    return rsimd_reduce_finish(RSIMD_RED_DIST, RSIMD_F64, b.n, &r[0], &o);
  }
  for (k = 1; k < 3; k++) {
    rsimd_reduce_settle(&r[k], o.precision);
    if (r[k].saw_na) r[0].saw_na = 1;
    if (r[k].saw_nan) r[0].saw_nan = 1;
  }
  r[0].f64 = r[0].f64 / (sqrt(r[1].f64) * sqrt(r[2].f64));
  return rsimd_reduce_finish(RSIMD_RED_COSINE, RSIMD_F64, b.n, &r[0], &o);
}

SEXP C_simd_dot(SEXP x, SEXP y, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_dot_impl(x, y, op, na_rm, na_check, precision));
}

/* var (sd = FALSE) or sd (sd = TRUE) of x, as base R: the mean (the sum
   in the precision mode, divided by the count, then refined for doubles by
   the mean of the deviations from it) and then the sum of squared
   deviations from it, divided by n - 1. Fewer than two elements (after
   na.rm), or any missing value without na.rm, give NA. */
static SEXP simd_var_impl(SEXP x, SEXP sd, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r, d;
  rsimd_opts o;
  rsimd_in in;
  int is_sd = rsimd_arg_lgl1(sd, "sd");
  int op = is_sd ? RSIMD_RED_SD : RSIMD_RED_VAR;
  double m, v;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->sum_f64(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->sum_i32(px, len, &r, &o);
      RSIMD_STOP_AT_NA(r);
    });
    break;
  default: bad_type(in.type);
  }
  if (r.count < 2 || (!o.na_rm && (r.saw_na || r.saw_nan))) {
    return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
  }
  /* Refined in every mode, as base R does: an error d in the mean adds
     n d^2 / (n - 1) to the result, as large as the variance itself when the
     spread is small next to the mean. */
  m = mean_value(&in, &o, &r, 1);
  if (!R_FINITE(m)) {
    /* Only for an infinite element (mean_value() rescales a sum of finite
       elements that overflowed), whose squared deviation is NaN; settle
       that here, because a vector tier's tail lanes (filled with m) would
       add NaN. */
    int inf = 0;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      if (rsimd_active->pred_f64(RSIMD_PRED_INFINITE, px, len, RSIMD_PRED_ANY, NULL)) {
        inf = 1;
        break;
      }
    });
    return Rf_ScalarReal(inf ? R_NaN : R_PosInf);
  }
  rsimd_reduce_result_init(&d, RSIMD_RED_SUM);
  if (in.type == RSIMD_F64) {
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        { rsimd_active->var_pass2_f64(px, len, m, &d, &o); });
  } else {
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off,
                        { rsimd_active->var_pass2_i32(px, len, m, &d, &o); });
  }
  v = rsimd_reduce_value(&d, o.precision) / (double) (r.count - 1);
  return Rf_ScalarReal(is_sd ? sqrt(v) : v);
}

SEXP C_simd_var(SEXP x, SEXP sd, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_var_impl(x, sd, na_rm, na_check, precision));
}

/* ---- Scans ---------------------------------------------------------------- */

/* Fills out[i], out[i + 1] ... of a double scan whose input x has its
   first NaN at i: NaN up to the first NA, NA from there on, as base R. */
static void fill_after_nan(const rsimd_in *in, R_xlen_t i, double *out) {
  double buf[RSIMD_CHUNK];
  R_xlen_t len, j, tick = 0;
  int na = 0;
  while (i < in->n) {
    const double *px = (const double *) rsimd_in_region(in, i, &len, buf);
    if (len > rsimd_stride) len = rsimd_stride;
    for (j = 0; j < len; j++) {
      if (!na && rsimd_is_na_f64(px[j])) na = 1;
      out[i + j] = na ? NA_REAL : R_NaN;
    }
    i += len;
    tick += len;
    if (tick >= rsimd_stride) {
      tick = 0;
      rsimd_check_interrupt();
    }
  }
}

/* cumsum (op 0), cumprod (op 1), cummin (op 2) or cummax (op 3) of x, with
   base R's result types: double for double x and for cumprod, integer for
   the others on integer and logical x, integer64 for cumsum, cummin and
   cummax of integer64 x. Everything from the first missing value on is
   missing; an integer cumsum that leaves the int32 range is NA from there
   on, with base R's warning, and an integer64 one that leaves int64 with
   bit64's. cumsum and cumprod of complex x are complex, computed and
   NA-fixed as base R does. */
static SEXP simd_scan_impl(SEXP x, SEXP op, SEXP precision) {
  rsimd_scan_state s;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op");
  R_xlen_t stop = -1, i;
  SEXP out;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, Rf_ScalarLogical(TRUE), in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  s.f64 = which == 1 ? 1.0 : which == 2 ? R_PosInf : which == 3 ? R_NegInf : 0.0;
  s.comp = 0.0;
  s.i64 = which == 2 ? INT64_MAX : which == 3 ? -INT64_MAX : 0;
  s.overflow = 0;
  if (in.type == RSIMD_I64 && which != 1) {
    int64_t *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_I64, in.n));
    po = (int64_t *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      R_xlen_t k_ = rsimd_active->scan_i64(which, (const int64_t *) px, len, po + off, &s);
      if (k_ >= 0) {
        stop = off + k_;
        break;
      }
    });
    for (i = stop < 0 ? in.n : stop; i < in.n; i++) po[i] = RSIMD_NA_I64;
    rsimd_set_i64_class(out);
    if (s.overflow) rsimd_warn_i64_overflow();
    UNPROTECT(1);
    return out;
  }
  if (in.type == RSIMD_C128 && which <= 1) return rsimd_c128_scan(&in, which);
  if (!is_numeric(in.type)) bad_type(in.type);
  out = PROTECT(rsimd_alloc_like(in.type == RSIMD_F64 || which == 1 ? RSIMD_F64 : RSIMD_I32, in.n));

#define RSIMD_SCAN_LOOP_(T, call)                                                        \
  RSIMD_FOREACH_CHUNK(&in, T, px, len, off, {                                            \
    R_xlen_t k_ = (call);                                                                \
    if (k_ >= 0) {                                                                       \
      stop = off + k_;                                                                   \
      break;                                                                             \
    }                                                                                    \
  })
  if (in.type == RSIMD_F64) {
    double *po = (double *) rsimd_out_ptr(out);
    switch (which) {
    case 0: RSIMD_SCAN_LOOP_(double, rsimd_active->cumsum_f64(px, len, po + off, &s, &o)); break;
    case 1: RSIMD_SCAN_LOOP_(double, rsimd_active->cumprod_f64(px, len, po + off, &s)); break;
    default:
      RSIMD_SCAN_LOOP_(double, rsimd_active->cumminmax_f64(px, len, which == 3, po + off, &s));
      break;
    }
    if (stop >= 0) fill_after_nan(&in, stop, po);
  } else if (which == 1) {
    double *po = (double *) rsimd_out_ptr(out);
    RSIMD_SCAN_LOOP_(int, rsimd_active->cumprod_i32(px, len, po + off, &s));
    for (i = stop < 0 ? in.n : stop; i < in.n; i++) po[i] = NA_REAL;
  } else {
    int *po = (int *) rsimd_out_ptr(out);
    if (which == 0) {
      RSIMD_SCAN_LOOP_(int, rsimd_active->cumsum_i32(px, len, po + off, &s));
    } else {
      RSIMD_SCAN_LOOP_(int, rsimd_active->cumminmax_i32(px, len, which == 3, po + off, &s));
    }
    for (i = stop < 0 ? in.n : stop; i < in.n; i++) po[i] = NA_INTEGER;
    if (s.overflow) rsimd_warn("integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'");
  }
#undef RSIMD_SCAN_LOOP_
  UNPROTECT(1);
  return out;
}

SEXP C_simd_scan(SEXP x, SEXP op, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(
    rsimd_sv_result(simd_scan_impl(x, op, precision), rsimd_arg_int1(op, "op") >= 2));
}
