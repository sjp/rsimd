/* .Call entry points of the reductions: classify the input, fold it chunk
   by chunk through the active implementation, and finish with base R's
   result types (rvec.h). Empty inputs run no kernel: the identity of the
   reduction is finished directly. */

#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

static void bad_type(rsimd_etype type) {
  Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[type]);
}

/* sum(x, na.rm): double for double x; integer for integer and logical x,
   or double when the total does not fit in an integer. precision is the
   integer code RSIMD_PREC_*. */
SEXP C_simd_sum(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->sum_f64(px, len, &r, &o); });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->sum_i32(px, len, &r, &o); });
    break;
  default: bad_type(in.type);
  }
  return rsimd_reduce_finish(RSIMD_RED_SUM, in.type, in.n, &r, &o);
}

/* prod(x, na.rm): always double. The precision mode does not apply. */
SEXP C_simd_prod(SEXP x, SEXP na_rm, SEXP na_check) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, R_NilValue, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_PROD);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->prod_f64(px, len, &r, &o); });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->prod_i32(px, len, &r, &o); });
    break;
  default: bad_type(in.type);
  }
  return rsimd_reduce_finish(RSIMD_RED_PROD, in.type, in.n, &r, &o);
}

/* The mean of x, from the sum fold r (sum_f64 or sum_i32 over every chunk,
   with r->count > 0 elements left): the exact integer sum divided in long
   double, or the double sum divided by the count, refined in pairwise and
   compensated modes by the mean of the deviations from it when finite. */
static double mean_value(const rsimd_in *in, const rsimd_opts *o, const rsimd_reduce_result *r) {
  double m;
  if (in->type != RSIMD_F64) {
    long double s = (long double) r->i64;
    if (r->overflow) s += (long double) r->f64;
    return (double) (s / (long double) r->count);
  }
  m = rsimd_reduce_value(r, o->precision) / (double) r->count;
  if (o->precision != RSIMD_PREC_FAST && R_FINITE(m)) {
    rsimd_reduce_result d;
    rsimd_reduce_result_init(&d, RSIMD_RED_SUM);
    RSIMD_FOREACH_CHUNK(in, double, px, len, off,
                        { rsimd_active->sum_dev_f64(px, len, m, &d, o); });
    m += rsimd_reduce_value(&d, o->precision) / (double) r->count;
  }
  return m;
}

/* mean(x, na.rm): always double. Doubles: the sum divided by the count; in
   pairwise and compensated modes refined, as base R does, by adding the
   mean of the deviations from it (a second pass) when it is finite.
   Integers and logicals: the exact sum divided in long double. */
SEXP C_simd_mean(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  double m;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_MEAN);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->sum_f64(px, len, &r, &o); });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->sum_i32(px, len, &r, &o); });
    break;
  default: bad_type(in.type);
  }
  /* Missing values and empty input. */
  if (r.count == 0 || (!o.na_rm && (r.saw_na || r.saw_nan))) {
    return rsimd_reduce_finish(RSIMD_RED_MEAN, in.type, in.n, &r, &o);
  }
  m = mean_value(&in, &o, &r);
  return Rf_ScalarReal(m);
}

/* min (op 0), max (op 1) or range (op 2) of x: both extrema come from one
   pass. Integer and logical results are integer, except for empty input
   (after na.rm), which gives Inf and -Inf with base R's warnings. */
SEXP C_simd_minmax(SEXP x, SEXP op, SEXP na_rm, SEXP na_check) {
  rsimd_reduce_result r, rmax;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op");
  SEXP lo, hi, out;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, R_NilValue, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_MIN);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->minmax_f64(px, len, &r, &o); });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->minmax_i32(px, len, &r, &o); });
    break;
  default: bad_type(in.type);
  }
  rmax = r;
  rmax.f64 = r.f64_hi;
  rmax.i64 = r.i64_hi;
  if (which == 0) return rsimd_reduce_finish(RSIMD_RED_MIN, in.type, in.n, &r, &o);
  if (which == 1) return rsimd_reduce_finish(RSIMD_RED_MAX, in.type, in.n, &rmax, &o);
  lo = PROTECT(rsimd_reduce_finish(RSIMD_RED_MIN, in.type, in.n, &r, &o));
  hi = PROTECT(rsimd_reduce_finish(RSIMD_RED_MAX, in.type, in.n, &rmax, &o));
  out = rsimd_range_pair(lo, hi);
  UNPROTECT(2);
  return out;
}

/* which.min (max = FALSE) or which.max (max = TRUE): missing values are
   ignored. The extremum comes from the min/max kernels, then the first
   element equal to it is looked up (so 0 and -0 tie, as in base R). Raw
   vectors use a one-pass kernel. */
SEXP C_simd_which(SEXP x, SEXP max) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int dir = rsimd_arg_lgl1(max, "max");
  int op = dir ? RSIMD_RED_WHICH_MAX : RSIMD_RED_WHICH_MIN;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, R_NilValue, R_NilValue, in.no_na_hint);
  o.na_rm = 1;
  rsimd_reduce_result_init(&r, op);
  switch (in.type) {
  case RSIMD_F64: {
    double v;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->minmax_f64(px, len, &r, &o); });
    if (r.count == 0) break;
    v = dir ? r.f64_hi : r.f64;
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      R_xlen_t k = rsimd_active->find_f64(px, len, v);
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
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->minmax_i32(px, len, &r, &o); });
    if (r.count == 0) break;
    v = (int) (dir ? r.i64_hi : r.i64);
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      R_xlen_t k = rsimd_active->find_i32(px, len, v);
      if (k >= 0) {
        r.idx = off + k;
        break;
      }
    });
    break;
  }
  case RSIMD_U8:
    RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off,
                        { rsimd_active->which_u8(px, len, off, dir, &r); });
    break;
  default: bad_type(in.type);
  }
  return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
}

/* any (all = FALSE) or all (all = TRUE) with three-valued logic. Double and
   raw inputs are read as logical values without being converted, with
   base R's warning; reading stops once the answer is known. */
SEXP C_simd_anyall(SEXP x, SEXP all, SEXP na_rm) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int is_all = rsimd_arg_lgl1(all, "all");
  int op = is_all ? RSIMD_RED_ALL : RSIMD_RED_ANY;
  int stop = is_all ? RSIMD_STOP_FALSE : RSIMD_STOP_TRUE;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, R_NilValue, R_NilValue, in.no_na_hint);
  rsimd_reduce_result_init(&r, op);
#define RSIMD_DONE_ (is_all ? r.any_false : r.any_true)
  switch (in.type) {
  case RSIMD_F64:
  case RSIMD_U8:
    if (in.n > 0) {
      Rf_warning("coercing argument of type '%s' to logical",
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
  default: bad_type(in.type);
  }
#undef RSIMD_DONE_
  return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
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
  default: bad_type(in->type);
  }
}

/* any_na (mode 0, logical), count_na (mode 1, double) or which_na (mode 2,
   the 1-based indices of the missing elements: integer, or double when one
   of them exceeds INT_MAX). NaN counts as missing. Raw vectors and inputs R knows
   to be NA-free are not read. */
SEXP C_simd_na(SEXP x, SEXP mode) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int m = rsimd_arg_int1(mode, "mode");
  SEXP out;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, R_NilValue, R_NilValue, in.no_na_hint);
  rsimd_reduce_result_init(&r, m == 0 ? RSIMD_RED_ANY_NA : RSIMD_RED_COUNT_NA);
  if (in.type == RSIMD_I64) bad_type(in.type);
  if (!in.no_na_hint) na_scan(&in, m == 0 ? RSIMD_NAMODE_ANY : RSIMD_NAMODE_COUNT, NULL, &r);
  if (m == 0) return rsimd_reduce_finish(RSIMD_RED_ANY_NA, in.type, in.n, &r, &o);
  if (m == 1) return rsimd_reduce_finish(RSIMD_RED_COUNT_NA, in.type, in.n, &r, &o);
  /* which_na: count, then write the indices; doubles for a long vector,
     converted to integer if the last index fits after all. */
  {
    R_xlen_t count = (R_xlen_t) r.i64, i;
    rsimd_etype t = in.n > INT_MAX && count > 0 ? RSIMD_F64 : RSIMD_I32;
    out = PROTECT(rsimd_alloc_like(t, count));
    if (count > 0) {
      r.i64 = 0;
      na_scan(&in, t == RSIMD_F64 ? RSIMD_NAMODE_WHICH_F64 : RSIMD_NAMODE_WHICH_I32,
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
  }
  return out;
}

/* ---- Sums of squares, norms, distances, variance -------------------------- */

static int is_numeric(rsimd_etype t) {
  return t == RSIMD_F64 || t == RSIMD_I32 || t == RSIMD_LGL;
}

/* sum_sq (op 0), norm (op 1) or sum_abs (op 2) of x. sum_abs is a sum of
   |x| with simd_sum's result types (integer for integer and logical x,
   double when the total does not fit); the others are double. */
SEXP C_simd_sum_sq(SEXP x, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op");

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    if (which == 2) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                          { rsimd_active->sumabs_f64(px, len, &r, &o); });
    } else {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->sumsq_f64(px, len, &r, &o); });
    }
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    if (which == 2) {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->sumabs_i32(px, len, &r, &o); });
    } else {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->sumsq_i32(px, len, &r, &o); });
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

/* dot (op 0), dist (op 1) or cosine (op 2) of x and y, which must have the
   same length (no broadcast). Double, integer and logical operands mix
   without conversion: the kernels read int32 elements as doubles. A pair
   with a missing element is missing (removed under na.rm). cosine is
   dot / (norm(x) * norm(y)), NaN for a zero vector or any Inf. */
SEXP C_simd_dot(SEXP x, SEXP y, SEXP op, SEXP na_rm, SEXP na_check, SEXP precision) {
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
  rsimd_opts_init(&o, which == 2 ? R_NilValue : na_rm, na_check, precision,
                  b.x.no_na_hint && b.y.no_na_hint);
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

/* var (sd = FALSE) or sd (sd = TRUE) of x: two passes as base R, the mean
   (as simd_mean, in the same precision mode) and then the sum of squared
   deviations from it, divided by n - 1. Fewer than two elements (after
   na.rm), or any missing value without na.rm, give NA. */
SEXP C_simd_var(SEXP x, SEXP sd, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_reduce_result r, d;
  rsimd_opts o;
  rsimd_in in;
  int is_sd = rsimd_arg_lgl1(sd, "sd");
  int op = is_sd ? RSIMD_RED_SD : RSIMD_RED_VAR;
  double m, v;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
  rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, { rsimd_active->sum_f64(px, len, &r, &o); });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, { rsimd_active->sum_i32(px, len, &r, &o); });
    break;
  default: bad_type(in.type);
  }
  if (r.count < 2 || (!o.na_rm && (r.saw_na || r.saw_nan))) {
    return rsimd_reduce_finish(op, in.type, in.n, &r, &o);
  }
  m = mean_value(&in, &o, &r);
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
   the others on integer and logical x. Everything from the first missing
   value on is missing; an integer cumsum that leaves the int32 range is NA
   from there on, with base R's warning. */
SEXP C_simd_scan(SEXP x, SEXP op, SEXP precision) {
  rsimd_scan_state s;
  rsimd_opts o;
  rsimd_in in;
  int which = rsimd_arg_int1(op, "op");
  R_xlen_t stop = -1, i;
  SEXP out;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, R_NilValue, R_NilValue, precision, in.no_na_hint);
  if (!is_numeric(in.type)) bad_type(in.type);
  s.f64 = which == 1 ? 1.0 : which == 2 ? R_PosInf : which == 3 ? R_NegInf : 0.0;
  s.comp = 0.0;
  s.overflow = 0;
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
    if (s.overflow) Rf_warning("integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'");
  }
#undef RSIMD_SCAN_LOOP_
  UNPROTECT(1);
  return out;
}
