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
  if (in.type != RSIMD_F64) {
    long double s = (long double) r.i64;
    if (r.overflow) s += (long double) r.f64;
    return Rf_ScalarReal((double) (s / (long double) r.count));
  }
  m = rsimd_reduce_value(&r, o.precision) / (double) r.count;
  if (o.precision != RSIMD_PREC_FAST && R_FINITE(m)) {
    rsimd_reduce_result d;
    rsimd_reduce_result_init(&d, RSIMD_RED_SUM);
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                        { rsimd_active->sum_dev_f64(px, len, m, &d, &o); });
    m += rsimd_reduce_value(&d, o.precision) / (double) r.count;
  }
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
