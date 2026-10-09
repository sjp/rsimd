/* .Call entry points of softmax and log-softmax (sigmoid is a unary
   elementary function, C_simd_math1's "sigmoid"). Both copy x into the
   double result (converting integer and logical input) while taking its
   maximum m, then run the passes of the softmax_f64 kernel over the result
   in place:
     softmax      out = exp(x - m), s = sum(out), out = out / s
     log-softmax  s = sum(exp(x - m)), out = (x - m) - log(s)
   with s summed in the precision mode, so no memory beyond the result is
   used. Missing and infinite values give what base R's formulas
   exp(x - max(x)) / sum(exp(x - max(x))) and
   x - max(x) - log(sum(exp(x - max(x)))) give, without their order
   dependence: NA everywhere if x has an NA, else NaN everywhere if it has
   a NaN or m is infinite (Inf - Inf, or -Inf - -Inf when every element is
   -Inf). Neither warns. Results are bare double vectors. */

#include <math.h>
#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

/* Runs softmax_f64 op `op` over the n elements of out in place, chunk by
   chunk with interrupt checks; chunks start at multiples of rsimd_stride,
   as pairwise summation needs. */
static void run_pass(int op, double *out, R_xlen_t n, double m, double c, rsimd_reduce_result *r,
                     const rsimd_opts *o) {
  R_xlen_t off, len;
  for (off = 0; off < n; off += len) {
    len = n - off < rsimd_stride ? n - off : rsimd_stride;
    rsimd_active->softmax_f64(op, out + off, len, m, c, out + off, r, o);
    if (off + len < n) rsimd_check_interrupt();
  }
}

static void fill(double *out, R_xlen_t n, double v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) out[i] = v;
}

static SEXP softmax(SEXP x, SEXP precision, int log_out) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  SEXP out;
  double *po, m, s;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init_fixed(&o, 0, 1, 0);
  o.precision = rsimd_arg_precision(precision);
  o.extrema = RSIMD_EXT_MAX;
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
  po = (double *) rsimd_out_ptr(out);

  /* The copy, and the maximum of the elements that are not missing, with
     saw_na and saw_nan; nothing after an NA matters. */
  rsimd_reduce_result_init(&r, RSIMD_RED_MIN);
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      memcpy(po + off, px, (size_t) len * sizeof(double));
      rsimd_active->minmax_f64(po + off, len, 0, &r, &o);
      if (r.saw_na) break;
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->convert(RSIMD_CVT_I32_F64, RSIMD_CVT_CHECKED, px, len, po + off);
      rsimd_active->minmax_f64(po + off, len, 0, &r, &o);
      if (r.saw_na) break;
    });
    break;
  default: rsimd_type_error(in.type, "x", "non-numeric argument to mathematical function");
  }
  m = r.f64_hi;
  if (in.n == 0) {
    /* numeric(0) */
  } else if (r.saw_na) {
    fill(po, in.n, NA_REAL);
  } else if (r.saw_nan || !R_FINITE(m)) {
    fill(po, in.n, R_NaN);
  } else {
    rsimd_reduce_result_init(&r, RSIMD_RED_SUM);
    run_pass(log_out ? RSIMD_SOFTMAX_SUM : RSIMD_SOFTMAX_EXP_SUM, po, in.n, m, 0.0, &r, &o);
    /* At least 1, the term of the maximum. */
    s = rsimd_reduce_value(&r, o.precision);
    if (log_out) run_pass(RSIMD_SOFTMAX_LOG, po, in.n, m, log(s), &r, &o);
    else run_pass(RSIMD_SOFTMAX_DIV, po, in.n, m, s, &r, &o);
  }
  UNPROTECT(1);
  return out;
}

static SEXP simd_softmax_impl(SEXP x, SEXP precision) {
  return softmax(x, precision, 0);
}

SEXP C_simd_softmax(SEXP x, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_softmax_impl(x, precision), 0));
}

static SEXP simd_log_softmax_impl(SEXP x, SEXP precision) {
  return softmax(x, precision, 1);
}

SEXP C_simd_log_softmax(SEXP x, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_log_softmax_impl(x, precision), 0));
}
