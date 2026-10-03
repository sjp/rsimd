/* .Call entry points of the reductions: classify the input, fold it chunk
   by chunk through the active implementation, and finish with base R's
   result types (rvec.h). */

#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

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
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  return rsimd_reduce_finish(RSIMD_RED_SUM, in.type, in.n, &r, &o);
}
