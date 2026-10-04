#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rsimd.h" /* declares every C_simd_* entry point */
#include "cpu_features.h"
#include "dispatch.h"
#include "rvec.h"

/* The cast through void (*)(void), the generic function type, avoids
   -Wcast-function-type: R calls each routine with its declared arguments. */
#define CALLDEF(name, n) {#name, (DL_FUNC) (void (*)(void)) &name, n}

static const R_CallMethodDef CallEntries[] = {
  CALLDEF(C_simd_version, 0),
  CALLDEF(C_simd_cpu_features, 0),
  CALLDEF(C_simd_cpu_tiers, 0),
  CALLDEF(C_simd_compiled_tiers, 0),
  CALLDEF(C_simd_available, 0),
  CALLDEF(C_simd_select, 1),
  CALLDEF(C_simd_current, 0),
  CALLDEF(C_simd_kernel_tiers, 1),
  CALLDEF(C_simd_probe_slots, 0),
  CALLDEF(C_simd_promote, 2),
  CALLDEF(C_simd_debug_regions, 1),
  CALLDEF(C_simd_debug_copy, 2),
  CALLDEF(C_simd_debug_bin, 2),
  CALLDEF(C_simd_debug_finish, 6),
  CALLDEF(C_simd_debug_opts, 4),
  CALLDEF(C_simd_sum, 4),
  CALLDEF(C_simd_prod, 3),
  CALLDEF(C_simd_mean, 4),
  CALLDEF(C_simd_minmax, 4),
  CALLDEF(C_simd_which, 2),
  CALLDEF(C_simd_anyall, 3),
  CALLDEF(C_simd_na, 2),
  CALLDEF(C_simd_sum_sq, 5),
  CALLDEF(C_simd_dot, 6),
  CALLDEF(C_simd_var, 5),
  CALLDEF(C_simd_scan, 3),
  CALLDEF(C_simd_ew1, 2),
  CALLDEF(C_simd_ew2, 4),
  CALLDEF(C_simd_ew3, 5),
  CALLDEF(C_simd_round_digits, 2),
  CALLDEF(C_simd_pred, 3),
  CALLDEF(C_simd_cmp, 3),
  CALLDEF(C_simd_logic, 3),
  CALLDEF(C_simd_bit, 5),
  CALLDEF(C_simd_popcount_total, 3),
  CALLDEF(C_simd_convert, 3),
  CALLDEF(C_simd_math1, 2),
  CALLDEF(C_simd_log, 2),
  CALLDEF(C_simd_math2, 3),
  CALLDEF(C_simd_sincos, 1),
  CALLDEF(C_simd_ulp_dist, 2),
  CALLDEF(C_simd_softmax, 2),
  CALLDEF(C_simd_log_softmax, 2),
  CALLDEF(C_simd_debug_fold, 6),
  CALLDEF(C_simd_debug_lgl, 4),
  CALLDEF(C_simd_debug_arith, 4),
  {NULL, NULL, 0}
};

void R_init_rsimd(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  rsimd_cpu_init();
  rsimd_dispatch_init();
  rsimd_rvec_init();
}
