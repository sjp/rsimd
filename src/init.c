#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rsimd.h" /* declares every C_simd_* entry point */
#include "cpu_features.h"

static const R_CallMethodDef CallEntries[] = {
  {"C_simd_version", (DL_FUNC) &C_simd_version, 0},
  {"C_simd_cpu_features", (DL_FUNC) &C_simd_cpu_features, 0},
  {"C_simd_cpu_tiers", (DL_FUNC) &C_simd_cpu_tiers, 0},
  {"C_simd_compiled_tiers", (DL_FUNC) &C_simd_compiled_tiers, 0},
  {NULL, NULL, 0}
};

void R_init_rsimd(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  rsimd_cpu_init();
}
