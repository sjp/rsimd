#include <R.h>
#include <Rinternals.h>
#include "rsimd.h"

/* Version of the package's native ABI. Bumped whenever the compiled code
   changes in a way the R code must know about. */
#define RSIMD_NATIVE_VERSION "0"

SEXP C_simd_version(void) {
  return Rf_mkString(RSIMD_NATIVE_VERSION);
}
