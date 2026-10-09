#include <R.h>
#include <Rinternals.h>
#include "rsimd.h"
#include "rsimd_config.h"

/* Clang's __VERSION__ names the compiler ("Clang 18.1.8", "Apple LLVM
   ..."); GCC's is the version number alone. */
#if defined(__clang__)
#define RSIMD_COMPILER __VERSION__
#elif defined(__GNUC__)
#define RSIMD_COMPILER "GCC " __VERSION__
#else
#define RSIMD_COMPILER ""
#endif

static SEXP version_string(const char *v) {
  return v[0] == '\0' ? Rf_ScalarString(NA_STRING) : Rf_mkString(v);
}

/* The bundled SIMDe and SLEEF versions (SLEEF NA when no tier uses it) and
   the compiler that built the package (NA if unknown), as a named list. */
SEXP C_simd_build_info(void) {
  const char *names[] = {"simde", "sleef", "compiler", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, version_string(RSIMD_SIMDE_VERSION));
  SET_VECTOR_ELT(out, 1, version_string(RSIMD_SLEEF_VERSION));
  SET_VECTOR_ELT(out, 2, version_string(RSIMD_COMPILER));
  UNPROTECT(1);
  return out;
}
