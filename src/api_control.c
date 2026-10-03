#include <R.h>
#include <Rinternals.h>
#include "rsimd.h"
#include "cpu_features.h"

SEXP C_simd_cpu_features(void) {
  const rsimd_cpu_features *f = rsimd_cpu();
  const rsimd_cpu_features *raw = rsimd_cpu_unmasked();
  const char *names[] = {"arch", "os", "method", "features", "sve_vector_length_bits", "masked", ""};
  SEXP out, feat, feat_names, masked;
  int i, n_masked = 0, k = 0;

  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, Rf_mkString(f->arch));
  SET_VECTOR_ELT(out, 1, Rf_mkString(f->os));
  SET_VECTOR_ELT(out, 2, Rf_mkString(f->detection_method));

  feat = PROTECT(Rf_allocVector(LGLSXP, RSIMD_CPU_FEATURE_COUNT));
  feat_names = PROTECT(Rf_allocVector(STRSXP, RSIMD_CPU_FEATURE_COUNT));
  for (i = 0; i < RSIMD_CPU_FEATURE_COUNT; i++) {
    int has = rsimd_cpu_has(f, (rsimd_cpu_feature) i);
    LOGICAL(feat)[i] = has;
    SET_STRING_ELT(feat_names, i, Rf_mkChar(rsimd_cpu_feature_names[i]));
    if (!has && rsimd_cpu_has(raw, (rsimd_cpu_feature) i)) n_masked++;
  }
  Rf_setAttrib(feat, R_NamesSymbol, feat_names);
  SET_VECTOR_ELT(out, 3, feat);

  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(f->sve_vector_length_bits));

  masked = Rf_allocVector(STRSXP, n_masked);
  SET_VECTOR_ELT(out, 5, masked);
  for (i = 0; i < RSIMD_CPU_FEATURE_COUNT; i++) {
    if (!rsimd_cpu_has(f, (rsimd_cpu_feature) i) && rsimd_cpu_has(raw, (rsimd_cpu_feature) i)) {
      SET_STRING_ELT(masked, k++, Rf_mkChar(rsimd_cpu_feature_names[i]));
    }
  }

  UNPROTECT(3);
  return out;
}

SEXP C_simd_cpu_tiers(void) {
  SEXP out = PROTECT(Rf_allocVector(LGLSXP, RSIMD_TIER_COUNT));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, RSIMD_TIER_COUNT));
  int i;
  for (i = 0; i < RSIMD_TIER_COUNT; i++) {
    LOGICAL(out)[i] = rsimd_cpu_supports((rsimd_tier) i);
    SET_STRING_ELT(nms, i, Rf_mkChar(rsimd_tier_names[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}
