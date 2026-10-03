#include <R.h>
#include <Rinternals.h>
#include "rsimd.h"
#include "cpu_features.h"
#include "rsimd_config.h"

/* Defined in each tier's translation unit (kernels/tier_info.inc.c). */
const char *rsimd_tier_name_none(void);
#ifdef RSIMD_HAVE_SSE2
const char *rsimd_tier_name_sse2(void);
#endif
#ifdef RSIMD_HAVE_AVX2
const char *rsimd_tier_name_avx2(void);
#endif
#ifdef RSIMD_HAVE_AVX512
const char *rsimd_tier_name_avx512(void);
#endif
#ifdef RSIMD_HAVE_NEON
const char *rsimd_tier_name_neon(void);
#endif
#ifdef RSIMD_HAVE_SVE
const char *rsimd_tier_name_sve(void);
#endif
#ifdef RSIMD_HAVE_SVE2
const char *rsimd_tier_name_sve2(void);
#endif

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

/* Tiers compiled into the package, in simd_tiers() order. Each tier object
   names itself when the CPU can run it (tier code is never entered
   otherwise; referencing the function still proves the object is linked).
   Attributes: "configured", configure's list, and
   "disabled", the tiers excluded by RSIMD_DISABLE_TIERS at build time. */
SEXP C_simd_compiled_tiers(void) {
  static const char *(*const name_fns[RSIMD_TIER_COUNT])(void) = {
    [RSIMD_TIER_NONE] = rsimd_tier_name_none,
#ifdef RSIMD_HAVE_SSE2
    [RSIMD_TIER_SSE2] = rsimd_tier_name_sse2,
#endif
#ifdef RSIMD_HAVE_AVX2
    [RSIMD_TIER_AVX2] = rsimd_tier_name_avx2,
#endif
#ifdef RSIMD_HAVE_AVX512
    [RSIMD_TIER_AVX512] = rsimd_tier_name_avx512,
#endif
#ifdef RSIMD_HAVE_NEON
    [RSIMD_TIER_NEON] = rsimd_tier_name_neon,
#endif
#ifdef RSIMD_HAVE_SVE
    [RSIMD_TIER_SVE] = rsimd_tier_name_sve,
#endif
#ifdef RSIMD_HAVE_SVE2
    [RSIMD_TIER_SVE2] = rsimd_tier_name_sve2,
#endif
  };
  SEXP out;
  int i, n = 0, k = 0;
  for (i = 0; i < RSIMD_TIER_COUNT; i++) n += name_fns[i] != NULL;
  out = PROTECT(Rf_allocVector(STRSXP, n));
  for (i = 0; i < RSIMD_TIER_COUNT; i++) {
    if (name_fns[i] == NULL) continue;
    SET_STRING_ELT(out, k++, Rf_mkChar(rsimd_cpu_supports((rsimd_tier) i)
                                         ? name_fns[i]()
                                         : rsimd_tier_names[i]));
  }
  Rf_setAttrib(out, Rf_install("configured"), Rf_mkString(RSIMD_CONFIG_TIERS));
  Rf_setAttrib(out, Rf_install("disabled"), Rf_mkString(RSIMD_CONFIG_DISABLED));
  UNPROTECT(1);
  return out;
}
