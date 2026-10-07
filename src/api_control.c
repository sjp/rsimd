#include <R.h>
#include <Rinternals.h>
#include "rsimd.h"
#include "cpu_features.h"
#include "rsimd_config.h"
#include "dispatch.h"
#include "rvec.h"

SEXP C_simd_cpu_features(void) {
  const rsimd_cpu_features *f = rsimd_cpu();
  const rsimd_cpu_features *raw = rsimd_cpu_unmasked();
  const char *names[] = {"arch",   "os",        "method", "features", "sve_vector_length_bits",
                         "masked", "apple_core", ""};
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
  SET_VECTOR_ELT(out, 6, Rf_ScalarLogical(f->apple_core));

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

static void set_config_attr(SEXP out, const char *name, const char *value) {
  SEXP v = PROTECT(Rf_mkString(value));
  Rf_setAttrib(out, Rf_install(name), v);
  UNPROTECT(1);
}

/* Tiers compiled into the package, in simd_tiers() order. Each tier's
   table names its tier through the tier_name slot when the CPU can run it
   (tier code is never entered otherwise; referencing the table still proves
   the object is linked). Attributes: "configured", configure's list,
   "sleef", the tiers built with SLEEF elementary functions, "disabled",
   the tiers excluded by RSIMD_DISABLE_TIERS at build time, and
   "test_hole", the slot emptied by RSIMD_TEST_HOLE ("" if none). */
SEXP C_simd_compiled_tiers(void) {
  SEXP out;
  int i, n = 0, k = 0;
  for (i = 0; i < RSIMD_TIER_COUNT; i++) n += rsimd_tier_compiled((rsimd_tier) i);
  out = PROTECT(Rf_allocVector(STRSXP, n));
  for (i = 0; i < RSIMD_TIER_COUNT; i++) {
    const struct rsimd_kernels *own = rsimd_tier_table((rsimd_tier) i);
    if (own == NULL) continue;
    SET_STRING_ELT(out, k++, Rf_mkChar(rsimd_cpu_supports((rsimd_tier) i) && own->tier_name != NULL
                                         ? own->tier_name()
                                         : rsimd_tier_names[i]));
  }
  /* Each value is protected while its symbol is installed: argument
     evaluation order is unspecified, and Rf_install() can allocate. */
  set_config_attr(out, "configured", RSIMD_CONFIG_TIERS);
  set_config_attr(out, "sleef", RSIMD_CONFIG_SLEEF);
  set_config_attr(out, "disabled", RSIMD_CONFIG_DISABLED);
  set_config_attr(out, "test_hole", RSIMD_CONFIG_TEST_HOLE);
  UNPROTECT(1);
  return out;
}

/* Available tiers in preference order, ending with "none". */
SEXP C_simd_available(void) {
  SEXP out;
  int p, n = 0, k = 0;
  for (p = 0; p < rsimd_tier_preference_count; p++) {
    n += rsimd_tier_available(rsimd_tier_preference[p]);
  }
  out = PROTECT(Rf_allocVector(STRSXP, n));
  for (p = 0; p < rsimd_tier_preference_count; p++) {
    rsimd_tier t = rsimd_tier_preference[p];
    if (rsimd_tier_available(t)) SET_STRING_ELT(out, k++, Rf_mkChar(rsimd_tier_names[t]));
  }
  UNPROTECT(1);
  return out;
}

/* Selects "auto" or a tier; returns 0, -1 (unknown) or -2 (unavailable). */
SEXP C_simd_select(SEXP name) {
  if (!Rf_isString(name) || XLENGTH(name) != 1 || STRING_ELT(name, 0) == NA_STRING) {
    Rf_error("'impl' must be a single string");
  }
  rsimd_impl_selected();
  return Rf_ScalarInteger(rsimd_select(CHAR(STRING_ELT(name, 0))));
}

/* The active tier. */
SEXP C_simd_current(void) {
  return Rf_mkString(rsimd_tier_names[rsimd_active_tier]);
}

/* For an available tier (NULL: the active one), the tier whose kernel each
   slot runs after fill-down, named by slot. */
SEXP C_simd_kernel_tiers(SEXP tier) {
  rsimd_tier t = rsimd_active_tier;
  SEXP out, nms;
  int s;
  if (!Rf_isNull(tier)) {
    if (!Rf_isString(tier) || XLENGTH(tier) != 1 || STRING_ELT(tier, 0) == NA_STRING) {
      Rf_error("'tier' must be NULL or a single string");
    }
    t = rsimd_tier_from_name(CHAR(STRING_ELT(tier, 0)));
    if (!rsimd_tier_available(t)) {
      Rf_error("tier '%s' is not available", CHAR(STRING_ELT(tier, 0)));
    }
  }
  out = PROTECT(Rf_allocVector(STRSXP, RSIMD_SLOT_COUNT));
  nms = PROTECT(Rf_allocVector(STRSXP, RSIMD_SLOT_COUNT));
  for (s = 0; s < RSIMD_SLOT_COUNT; s++) {
    SET_STRING_ELT(out, s, Rf_mkChar(rsimd_tier_names[rsimd_slot_source(t, s)]));
    SET_STRING_ELT(nms, s, Rf_mkChar(rsimd_slot_names[s]));
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

/* Calls the internal slots through the active table: what the tier_name
   and fill_probe kernels it holds actually return. */
SEXP C_simd_probe_slots(void) {
  const char *names[] = {"tier_name", "fill_probe", ""};
  SEXP out;
  rsimd_entry(); /* a call that errored may have left a pinned table active */
  out = PROTECT(Rf_mkNamed(STRSXP, names));
  SET_STRING_ELT(out, 0, Rf_mkChar(rsimd_active->tier_name()));
  SET_STRING_ELT(out, 1, Rf_mkChar(rsimd_active->fill_probe()));
  UNPROTECT(1);
  return out;
}
