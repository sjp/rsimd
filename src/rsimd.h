#ifndef RSIMD_H
#define RSIMD_H

#include <Rinternals.h>

/* .Call entry points, registered in init.c. */

/* control */
SEXP C_simd_version(void);
SEXP C_simd_cpu_features(void);
SEXP C_simd_cpu_tiers(void);
SEXP C_simd_compiled_tiers(void);
SEXP C_simd_available(void);
SEXP C_simd_select(SEXP name);
SEXP C_simd_current(void);
SEXP C_simd_kernel_tiers(SEXP tier);
SEXP C_simd_probe_slots(void);

/* access layer */
SEXP C_simd_promote(SEXP x, SEXP y);
SEXP C_simd_debug_regions(SEXP x);
SEXP C_simd_debug_copy(SEXP x, SEXP no_na);
SEXP C_simd_debug_bin(SEXP x, SEXP y);
SEXP C_simd_debug_finish(SEXP op, SEXP type, SEXP n, SEXP fields, SEXP precision, SEXP na_rm);
SEXP C_simd_debug_opts(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision);

/* self-test kernels (NA, precision and overflow rules) */
SEXP C_simd_debug_fold(SEXP x, SEXP y, SEXP term, SEXP na_rm, SEXP na_check, SEXP precision);
SEXP C_simd_debug_lgl(SEXP x, SEXP op, SEXP na_rm, SEXP na_check);
SEXP C_simd_debug_arith(SEXP x, SEXP y, SEXP op, SEXP na_check);

#endif /* RSIMD_H */
