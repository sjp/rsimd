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

#endif /* RSIMD_H */
