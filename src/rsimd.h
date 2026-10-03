#ifndef RSIMD_H
#define RSIMD_H

#include <Rinternals.h>

/* .Call entry points, registered in init.c. */

/* control */
SEXP C_simd_version(void);
SEXP C_simd_cpu_features(void);
SEXP C_simd_cpu_tiers(void);
SEXP C_simd_compiled_tiers(void);

#endif /* RSIMD_H */
