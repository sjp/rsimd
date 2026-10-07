#ifndef RSIMD_CPU_FEATURES_H
#define RSIMD_CPU_FEATURES_H

#include "tiers.h"

/* Runtime CPU and OS feature detection. Plain C with no R dependencies.

   Every feature bit, in the fixed order used by the struct, by
   simd_cpu_features() and by RSIMD_CPU_FEATURES_MASK. os_avx and os_avx512
   mean the OS has enabled saving of the AVX and AVX-512 register state. */
#define RSIMD_CPU_FEATURE_LIST(X) \
  X(sse2) X(sse3) X(ssse3) X(sse4_1) X(sse4_2) X(avx) X(avx2) X(fma) \
  X(avx512f) X(avx512bw) X(avx512dq) X(avx512vl) X(os_avx) X(os_avx512) \
  X(neon) X(fp16) X(dotprod) X(i8mm) X(bf16) X(sve) X(sve2)

typedef enum {
#define RSIMD_CPU_X(name) RSIMD_CPU_##name,
  RSIMD_CPU_FEATURE_LIST(RSIMD_CPU_X)
#undef RSIMD_CPU_X
  RSIMD_CPU_FEATURE_COUNT
} rsimd_cpu_feature;

typedef struct rsimd_cpu_features {
#define RSIMD_CPU_X(name) unsigned name : 1;
  RSIMD_CPU_FEATURE_LIST(RSIMD_CPU_X)
#undef RSIMD_CPU_X
  int sve_vector_length_bits;   /* 0 if unknown or no SVE */
  int apple_core;               /* an Apple arm64 core (implementer 0x61) */
  const char *arch;             /* "x86_64", "i686", "aarch64", "armv7", "other" */
  const char *os;               /* "linux", "darwin", "windows", "freebsd", "other" */
  const char *detection_method; /* "cpuid", "getauxval", "sysctl", "win32", "elf_aux_info", "none" */
} rsimd_cpu_features;

/* Feature names in list order. */
extern const char *const rsimd_cpu_feature_names[RSIMD_CPU_FEATURE_COUNT];

/* Detects features once and caches them. rsimd_cpu_init() is called when
   the shared library is loaded; rsimd_cpu() initialises on first use if
   needed. After initialisation the structs never change, so they may be
   read from any context without synchronisation.

   Detected bits are made consistent (a feature is reported only when the
   features it builds on are, e.g. avx2 needs avx and os_avx), then the
   comma-separated feature names in the environment variable
   RSIMD_CPU_FEATURES_MASK are forced off, together with every feature that
   builds on them. The mask is for testing fallback paths; it never turns a
   feature on. */
void rsimd_cpu_init(void);
const rsimd_cpu_features *rsimd_cpu(void);

/* The features as detected, before RSIMD_CPU_FEATURES_MASK was applied. */
const rsimd_cpu_features *rsimd_cpu_unmasked(void);

int rsimd_cpu_has(const rsimd_cpu_features *f, rsimd_cpu_feature which);

/* 1 if the CPU and OS can run tier t (whether or not it was compiled in).
   "none" is always supported; reserved tiers never are. */
int rsimd_cpu_supports(rsimd_tier t);

#endif /* RSIMD_CPU_FEATURES_H */
