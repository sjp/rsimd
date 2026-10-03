#ifndef RSIMD_TIERS_H
#define RSIMD_TIERS_H

/* Implementation tiers. The enum order is fixed and matches simd_tiers();
   it is not the preference order used by "auto". This header has no R
   dependencies so plain C code (CPU detection) can use it. */
typedef enum {
  RSIMD_TIER_NONE = 0,
  RSIMD_TIER_SSE2,
  RSIMD_TIER_AVX2,
  RSIMD_TIER_AVX512,
  RSIMD_TIER_NEON,
  RSIMD_TIER_SVE,
  RSIMD_TIER_SVE2,
  RSIMD_TIER_RVV,     /* reserved: never compiled */
  RSIMD_TIER_WASM128, /* reserved: never compiled */
  RSIMD_TIER_COUNT
} rsimd_tier;

/* "none", "sse2", ... in enum order. */
extern const char *const rsimd_tier_names[RSIMD_TIER_COUNT];

/* Returns RSIMD_TIER_COUNT for an unknown name. "auto" is not a tier. */
rsimd_tier rsimd_tier_from_name(const char *s);

#endif /* RSIMD_TIERS_H */
