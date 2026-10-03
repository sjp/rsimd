#include <string.h>
#include "tiers.h"

const char *const rsimd_tier_names[RSIMD_TIER_COUNT] = {
  "none", "sse2", "avx2", "avx512", "neon", "sve", "sve2", "rvv", "wasm128"
};

rsimd_tier rsimd_tier_from_name(const char *s) {
  int i;
  if (s == NULL) return RSIMD_TIER_COUNT;
  for (i = 0; i < RSIMD_TIER_COUNT; i++) {
    if (strcmp(s, rsimd_tier_names[i]) == 0) return (rsimd_tier) i;
  }
  return RSIMD_TIER_COUNT;
}
