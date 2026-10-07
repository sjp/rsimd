#ifndef RSIMD_TIERS_H
#define RSIMD_TIERS_H

/* Target architecture, from the compiler's predefined macros so it always
   matches the code being compiled. RSIMD_ARCH_X86 and RSIMD_ARCH_ARM cover
   both the 64-bit and the 32-bit variants. */
#if defined(__x86_64__) || defined(_M_X64)
#define RSIMD_ARCH_X86 1
#define RSIMD_ARCH_X86_64 1
#elif defined(__i386__) || defined(_M_IX86)
#define RSIMD_ARCH_X86 1
#define RSIMD_ARCH_X86_32 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#define RSIMD_ARCH_ARM 1
#define RSIMD_ARCH_ARM64 1
#elif defined(__arm__) || defined(_M_ARM)
#define RSIMD_ARCH_ARM 1
#define RSIMD_ARCH_ARM32 1
#endif

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

/* 1 on an Apple arm64 core (CPU implementer 0x61), set when the CPU is
   detected; a kernel may take a path that is faster on those cores. */
extern int rsimd_apple_core;

#endif /* RSIMD_TIERS_H */
