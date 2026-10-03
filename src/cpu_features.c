/* Runtime CPU and OS feature detection.

   Hand-rolled on purpose: IFUNC, target_clones and __builtin_cpu_supports
   are unavailable or unreliable on some of R's toolchains (MinGW, musl,
   Apple clang on arm64). This file contains no SIMD intrinsics and is built
   with R's default compiler flags. */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "cpu_features.h"

#if defined(__x86_64__) || defined(_M_X64)
#define RSIMD_CPU_X86 1
#define RSIMD_CPU_ARCH "x86_64"
#elif defined(__i386__) || defined(_M_IX86)
#define RSIMD_CPU_X86 1
#define RSIMD_CPU_ARCH "i686"
#elif defined(__aarch64__) || defined(_M_ARM64)
#define RSIMD_CPU_ARM64 1
#define RSIMD_CPU_ARCH "aarch64"
#elif defined(__arm__) || defined(_M_ARM)
#define RSIMD_CPU_ARM32 1
#define RSIMD_CPU_ARCH "armv7"
#else
#define RSIMD_CPU_ARCH "other"
#endif

#if defined(__linux__)
#define RSIMD_CPU_OS "linux"
#elif defined(__APPLE__)
#define RSIMD_CPU_OS "darwin"
#elif defined(_WIN32)
#define RSIMD_CPU_OS "windows"
#elif defined(__FreeBSD__)
#define RSIMD_CPU_OS "freebsd"
#else
#define RSIMD_CPU_OS "other"
#endif

#if defined(RSIMD_CPU_X86)
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#if (defined(RSIMD_CPU_ARM64) || defined(RSIMD_CPU_ARM32)) && \
    (defined(__linux__) || defined(__FreeBSD__))
#include <sys/auxv.h>
#define RSIMD_CPU_AUXV 1
#endif

#if defined(RSIMD_CPU_ARM64) && defined(__linux__)
#include <sys/prctl.h>
#ifndef PR_SVE_GET_VL
#define PR_SVE_GET_VL 51
#endif
#ifndef PR_SVE_VL_LEN_MASK
#define PR_SVE_VL_LEN_MASK 0xffff
#endif
#endif

/* Defined in the SVE tier's translation unit when that tier is compiled. */
#if defined(RSIMD_HAVE_SVE)
int rsimd_sve_vl_bits_sve(void);
#endif

const char *const rsimd_cpu_feature_names[RSIMD_CPU_FEATURE_COUNT] = {
#define RSIMD_CPU_X(name) #name,
  RSIMD_CPU_FEATURE_LIST(RSIMD_CPU_X)
#undef RSIMD_CPU_X
};

int rsimd_cpu_has(const rsimd_cpu_features *f, rsimd_cpu_feature which) {
  switch (which) {
#define RSIMD_CPU_X(name) case RSIMD_CPU_##name: return (int) f->name;
    RSIMD_CPU_FEATURE_LIST(RSIMD_CPU_X)
#undef RSIMD_CPU_X
  default: return 0;
  }
}

static void cpu_clear(rsimd_cpu_features *f, rsimd_cpu_feature which) {
  switch (which) {
#define RSIMD_CPU_X(name) case RSIMD_CPU_##name: f->name = 0; break;
    RSIMD_CPU_FEATURE_LIST(RSIMD_CPU_X)
#undef RSIMD_CPU_X
  default: break;
  }
}

/* A feature is only reported when every feature it builds on is. Used both
   to tidy up what the hardware reports (e.g. AVX without OS support for
   the AVX state is unusable) and to cascade RSIMD_CPU_FEATURES_MASK. */
static const struct {
  rsimd_cpu_feature feature, needs;
} cpu_deps[] = {
  {RSIMD_CPU_sse3, RSIMD_CPU_sse2},
  {RSIMD_CPU_ssse3, RSIMD_CPU_sse3},
  {RSIMD_CPU_sse4_1, RSIMD_CPU_ssse3},
  {RSIMD_CPU_sse4_2, RSIMD_CPU_sse4_1},
  {RSIMD_CPU_avx, RSIMD_CPU_sse4_2},
  {RSIMD_CPU_avx, RSIMD_CPU_os_avx},
  {RSIMD_CPU_avx2, RSIMD_CPU_avx},
  {RSIMD_CPU_fma, RSIMD_CPU_avx},
  {RSIMD_CPU_os_avx512, RSIMD_CPU_os_avx},
  {RSIMD_CPU_avx512f, RSIMD_CPU_avx2},
  {RSIMD_CPU_avx512f, RSIMD_CPU_fma},
  {RSIMD_CPU_avx512f, RSIMD_CPU_os_avx512},
  {RSIMD_CPU_avx512bw, RSIMD_CPU_avx512f},
  {RSIMD_CPU_avx512dq, RSIMD_CPU_avx512f},
  {RSIMD_CPU_avx512vl, RSIMD_CPU_avx512f},
  {RSIMD_CPU_fp16, RSIMD_CPU_neon},
  {RSIMD_CPU_dotprod, RSIMD_CPU_neon},
  {RSIMD_CPU_i8mm, RSIMD_CPU_neon},
  {RSIMD_CPU_bf16, RSIMD_CPU_neon},
  {RSIMD_CPU_sve, RSIMD_CPU_neon},
  {RSIMD_CPU_sve2, RSIMD_CPU_sve}
};

static void cpu_make_consistent(rsimd_cpu_features *f) {
  size_t i;
  int changed = 1;
  while (changed) {
    changed = 0;
    for (i = 0; i < sizeof cpu_deps / sizeof cpu_deps[0]; i++) {
      if (rsimd_cpu_has(f, cpu_deps[i].feature) && !rsimd_cpu_has(f, cpu_deps[i].needs)) {
        cpu_clear(f, cpu_deps[i].feature);
        changed = 1;
      }
    }
  }
  if (!f->sve || f->sve_vector_length_bits < 128 || f->sve_vector_length_bits > 2048 ||
      f->sve_vector_length_bits % 128 != 0) {
    f->sve_vector_length_bits = 0;
  }
}

static int ascii_lower(int c) {
  return (c >= 'A' && c <= 'Z') ? c - 'A' + 'a' : c;
}

static int is_space(int c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

/* Clears every feature named in the comma-separated list `mask`. Names are
   matched case-insensitively after trimming whitespace; unknown names are
   ignored here (the R side warns about them). */
static void cpu_apply_mask(rsimd_cpu_features *f, const char *mask) {
  const char *p = mask;
  if (mask == NULL) return;
  while (*p != '\0') {
    const char *start, *end;
    size_t len;
    int i;
    while (*p == ',' || is_space((unsigned char) *p)) p++;
    start = p;
    while (*p != '\0' && *p != ',') p++;
    end = p;
    while (end > start && is_space((unsigned char) end[-1])) end--;
    len = (size_t) (end - start);
    if (len == 0) continue;
    for (i = 0; i < RSIMD_CPU_FEATURE_COUNT; i++) {
      const char *name = rsimd_cpu_feature_names[i];
      size_t k;
      if (strlen(name) != len) continue;
      for (k = 0; k < len; k++) {
        if (ascii_lower((unsigned char) start[k]) != name[k]) break;
      }
      if (k == len) cpu_clear(f, (rsimd_cpu_feature) i);
    }
  }
}

/* ---- x86 and x86-64: cpuid and xgetbv ---------------------------------- */

#if defined(RSIMD_CPU_X86)

static unsigned cpuid_max_leaf(void) {
#if defined(_MSC_VER)
  int r[4];
  __cpuid(r, 0);
  return (unsigned) r[0];
#else
  /* Returns 0 when the cpuid instruction itself is unavailable (old i386). */
  return __get_cpuid_max(0, NULL);
#endif
}

static void cpuid(unsigned leaf, unsigned subleaf, unsigned r[4]) {
#if defined(_MSC_VER)
  int v[4];
  __cpuidex(v, (int) leaf, (int) subleaf);
  r[0] = (unsigned) v[0]; r[1] = (unsigned) v[1];
  r[2] = (unsigned) v[2]; r[3] = (unsigned) v[3];
#else
  /* <cpuid.h> preserves ebx when it is the PIC register on i686. */
  __cpuid_count(leaf, subleaf, r[0], r[1], r[2], r[3]);
#endif
}

/* Only valid when CPUID.1:ECX.OSXSAVE is set. Encoded as bytes so it also
   assembles with assemblers that predate the mnemonic. */
static unsigned xgetbv0(void) {
#if defined(_MSC_VER)
  return (unsigned) _xgetbv(0);
#else
  unsigned eax, edx;
  __asm__ __volatile__(".byte 0x0f, 0x01, 0xd0" : "=a"(eax), "=d"(edx) : "c"(0));
  (void) edx;
  return eax;
#endif
}

#if defined(__APPLE__)
static int sysctl_flag(const char *name);
#endif

static void detect_x86(rsimd_cpu_features *f) {
  unsigned max_leaf = cpuid_max_leaf();
  unsigned r[4];
  f->detection_method = "cpuid";
  if (max_leaf < 1) return;

  cpuid(1, 0, r);
  f->sse2 = (r[3] >> 26) & 1u;
  f->sse3 = r[2] & 1u;
  f->ssse3 = (r[2] >> 9) & 1u;
  f->fma = (r[2] >> 12) & 1u;
  f->sse4_1 = (r[2] >> 19) & 1u;
  f->sse4_2 = (r[2] >> 20) & 1u;
  f->avx = (r[2] >> 28) & 1u;
  if ((r[2] >> 27) & 1u) { /* OSXSAVE */
    unsigned xcr0 = xgetbv0();
    f->os_avx = (xcr0 & 0x6u) == 0x6u;
    f->os_avx512 = (xcr0 & 0xE6u) == 0xE6u;
#if defined(__APPLE__)
    /* macOS enables the AVX-512 state lazily on first use, so XCR0 does not
       show it up front; the kernel reports support through sysctl. */
    if (f->os_avx && !f->os_avx512) f->os_avx512 = sysctl_flag("hw.optional.avx512f");
#endif
  }
#if defined(_WIN32)
  /* Sanity check: PF_XSAVE_ENABLED (17). */
  if (!IsProcessorFeaturePresent(17)) f->os_avx = f->os_avx512 = 0;
#endif

  if (max_leaf >= 7) {
    cpuid(7, 0, r);
    f->avx2 = (r[1] >> 5) & 1u;
    f->avx512f = (r[1] >> 16) & 1u;
    f->avx512dq = (r[1] >> 17) & 1u;
    f->avx512bw = (r[1] >> 30) & 1u;
    f->avx512vl = (r[1] >> 31) & 1u;
  }
}

#endif /* RSIMD_CPU_X86 */

/* ---- macOS: sysctl ----------------------------------------------------- */

#if defined(__APPLE__)
/* A missing key (older macOS) counts as unsupported. */
static int sysctl_flag(const char *name) {
  int v = 0;
  size_t len = sizeof v;
  if (sysctlbyname(name, &v, &len, NULL, 0) != 0) return 0;
  return v != 0;
}
#endif

/* ---- arm64 ------------------------------------------------------------- */

/* Linux and FreeBSD hwcap bits, local names so they never clash with (or
   depend on) the system headers. */
#define RSIMD_HWCAP_FPHP (1ul << 9)
#define RSIMD_HWCAP_ASIMDHP (1ul << 10)
#define RSIMD_HWCAP_ASIMDDP (1ul << 20)
#define RSIMD_HWCAP_SVE (1ul << 22)
#define RSIMD_HWCAP2_SVE2 (1ul << 1)
#define RSIMD_HWCAP2_I8MM (1ul << 13)
#define RSIMD_HWCAP2_BF16 (1ul << 14)
#define RSIMD_HWCAP_ARM_NEON (1ul << 12) /* 32-bit arm */

#if defined(RSIMD_CPU_AUXV)
static unsigned long auxv(int which) {
#if defined(__linux__)
  return getauxval((unsigned long) which);
#else
  unsigned long v = 0;
  if (elf_aux_info(which, &v, (int) sizeof v) != 0) return 0;
  return v;
#endif
}
#if defined(__linux__)
#define RSIMD_AUXV_METHOD "getauxval"
#else
#define RSIMD_AUXV_METHOD "elf_aux_info"
#endif
#endif /* RSIMD_CPU_AUXV */

#if defined(RSIMD_CPU_ARM64)

static int sve_vector_length_bits(void) {
#if defined(RSIMD_HAVE_SVE)
  return rsimd_sve_vl_bits_sve();
#elif defined(__linux__)
  int r = prctl(PR_SVE_GET_VL, 0, 0, 0, 0);
  if (r < 0) return 0;
  return (r & PR_SVE_VL_LEN_MASK) * 8;
#else
  return 0;
#endif
}

static void detect_arm64(rsimd_cpu_features *f) {
  f->neon = 1; /* Advanced SIMD is architectural on AArch64 */
#if defined(RSIMD_CPU_AUXV)
  {
    unsigned long hw = auxv(AT_HWCAP);
#if defined(AT_HWCAP2)
    unsigned long hw2 = auxv(AT_HWCAP2);
#else
    unsigned long hw2 = 0;
#endif
    f->detection_method = RSIMD_AUXV_METHOD;
    f->fp16 = (hw & RSIMD_HWCAP_FPHP) && (hw & RSIMD_HWCAP_ASIMDHP);
    f->dotprod = (hw & RSIMD_HWCAP_ASIMDDP) != 0;
    f->sve = (hw & RSIMD_HWCAP_SVE) != 0;
    f->sve2 = (hw2 & RSIMD_HWCAP2_SVE2) != 0;
    f->i8mm = (hw2 & RSIMD_HWCAP2_I8MM) != 0;
    f->bf16 = (hw2 & RSIMD_HWCAP2_BF16) != 0;
  }
#elif defined(__APPLE__)
  /* No Apple Silicon has SVE, and there is no sysctl key for it. */
  f->detection_method = "sysctl";
  f->fp16 = sysctl_flag("hw.optional.arm.FEAT_FP16");
  f->dotprod = sysctl_flag("hw.optional.arm.FEAT_DotProd");
  f->i8mm = sysctl_flag("hw.optional.arm.FEAT_I8MM");
  f->bf16 = sysctl_flag("hw.optional.arm.FEAT_BF16");
#elif defined(_WIN32)
  /* Numeric PF_* constants: MinGW headers may predate the names. */
  f->detection_method = "win32";
  f->dotprod = IsProcessorFeaturePresent(43) != 0; /* PF_ARM_V82_DP_INSTRUCTIONS_AVAILABLE */
  f->sve = IsProcessorFeaturePresent(46) != 0;     /* PF_ARM_SVE_INSTRUCTIONS_AVAILABLE */
  f->sve2 = IsProcessorFeaturePresent(47) != 0;    /* PF_ARM_SVE2_INSTRUCTIONS_AVAILABLE */
  f->i8mm = IsProcessorFeaturePresent(66) != 0;    /* PF_ARM_V82_I8MM_INSTRUCTIONS_AVAILABLE */
  f->fp16 = IsProcessorFeaturePresent(67) != 0;    /* PF_ARM_V82_FP16_INSTRUCTIONS_AVAILABLE */
  f->bf16 = IsProcessorFeaturePresent(68) != 0;    /* PF_ARM_V86_BF16_INSTRUCTIONS_AVAILABLE */
#else
  f->detection_method = "none";
#endif
  if (f->sve) f->sve_vector_length_bits = sve_vector_length_bits();
}

#endif /* RSIMD_CPU_ARM64 */

#if defined(RSIMD_CPU_ARM32)
static void detect_arm32(rsimd_cpu_features *f) {
#if defined(RSIMD_CPU_AUXV)
  f->detection_method = RSIMD_AUXV_METHOD;
  f->neon = (auxv(AT_HWCAP) & RSIMD_HWCAP_ARM_NEON) != 0;
#else
  f->detection_method = "none";
#endif
}
#endif

/* ---- initialisation ---------------------------------------------------- */

static rsimd_cpu_features cpu_unmasked;
static rsimd_cpu_features cpu_final;
static int cpu_done = 0;

void rsimd_cpu_init(void) {
  rsimd_cpu_features f;
  if (cpu_done) return;

  memset(&f, 0, sizeof f);
  f.arch = RSIMD_CPU_ARCH;
  f.os = RSIMD_CPU_OS;
  f.detection_method = "none";
#if defined(RSIMD_CPU_X86)
  detect_x86(&f);
#elif defined(RSIMD_CPU_ARM64)
  detect_arm64(&f);
#elif defined(RSIMD_CPU_ARM32)
  detect_arm32(&f);
#endif
  cpu_make_consistent(&f);
  cpu_unmasked = f;

  cpu_apply_mask(&f, getenv("RSIMD_CPU_FEATURES_MASK"));
  cpu_make_consistent(&f);
  cpu_final = f;
  cpu_done = 1;
}

const rsimd_cpu_features *rsimd_cpu(void) {
  if (!cpu_done) rsimd_cpu_init();
  return &cpu_final;
}

const rsimd_cpu_features *rsimd_cpu_unmasked(void) {
  if (!cpu_done) rsimd_cpu_init();
  return &cpu_unmasked;
}

int rsimd_cpu_supports(rsimd_tier t) {
  const rsimd_cpu_features *f = rsimd_cpu();
  switch (t) {
  case RSIMD_TIER_NONE: return 1;
  case RSIMD_TIER_SSE2: return f->sse2;
  case RSIMD_TIER_AVX2: return f->avx2 && f->fma && f->os_avx;
  case RSIMD_TIER_AVX512:
    return f->avx512f && f->avx512bw && f->avx512dq && f->avx512vl && f->os_avx512;
  case RSIMD_TIER_NEON: return f->neon;
  case RSIMD_TIER_SVE: return f->sve;
  case RSIMD_TIER_SVE2: return f->sve && f->sve2;
  default: return 0; /* reserved tiers */
  }
}
