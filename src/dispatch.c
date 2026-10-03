#include <stddef.h>
#include <string.h>
#include <R.h>
#include "rsimd_config.h"
#include "cpu_features.h"
#include "dispatch.h"

/* The tables defined by the tier objects (kernels/table.inc.c). Only tiers
   that configure built are referenced, so the others need not exist. */
extern const struct rsimd_kernels rsimd_kernels_none;
#ifdef RSIMD_HAVE_SSE2
extern const struct rsimd_kernels rsimd_kernels_sse2;
#endif
#ifdef RSIMD_HAVE_AVX2
extern const struct rsimd_kernels rsimd_kernels_avx2;
#endif
#ifdef RSIMD_HAVE_AVX512
extern const struct rsimd_kernels rsimd_kernels_avx512;
#endif
#ifdef RSIMD_HAVE_NEON
extern const struct rsimd_kernels rsimd_kernels_neon;
#endif
#ifdef RSIMD_HAVE_SVE
extern const struct rsimd_kernels rsimd_kernels_sve;
#endif
#ifdef RSIMD_HAVE_SVE2
extern const struct rsimd_kernels rsimd_kernels_sve2;
#endif

static const struct rsimd_kernels *const tier_tables[RSIMD_TIER_COUNT] = {
  [RSIMD_TIER_NONE] = &rsimd_kernels_none,
#ifdef RSIMD_HAVE_SSE2
  [RSIMD_TIER_SSE2] = &rsimd_kernels_sse2,
#endif
#ifdef RSIMD_HAVE_AVX2
  [RSIMD_TIER_AVX2] = &rsimd_kernels_avx2,
#endif
#ifdef RSIMD_HAVE_AVX512
  [RSIMD_TIER_AVX512] = &rsimd_kernels_avx512,
#endif
#ifdef RSIMD_HAVE_NEON
  [RSIMD_TIER_NEON] = &rsimd_kernels_neon,
#endif
#ifdef RSIMD_HAVE_SVE
  [RSIMD_TIER_SVE] = &rsimd_kernels_sve,
#endif
#ifdef RSIMD_HAVE_SVE2
  [RSIMD_TIER_SVE2] = &rsimd_kernels_sve2,
#endif
};

const char *const rsimd_slot_names[RSIMD_SLOT_COUNT] = {
#define RSIMD_OP(name, ret, args) #name,
#include "kernel_list.h"
#undef RSIMD_OP
};

#if defined(RSIMD_ARCH_X86)
const rsimd_tier rsimd_tier_preference[] = {
  RSIMD_TIER_AVX512, RSIMD_TIER_AVX2, RSIMD_TIER_SSE2, RSIMD_TIER_NONE
};
#elif defined(RSIMD_ARCH_ARM)
const rsimd_tier rsimd_tier_preference[] = {
  RSIMD_TIER_SVE2, RSIMD_TIER_SVE, RSIMD_TIER_NEON, RSIMD_TIER_NONE
};
#else
const rsimd_tier rsimd_tier_preference[] = {RSIMD_TIER_NONE};
#endif
const int rsimd_tier_preference_count =
  (int) (sizeof(rsimd_tier_preference) / sizeof(rsimd_tier_preference[0]));

/* Fill-down treats a table as an array of slots. All members are function
   pointers, which may be stored in any function pointer type and converted
   back, so the union gives a cast-free view of the slots. */
typedef void (*rsimd_any_fn)(void);
typedef union {
  struct rsimd_kernels k;
  rsimd_any_fn slot[RSIMD_SLOT_COUNT];
} rsimd_table;
_Static_assert(sizeof(struct rsimd_kernels) == RSIMD_SLOT_COUNT * sizeof(rsimd_any_fn),
               "struct rsimd_kernels must contain only function pointers");

static rsimd_table resolved[RSIMD_TIER_COUNT];
static unsigned char slot_source[RSIMD_TIER_COUNT][RSIMD_SLOT_COUNT];
static int available[RSIMD_TIER_COUNT];
static int initialised = 0;

const struct rsimd_kernels *rsimd_active = &rsimd_kernels_none;
rsimd_tier rsimd_active_tier = RSIMD_TIER_NONE;
int rsimd_active_is_auto = 1;

int rsimd_tier_compiled(rsimd_tier t) {
  return t >= 0 && t < RSIMD_TIER_COUNT && tier_tables[t] != NULL;
}

const struct rsimd_kernels *rsimd_tier_table(rsimd_tier t) {
  return rsimd_tier_compiled(t) ? tier_tables[t] : NULL;
}

static void fill_down(void) {
  int p, q, s;
  for (p = 0; p < rsimd_tier_preference_count; p++) {
    rsimd_tier t = rsimd_tier_preference[p];
    if (!available[t]) continue;
    memcpy(&resolved[t].k, tier_tables[t], sizeof(struct rsimd_kernels));
    for (s = 0; s < RSIMD_SLOT_COUNT; s++) {
      slot_source[t][s] = (unsigned char) t;
      /* Empty: take the first lower available tier that has the kernel. */
      for (q = p + 1; resolved[t].slot[s] == NULL && q < rsimd_tier_preference_count; q++) {
        rsimd_tier u = rsimd_tier_preference[q];
        rsimd_table own;
        if (!available[u]) continue;
        memcpy(&own.k, tier_tables[u], sizeof(struct rsimd_kernels));
        resolved[t].slot[s] = own.slot[s];
        slot_source[t][s] = (unsigned char) u;
      }
    }
  }
}

#ifdef RSIMD_DEBUG
/* The none table is the last resort of the fill-down: every slot must be
   set. Missing kernels are a compile error in normal builds; this catches
   anything that slips through. */
static void check_none_table(void) {
  rsimd_table own;
  char missing[1024] = "";
  int s, n = 0;
  memcpy(&own.k, tier_tables[RSIMD_TIER_NONE], sizeof(struct rsimd_kernels));
  for (s = 0; s < RSIMD_SLOT_COUNT; s++) {
    if (own.slot[s] != NULL) continue;
    if (n++) strncat(missing, ", ", sizeof(missing) - strlen(missing) - 1);
    strncat(missing, rsimd_slot_names[s], sizeof(missing) - strlen(missing) - 1);
  }
  if (n) Rf_error("rsimd: the none tier has no kernel for %d slot(s): %s", n, missing);
}
#endif

void rsimd_dispatch_init(void) {
  int t;
  if (initialised) return;
#ifdef RSIMD_DEBUG
  check_none_table();
#endif
  for (t = 0; t < RSIMD_TIER_COUNT; t++) {
    available[t] = rsimd_tier_compiled((rsimd_tier) t) && rsimd_cpu_supports((rsimd_tier) t);
  }
  available[RSIMD_TIER_NONE] = 1;
  fill_down();
  initialised = 1;
  rsimd_select("auto");
}

int rsimd_tier_available(rsimd_tier t) {
  if (t < 0 || t >= RSIMD_TIER_COUNT) return 0;
  if (!initialised) rsimd_dispatch_init();
  return available[t];
}

rsimd_tier rsimd_slot_source(rsimd_tier t, int slot) {
  if (!rsimd_tier_available(t) || slot < 0 || slot >= RSIMD_SLOT_COUNT) {
    return RSIMD_TIER_COUNT;
  }
  return (rsimd_tier) slot_source[t][slot];
}

rsimd_tier rsimd_best_tier(void) {
  int p;
  for (p = 0; p < rsimd_tier_preference_count; p++) {
    if (rsimd_tier_available(rsimd_tier_preference[p])) return rsimd_tier_preference[p];
  }
  return RSIMD_TIER_NONE;
}

int rsimd_select(const char *name) {
  int is_auto = name != NULL && strcmp(name, "auto") == 0;
  rsimd_tier t = is_auto ? rsimd_best_tier() : rsimd_tier_from_name(name);
  if (t == RSIMD_TIER_COUNT) return -1;
  if (!rsimd_tier_available(t)) return -2;
  rsimd_active = &resolved[t].k;
  rsimd_active_tier = t;
  rsimd_active_is_auto = is_auto;
  return 0;
}
