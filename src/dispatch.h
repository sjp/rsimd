#ifndef RSIMD_DISPATCH_H
#define RSIMD_DISPATCH_H

/* Implementation selection.
 *
 * Every compiled tier (src/tier_<tier>.c) defines one constant table of
 * kernel pointers, rsimd_kernels_<tier>. At load the dispatcher builds a
 * resolved copy of each table the CPU can run, filling every empty slot from
 * the next available tier down the preference order and finally from none,
 * which implements every slot. rsimd_select() then makes one resolved table
 * active, and all kernels are called through rsimd_active.
 *
 * Because of that fill-down, the selected tier (simd_current()) is what was
 * requested, while an individual op may run the kernel of a lower tier.
 *
 * Threads: R evaluates on one thread, and rsimd_active is a plain global
 * written only from .Call entry points. If kernels are ever run on several
 * threads, each .Call must read rsimd_active once into a local on entry and
 * pass that table down. Kernels receive everything through their arguments
 * and read no globals, so the change stays inside the entry points.
 */

#include "tiers.h"

#define RSIMD_OP(name, ret, args) typedef ret (*rsimd_fn_##name) args;
#include "kernel_list.h"
#undef RSIMD_OP

struct rsimd_kernels {
#define RSIMD_OP(name, ret, args) rsimd_fn_##name name;
#include "kernel_list.h"
#undef RSIMD_OP
};

/* Slot indices, in struct order. */
enum {
#define RSIMD_OP(name, ret, args) RSIMD_SLOT_##name,
#include "kernel_list.h"
#undef RSIMD_OP
  RSIMD_SLOT_COUNT
};

/* Slot names, in struct order. */
extern const char *const rsimd_slot_names[RSIMD_SLOT_COUNT];

/* Tiers in "auto" preference order, best first, ending with none. */
extern const rsimd_tier rsimd_tier_preference[];
extern const int rsimd_tier_preference_count;

/* Compiled into the package (none always is). Reserved tiers never are. */
int rsimd_tier_compiled(rsimd_tier t);
/* Compiled and supported by the CPU and operating system. */
int rsimd_tier_available(rsimd_tier t);
/* The tier's own table, before fill-down; NULL unless compiled. */
const struct rsimd_kernels *rsimd_tier_table(rsimd_tier t);
/* The tier whose kernel slot `slot` of tier t's resolved table runs, or
   RSIMD_TIER_COUNT if t is not available. */
rsimd_tier rsimd_slot_source(rsimd_tier t, int slot);

/* Active table and its tier: never NULL once the library is loaded. */
extern const struct rsimd_kernels *rsimd_active;
extern rsimd_tier rsimd_active_tier;
extern int rsimd_active_is_auto; /* 1 if selected as "auto" */

/* Selects "auto" or a tier id. Returns 0 on success, -1 for an unknown
   name and -2 for a tier that is not available; the selection is unchanged
   on failure. Only sets the three globals above. */
int rsimd_select(const char *name);
/* First available tier in preference order. */
rsimd_tier rsimd_best_tier(void);

/* Builds the resolved tables and selects "auto". Called once at load,
   after CPU detection. */
void rsimd_dispatch_init(void);

#endif /* RSIMD_DISPATCH_H */
