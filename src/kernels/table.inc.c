/* Instantiates the tier's dispatch table, rsimd_kernels_<tier>. Included
 * once at the end of every tier translation unit, after the kernel sources.
 *
 * Each slot of struct rsimd_kernels (src/kernel_list.h) is set to
 * RSIMD_KERNEL(slot), unless the kernel sources defined RSIMD_SKIP_<slot>
 * as 1, in which case it is NULL and the dispatcher fills it from a lower
 * tier. The none tier must implement every slot: a skip there is a compile
 * error, and so is a slot whose kernel was never written.
 *
 * Test builds: RSIMD_TEST_HOLE_<slot> (from configure's RSIMD_TEST_HOLE)
 * empties the slot in every tier but none; RSIMD_TEST_NONE_HOLE_<slot>
 * empties it in none, to exercise the debug check of the none table.
 */

#include "rsimd_config.h"
#include "dispatch.h"

/* RSIMD_IS_SET(MACRO) is 1 if MACRO is defined as 1, else 0. */
#define RSIMD_PH_1 0,
#define RSIMD_TAKE2(a, b, ...) b
#define RSIMD_IS_SET_3(junk) RSIMD_TAKE2(junk 1, 0, ~)
#define RSIMD_IS_SET_2(val) RSIMD_IS_SET_3(RSIMD_PH_##val)
#define RSIMD_IS_SET(x) RSIMD_IS_SET_2(x)

#define RSIMD_OR(a, b) RSIMD_OR_(a, b)
#define RSIMD_OR_(a, b) RSIMD_OR_##a##b
#define RSIMD_OR_00 0
#define RSIMD_OR_01 1
#define RSIMD_OR_10 1
#define RSIMD_OR_11 1

/* RSIMD_IF(c)(then, else), for c 0 or 1; the branch not taken is dropped
   before compilation, so it may name a kernel that does not exist. */
#define RSIMD_IF(c) RSIMD_IF_(c)
#define RSIMD_IF_(c) RSIMD_IF_##c
#define RSIMD_IF_0(t, f) f
#define RSIMD_IF_1(t, f) t

#define RSIMD_SKIPPED(name) RSIMD_IS_SET(RSIMD_SKIP_##name)

#if RSIMD_TIER_IS(none)
#define RSIMD_OP(name, ret, args) \
  _Static_assert(!RSIMD_SKIPPED(name), "the none tier must implement slot " #name);
#include "kernel_list.h"
#undef RSIMD_OP
#define RSIMD_EMPTY(name) RSIMD_IS_SET(RSIMD_TEST_NONE_HOLE_##name)
#else
#define RSIMD_EMPTY(name) RSIMD_OR(RSIMD_SKIPPED(name), RSIMD_IS_SET(RSIMD_TEST_HOLE_##name))
#endif

extern const struct rsimd_kernels RSIMD_CAT(rsimd_kernels, RSIMD_TIER);
const struct rsimd_kernels RSIMD_CAT(rsimd_kernels, RSIMD_TIER) = {
#define RSIMD_OP(name, ret, args) .name = RSIMD_IF(RSIMD_EMPTY(name))(NULL, RSIMD_KERNEL(name)),
#include "kernel_list.h"
#undef RSIMD_OP
};
