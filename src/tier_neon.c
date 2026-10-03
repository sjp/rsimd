/* The neon tier: built only if configure's probe for it succeeds, compiled
   with the tier's flags, and entered only after a runtime CPU check. */
#define RSIMD_TIER neon
#include "kernels/common.inc.h"
#include "kernels/tier_info.inc.c"
/* Kernel sources (kernels/<family>.inc.c) are included here. */
#include "kernels/selftest.inc.c"
#include "kernels/table.inc.c"
