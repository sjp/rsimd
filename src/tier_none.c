/* The none tier: plain C, always built and always usable; the reference
   every other tier is tested against. */
#define RSIMD_TIER none
#include "kernels/common.inc.h"
#include "kernels/tier_info.inc.c"
/* Kernel sources (kernels/<family>.inc.c) are included here. */
#include "kernels/reduce.inc.c"
#include "kernels/scan.inc.c"
#include "kernels/selftest.inc.c"
#include "kernels/table.inc.c"
