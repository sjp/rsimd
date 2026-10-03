/* The avx512 tier: built only if configure's probe for it succeeds, compiled
   with the tier's flags, and entered only after a runtime CPU check. */
#define RSIMD_TIER avx512
#include "kernels/common.inc.h"
#include "kernels/tier_info.inc.c"
/* Kernel sources (kernels/<family>.inc.c) are included here, followed by
   the tier's dispatch table. */
