/* The none tier: plain C, always built and always usable; the reference
   every other tier is tested against. */
#define RSIMD_TIER none
#include "kernels/common.inc.h"
#include "kernels/tier_info.inc.c"
/* Kernel sources (kernels/<family>.inc.c) are included here. */
#include "kernels/reduce.inc.c"
#include "kernels/scan.inc.c"
#include "kernels/arith.inc.c"
/* After arith.inc.c, whose operand helpers they use. */
#include "kernels/predicates.inc.c"
#include "kernels/complex.inc.c"
#include "kernels/compare.inc.c"
#include "kernels/bitwise.inc.c"
#include "kernels/int64.inc.c"
#include "kernels/convert.inc.c"
#include "kernels/math.inc.c"
#include "kernels/ml.inc.c"
#include "kernels/selftest.inc.c"
#include "kernels/table.inc.c"
