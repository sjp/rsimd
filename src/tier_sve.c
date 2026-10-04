/* The sve tier: built only if configure's probe for it succeeds, compiled
   with the tier's flags, and entered only after a runtime CPU check. */
#define RSIMD_TIER sve
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
#include "kernels/convert.inc.c"
#include "kernels/math.inc.c"
#include "kernels/ml.inc.c"
#include "kernels/selftest.inc.c"
#include "kernels/table.inc.c"
