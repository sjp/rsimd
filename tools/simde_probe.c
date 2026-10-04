/*
 * Probe translation unit listing every SIMDe header rsimd uses.
 *
 * tools/vendor_simde.sh vendors exactly the include closure of this file, so
 * the subset in src/vendor/simde is derived from it rather than maintained by
 * hand. To use another SIMDe header in a kernel, add it here, re-run
 * tools/vendor_simde.sh and check the size budget it reports.
 *
 * Each #include line must name a header relative to the SIMDe root and stand
 * on a line of its own: the scripts in tools/ read this file textually.
 */

/* x86 base set: sse2, avx2 and neon tiers */
#include "x86/sse2.h"
#include "x86/sse4.1.h"
#include "x86/sse4.2.h"
#include "x86/avx.h"
#include "x86/avx2.h"
#include "x86/fma.h"

/* Curated AVX-512 subset: avx512 tier (f64, i32 and i64 lanes, masks).
 * SIMDe implements only part of AVX-512, so kernels must stay within what
 * these headers provide or guard with #if and fall back to 256-bit code. */
#include "x86/avx512/types.h"
#include "x86/avx512/mov.h"
#include "x86/avx512/loadu.h"
#include "x86/avx512/storeu.h"
#include "x86/avx512/set1.h"
#include "x86/avx512/setzero.h"
#include "x86/avx512/set.h"
#include "x86/avx512/setone.h"
#include "x86/avx512/broadcast.h"
#include "x86/avx512/add.h"
#include "x86/avx512/sub.h"
#include "x86/avx512/mul.h"
#include "x86/avx512/mullo.h"
#include "x86/avx512/div.h"
#include "x86/avx512/fmadd.h"
#include "x86/avx512/fmsub.h"
#include "x86/avx512/fnmadd.h"
#include "x86/avx512/max.h"
#include "x86/avx512/min.h"
#include "x86/avx512/abs.h"
#include "x86/avx512/sqrt.h"
#include "x86/avx512/and.h"
#include "x86/avx512/andnot.h"
#include "x86/avx512/or.h"
#include "x86/avx512/xor.h"
#include "x86/avx512/ternarylogic.h"
#include "x86/avx512/cmp.h"
#include "x86/avx512/cmpeq.h"
#include "x86/avx512/cmpneq.h"
#include "x86/avx512/cmpgt.h"
#include "x86/avx512/cmpge.h"
#include "x86/avx512/cmplt.h"
#include "x86/avx512/cmple.h"
#include "x86/avx512/blend.h"
#include "x86/avx512/mov_mask.h"
#include "x86/avx512/movm.h"
#include "x86/avx512/kand.h"
#include "x86/avx512/knot.h"
#include "x86/avx512/kxor.h"
#include "x86/avx512/test.h"
#include "x86/avx512/cvt.h"
#include "x86/avx512/cvtt.h"
#include "x86/avx512/cvts.h"
#include "x86/avx512/cast.h"
#include "x86/avx512/extract.h"
#include "x86/avx512/insert.h"
#include "x86/avx512/reduce.h"
#include "x86/avx512/round.h"
#include "x86/avx512/roundscale.h"
#include "x86/avx512/permutexvar.h"
#include "x86/avx512/permutex2var.h"
#include "x86/avx512/permutex2var.h"
#include "x86/avx512/shuffle.h"
#include "x86/avx512/unpacklo.h"
#include "x86/avx512/unpackhi.h"
#include "x86/avx512/popcnt.h"
#include "x86/avx512/lzcnt.h"
#include "x86/avx512/sll.h"
#include "x86/avx512/slli.h"
#include "x86/avx512/sllv.h"
#include "x86/avx512/srl.h"
#include "x86/avx512/srli.h"
#include "x86/avx512/srlv.h"
#include "x86/avx512/sra.h"
#include "x86/avx512/srai.h"
#include "x86/avx512/srav.h"
#include "x86/avx512/rol.h"
#include "x86/avx512/ror.h"
#include "x86/avx512/negate.h"
#include "x86/avx512/copysign.h"
#include "x86/avx512/xorsign.h"
#include "x86/avx512/fpclass.h"
#include "x86/avx512/compress.h"
#include "x86/avx512/expand.h"

int main(void) { return 0; }
