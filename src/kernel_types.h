#ifndef RSIMD_KERNEL_TYPES_H
#define RSIMD_KERNEL_TYPES_H

/* Plain C types shared by the kernels and the R access layer (rvec.h).
 *
 * Kernel signatures in kernel_list.h use only these types, scalars and
 * const pointers, so the tier translation units never see the SEXP-handling
 * code. <Rinternals.h> is included only for R_xlen_t; R_NO_REMAP keeps its
 * short names (length, error ...) out of the SIMD headers when this is the
 * first R header of a translation unit.
 */

#include <stdint.h>
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>

/* Accumulation modes for floating-point reductions (option rsimd.precision). */
enum { RSIMD_PREC_FAST = 0, RSIMD_PREC_PAIRWISE = 1, RSIMD_PREC_COMPENSATED = 2 };

typedef struct {
  int na_rm;     /* from na.rm= */
  int na_check;  /* from na_check= / option rsimd.na_check, cleared when the
                    inputs are known NA-free */
  int precision; /* RSIMD_PREC_FAST / PAIRWISE / COMPENSATED */
} rsimd_opts;

/* What the missing-value kernels (na_f64 ...) compute. */
enum {
  RSIMD_NAMODE_ANY = 0,
  RSIMD_NAMODE_COUNT = 1,
  RSIMD_NAMODE_WHICH_I32 = 2,
  RSIMD_NAMODE_WHICH_F64 = 3
};

/* Partial sums of pairwise summation: s[k] is the sum of 2^k leaves for
   every set bit k of `leaves` (see na.h). */
typedef struct {
  double s[64];
  uint64_t leaves;
} rsimd_pairwise;

/* The running state of a reduction. Entry points initialise it with
   rsimd_reduce_result_init(), call the kernel once per chunk to fold that
   chunk in, and convert it to the R result with rsimd_reduce_finish(). */
typedef struct {
  double f64;       /* accumulated value (fast and compensated modes); the
                       minimum for min/max kernels */
  double comp;      /* compensation term, compensated mode only */
  rsimd_pairwise pw; /* leaf sums, pairwise mode only */
  int64_t i64;      /* integer accumulator (sum/count), integer minimum */
  double f64_hi;    /* the maximum for min/max kernels */
  int64_t i64_hi;   /* integer maximum */
  R_xlen_t idx;     /* 0-based index of the extremum, -1 if none */
  R_xlen_t count;   /* elements that survived na.rm (all of them without) */
  /* saw_na: an NA was seen; saw_nan: a NaN (NA included) was seen;
     overflow: an integer sum left int64 and continues in f64. */
  unsigned saw_na : 1, saw_nan : 1, overflow : 1, any_true : 1, any_false : 1;
} rsimd_reduce_result;

/* Element types of the two operands of dot, dist and cosine: both double,
   both int32 (integer or logical), or double x with int32 y (entry points
   swap a pair the other way round, as these ops are symmetric). */
enum { RSIMD_PAIR_F64_F64 = 0, RSIMD_PAIR_I32_I32 = 1, RSIMD_PAIR_F64_I32 = 2 };

/* The running state of a prefix scan (cumsum ...) between chunks: the
   last value (the sum, product, minimum or maximum so far), the Neumaier
   compensation of a compensated cumsum, and whether an integer cumsum
   overflowed. Entry points start f64 at the identity: 0, 1, Inf or -Inf. */
typedef struct {
  double f64;
  double comp;
  unsigned overflow : 1;
} rsimd_scan_state;

/* Elementwise kernels (kernels/arith.inc.c) take an op code and flags.

   Binary op codes, for ew2_f64 (ADD .. COPYSIGN) and ew2_i32 (ADD .. MOD,
   PMIN .. PMAX_NUM and the _WRAP ops). PMIN and PMAX propagate missing
   values like base R's pmin/pmax: when both operands are missing the
   second one's kind (NA or NaN) wins. PMIN_NUM and PMAX_NUM ignore a
   missing operand (base R's na.rm = TRUE), and give the second operand
   when both are missing. */
enum {
  RSIMD_EW_ADD = 0,
  RSIMD_EW_SUB,
  RSIMD_EW_MUL,
  RSIMD_EW_DIV,
  RSIMD_EW_IDIV,
  RSIMD_EW_MOD,
  RSIMD_EW_PMIN,
  RSIMD_EW_PMAX,
  RSIMD_EW_PMIN_NUM,
  RSIMD_EW_PMAX_NUM,
  RSIMD_EW_COPYSIGN,
  RSIMD_EW_ADD_WRAP,
  RSIMD_EW_SUB_WRAP,
  RSIMD_EW_MUL_WRAP
};
/* Ternary op codes (ew3_f64, ew3_i32): fma(x, y, z) rounds once; mul_add
   is (x * y) + z and add_mul (x + y) * z, each step rounded (checked for
   integers, as base R's composition); lerp(x, y, t) is
   fma(t, y, (1 - t) * x); clamp(x, lo, hi) is pmin(pmax(x, lo), hi).
   ew3_i32 takes MUL_ADD, ADD_MUL and CLAMP. */
enum { RSIMD_EW_FMA = 0, RSIMD_EW_MUL_ADD, RSIMD_EW_ADD_MUL, RSIMD_EW_LERP, RSIMD_EW_CLAMP };
/* Unary op codes: ew1_f64 takes them all (ROUND is half to even);
   ew1_i32 takes NEG and ABS, which wrap and so map NA to itself. */
enum {
  RSIMD_EW_NEG = 0,
  RSIMD_EW_ABS,
  RSIMD_EW_SIGN,
  RSIMD_EW_RECIP,
  RSIMD_EW_SQRT,
  RSIMD_EW_FLOOR,
  RSIMD_EW_CEIL,
  RSIMD_EW_TRUNC,
  RSIMD_EW_ROUND
};

/* Operand flags of the elementwise kernels, for operand k (0 = x, 1 = y,
   2 = z): RSIMD_EW_SCALAR(k) broadcasts element 0 of the operand;
   RSIMD_EW_I32(k), for the f64 kernels, says the operand holds int32
   elements (integer or logical, NA becoming NA_real_) rather than
   doubles. */
#define RSIMD_EW_SCALAR(k) (1 << (k))
#define RSIMD_EW_I32(k) (8 << (k))

/* Status bits returned by the elementwise kernels: a checked integer op
   overflowed (outside the NA lanes), sqrt made NaN from a number, clamp
   saw lo > hi. */
enum { RSIMD_EW_OVERFLOW = 1, RSIMD_EW_NAN_PRODUCED = 2, RSIMD_EW_LO_GT_HI = 4 };

#endif /* RSIMD_KERNEL_TYPES_H */
