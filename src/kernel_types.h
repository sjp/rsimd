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

/* Elementary functions (kernels/math.inc.c). Unary op codes for
   math1_f64; LOGB is log(x) / p, the logarithm to a base whose natural
   logarithm is p; SIGMOID is the logistic function 1 / (1 + exp(-x)),
   computed as e = exp(-|x|), then 1 / (1 + e) for x >= 0 and e / (1 + e)
   otherwise, which is accurate in both tails (subnormal results down to
   x = -745). Binary op codes for math2_f64: pow(x, y) with base R's
   rules for x ^ y, atan2(y, x) and hypot(x, y). The kernels return
   RSIMD_EW_NAN_PRODUCED when a result is NaN where no operand was (pow
   never sets it, as base R's ^ never warns). */
enum {
  RSIMD_MATH_EXP = 0,
  RSIMD_MATH_EXP2,
  RSIMD_MATH_EXP10,
  RSIMD_MATH_EXPM1,
  RSIMD_MATH_LOG,
  RSIMD_MATH_LOG2,
  RSIMD_MATH_LOG10,
  RSIMD_MATH_LOG1P,
  RSIMD_MATH_LOGB,
  RSIMD_MATH_CBRT,
  RSIMD_MATH_SIN,
  RSIMD_MATH_COS,
  RSIMD_MATH_TAN,
  RSIMD_MATH_ASIN,
  RSIMD_MATH_ACOS,
  RSIMD_MATH_ATAN,
  RSIMD_MATH_SINPI,
  RSIMD_MATH_COSPI,
  RSIMD_MATH_TANPI,
  RSIMD_MATH_SINH,
  RSIMD_MATH_COSH,
  RSIMD_MATH_TANH,
  RSIMD_MATH_ASINH,
  RSIMD_MATH_ACOSH,
  RSIMD_MATH_ATANH,
  RSIMD_MATH_SIGMOID
};
enum { RSIMD_MATH_POW = 0, RSIMD_MATH_ATAN2, RSIMD_MATH_HYPOT };

/* The passes of softmax and log-softmax (kernels/ml.inc.c), op codes of
   softmax_f64 over x, which may be the same buffer as out:
     EXP_SUM  out[i] = exp(x[i] - m), and the sum of those added to r;
     SUM      only adds the sum of exp(x[i] - m) to r (out unused);
     DIV      out[i] = x[i] / c;
     LOG      out[i] = (x[i] - m) - c.
   The sums are in the precision mode, folding the exponentials in blocks
   of RSIMD_PAIRWISE_LEAF as sum_f64 would fold them. */
enum { RSIMD_SOFTMAX_EXP_SUM = 0, RSIMD_SOFTMAX_SUM, RSIMD_SOFTMAX_DIV, RSIMD_SOFTMAX_LOG };

/* Predicates (kernels/predicates.inc.c): is.na (NA or NaN), is.nan (NaN
   but not NA), is.finite, is.infinite, negative (sign bit set and not
   NaN; for integers x < 0 and not NA) and zero (x == 0, so also -0).
   pred_i32 takes NA, FINITE (not NA), NEGATIVE and ZERO. */
enum {
  RSIMD_PRED_NA = 0,
  RSIMD_PRED_NAN,
  RSIMD_PRED_FINITE,
  RSIMD_PRED_INFINITE,
  RSIMD_PRED_NEGATIVE,
  RSIMD_PRED_ZERO
};
/* What a predicate kernel computes: the logical vector, or whether any or
   all elements satisfy it. */
enum { RSIMD_PRED_ELT = 0, RSIMD_PRED_ANY, RSIMD_PRED_ALL };

/* Comparisons (kernels/compare.inc.c), x op y. */
enum { RSIMD_CMP_EQ = 0, RSIMD_CMP_NE, RSIMD_CMP_LT, RSIMD_CMP_LE, RSIMD_CMP_GT, RSIMD_CMP_GE };

/* Three-valued logic (kernels/bitwise.inc.c): x & y, x | y, xor(x, y), !x. */
enum { RSIMD_LOGIC_AND = 0, RSIMD_LOGIC_OR, RSIMD_LOGIC_XOR, RSIMD_LOGIC_NOT };

/* Bitwise ops (kernels/bitwise.inc.c). AND, OR and XOR are binary; the
   others take x and a count k: SHL, SHR (logical) and SAR (arithmetic)
   shift by k, ROTL and ROTR rotate by k, POPCNT, LZCNT and TZCNT count
   bits (k unused). */
enum {
  RSIMD_BIT_AND = 0,
  RSIMD_BIT_OR,
  RSIMD_BIT_XOR,
  RSIMD_BIT_NOT,
  RSIMD_BIT_SHL,
  RSIMD_BIT_SHR,
  RSIMD_BIT_SAR,
  RSIMD_BIT_ROTL,
  RSIMD_BIT_ROTR,
  RSIMD_BIT_POPCNT,
  RSIMD_BIT_LZCNT,
  RSIMD_BIT_TZCNT
};

/* Type conversions (kernels/convert.inc.c): source and destination
   element types (I32 is integer or logical storage, LGL the logical result
   0/1/NA), and the mode of the conversions to integer and raw: CHECKED
   (out of range gives NA or 00, and a warning), SATURATING (clamped) or
   TRUNCATING (two's complement wrap of the truncated value). */
enum {
  RSIMD_CVT_F64_I32 = 0,
  RSIMD_CVT_F64_U8,
  RSIMD_CVT_F64_LGL,
  RSIMD_CVT_I32_F64,
  RSIMD_CVT_I32_U8,
  RSIMD_CVT_I32_LGL,
  RSIMD_CVT_U8_F64,
  RSIMD_CVT_U8_I32,
  RSIMD_CVT_U8_LGL
};
enum { RSIMD_CVT_CHECKED = 0, RSIMD_CVT_SATURATING, RSIMD_CVT_TRUNCATING };
/* Status bits of the conversion kernels, for base R's coercion warnings
   (in checked mode only): a double outside the integer range ("NAs
   introduced by coercion to integer range"), a value not in 0..255 or
   missing ("out-of-range values treated as 0 in coercion to raw"). */
enum { RSIMD_CVT_WARN_INT = 1, RSIMD_CVT_WARN_RAW = 2 };

#endif /* RSIMD_KERNEL_TYPES_H */
