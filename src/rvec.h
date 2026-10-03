#ifndef RSIMD_RVEC_H
#define RSIMD_RVEC_H

/* R vector access layer.
 *
 * Every .Call entry point turns its SEXP arguments into kernel calls through
 * this layer: it classifies and validates inputs, reads them in contiguous
 * chunks whether or not they are ALTREP, applies the length-1 broadcast
 * rule, runs the chunk loop with interrupt checks, allocates results and
 * applies the attribute policy. Entry points and kernels never touch R
 * vector data or lengths directly; tools/lint_c.sh enforces that only
 * rvec.c uses the raw accessors (DATAPTR, REAL(), LENGTH() ...) and
 * R_CheckUserInterrupt.
 *
 * Errors are raised with Rf_error, which inside .Call reports the user's
 * call (the R function that made the .Call), so no call object is passed.
 */

#include <stddef.h>
#include "kernel_types.h"
#include "na.h"

/* ---- Element types ----------------------------------------------------- */

typedef enum {
  RSIMD_F64,  /* double */
  RSIMD_I32,  /* integer (factors rejected) */
  RSIMD_LGL,  /* logical, int storage */
  RSIMD_U8,   /* raw */
  RSIMD_I64,  /* bit64::integer64: a double vector holding int64_t bits */
  RSIMD_C128, /* complex */
  RSIMD_BAD
} rsimd_etype;

/* "double", "integer", "logical", "raw", "integer64", "complex",
   "unsupported", indexed by rsimd_etype. */
extern const char *const rsimd_etype_names[RSIMD_BAD + 1];

/* The element type of x, or RSIMD_BAD for anything that is not a plain
   atomic vector of a supported type (factors, data frames, lists, S4
   objects, character ...). integer64 is recognised by its class, so bit64
   need not be loaded. */
rsimd_etype rsimd_etype_of(SEXP x);
/* As rsimd_etype_of(), but errors for RSIMD_BAD, naming the argument:
   "'x' must be an atomic vector (double, integer, logical, raw, complex or
   integer64), not <class>". */
rsimd_etype rsimd_check_atomic(SEXP x, const char *arg);
/* Bytes per element. */
size_t rsimd_etype_size(rsimd_etype t);

/* Common type of the operands of a binary operation, following base R
   arithmetic: LGL -> I32 -> F64; I64 with I64, I32 or LGL -> I64 (bit64);
   I64 with F64 -> F64 (the R side warns "integer64 coerced to double");
   C128 with F64, I32 or LGL -> C128. U8 combines only with U8 (bitwise
   ops; arithmetic entry points reject raw themselves). Returns RSIMD_BAD
   when the pair cannot be combined.
   Two logical operands stay RSIMD_LGL: their storage is int, so kernels
   treat them as I32 without a copy, and entry points pick the result type
   (base R makes TRUE + TRUE an integer). Coercion itself happens R-side,
   only for an operand whose type differs from the result. */
rsimd_etype rsimd_promote(rsimd_etype a, rsimd_etype b);
/* 1 if combining a and b converts integer64 to double. */
int rsimd_promote_warns(rsimd_etype a, rsimd_etype b);

/* ---- Chunked read access ----------------------------------------------- */

/* Elements per region on the ALTREP path; the region buffer is at most
   64 KiB (complex). A multiple of RSIMD_PAIRWISE_LEAF. */
#define RSIMD_CHUNK 4096

typedef struct {
  SEXP sx;
  rsimd_etype type;
  R_xlen_t n;
  const void *ptr; /* non-NULL => contiguous fast path */
  int no_na_hint;  /* 1 if known NA-free (ALTREP *_NO_NA, or raw) */
} rsimd_in;

/* Classifies x (erroring as rsimd_check_atomic() with name `arg`) and finds
   its data pointer without materialising ALTREP objects. Returns 0. */
int rsimd_in_init(rsimd_in *v, SEXP x, const char *arg);

/* A pointer to elements [i, i + *len). On the contiguous path that is the
   data itself and *len = n - i; on the ALTREP path the elements are copied
   into `buf` (caller-owned, RSIMD_CHUNK elements of the type) and
   *len = min(RSIMD_CHUNK, n - i). Requires 0 <= i < n. Inputs are never
   materialised: the data pointer comes from DATAPTR_OR_NULL. */
const void *rsimd_in_region(const rsimd_in *v, R_xlen_t i, R_xlen_t *len, void *buf);

/* ---- Chunk loop with interrupts ----------------------------------------- */

/* Default number of elements between interrupt checks. */
#define RSIMD_INTERRUPT_STRIDE ((R_xlen_t) 1 << 20)

/* The stride in effect: RSIMD_INTERRUPT_STRIDE, or the value of the
   environment variable RSIMD_DEBUG_STRIDE (a positive integer) when the
   library was loaded, rounded up to a multiple of RSIMD_PAIRWISE_LEAF
   (128). A small stride exercises chunking in tests. Every chunk boundary
   is a multiple of 128 (RSIMD_CHUNK is one too), which pairwise summation
   relies on to give the same result however the input is chunked. */
extern R_xlen_t rsimd_stride;

/* Reads RSIMD_DEBUG_STRIDE; called once at load. */
void rsimd_rvec_init(void);

/* R_CheckUserInterrupt(). It may longjmp, so it is only called between
   chunks, where every SEXP the entry point allocated is PROTECTed (entry
   points PROTECT their result before the loop). */
void rsimd_check_interrupt(void);

/* Runs the statements in `...` once per chunk of input `in` (an rsimd_in
   pointer) with
     const T *px      the chunk's elements,
     R_xlen_t len     its length (1 <= len <= rsimd_stride, and
                      <= RSIMD_CHUNK on the ALTREP path),
     R_xlen_t off     the index of its first element.
   T must be the C type of in->type (double for F64 and I64, int for I32
   and LGL, Rbyte for U8, Rcomplex for C128). Reductions fold each chunk
   into their rsimd_reduce_result; elementwise ops write the result at off.
   The body may `break` out of the loop (early exit) but must not use
   `continue`. Interrupts are checked after every rsimd_stride elements. */
#define RSIMD_FOREACH_CHUNK(in, T, px, len, off, ...)                                  \
  do {                                                                                 \
    T rsimd_buf_[RSIMD_CHUNK];                                                         \
    R_xlen_t rsimd_tick_ = 0;                                                          \
    for (R_xlen_t off = 0; off < (in)->n;) {                                           \
      R_xlen_t len;                                                                    \
      const T *px = (const T *) rsimd_in_region((in), off, &len, rsimd_buf_);          \
      if (len > rsimd_stride) len = rsimd_stride;                                      \
      __VA_ARGS__                                                                      \
      off += len;                                                                      \
      rsimd_tick_ += len;                                                              \
      if (rsimd_tick_ >= rsimd_stride) {                                               \
        rsimd_tick_ = 0;                                                               \
        rsimd_check_interrupt();                                                       \
      }                                                                                \
    }                                                                                  \
  } while (0)

/* ---- Binary operands and broadcast -------------------------------------- */

typedef struct {
  rsimd_in x, y;
  R_xlen_t n;              /* result length */
  int y_scalar, x_scalar;  /* operand has length 1 and is broadcast */
} rsimd_bin;

/* Length-1 broadcast only: equal lengths give n = nx; ny == 1 sets
   y_scalar (n = nx); otherwise nx == 1 sets x_scalar (n = ny). Any other
   pair errors: "lengths of 'x' (<nx>) and 'y' (<ny>) must be equal or one
   of them must be 1". So (0, 0), (0, 1) and (1, 0) give n = 0, and (0, 3)
   is an error. Both operands must already have the common type (the R side
   promotes them). Returns 0. */
int rsimd_bin_init(rsimd_bin *b, SEXP x, SEXP y);

/* As RSIMD_FOREACH_CHUNK for both operands of `b` (an rsimd_bin pointer):
   px and py point at the chunk's elements of x and y. A scalar operand's
   pointer addresses its single element on every chunk (kernels read x[0]
   when x_scalar is set). */
#define RSIMD_FOREACH_CHUNK2(b, T, px, py, len, off, ...)                              \
  do {                                                                                 \
    T rsimd_bufx_[RSIMD_CHUNK], rsimd_bufy_[RSIMD_CHUNK];                              \
    R_xlen_t rsimd_tick_ = 0, rsimd_l_;                                                \
    const T *rsimd_sx_ =                                                               \
        (b)->x_scalar ? (const T *) rsimd_in_region(&(b)->x, 0, &rsimd_l_, rsimd_bufx_) \
                      : NULL;                                                          \
    const T *rsimd_sy_ =                                                               \
        (b)->y_scalar ? (const T *) rsimd_in_region(&(b)->y, 0, &rsimd_l_, rsimd_bufy_) \
                      : NULL;                                                          \
    for (R_xlen_t off = 0; off < (b)->n;) {                                            \
      R_xlen_t len = (b)->n - off;                                                     \
      const T *px = rsimd_sx_, *py = rsimd_sy_;                                        \
      if (len > rsimd_stride) len = rsimd_stride;                                      \
      if (px == NULL) {                                                                \
        px = (const T *) rsimd_in_region(&(b)->x, off, &rsimd_l_, rsimd_bufx_);        \
        if (rsimd_l_ < len) len = rsimd_l_;                                            \
      }                                                                                \
      if (py == NULL) {                                                                \
        py = (const T *) rsimd_in_region(&(b)->y, off, &rsimd_l_, rsimd_bufy_);        \
        if (rsimd_l_ < len) len = rsimd_l_;                                            \
      }                                                                                \
      __VA_ARGS__                                                                      \
      off += len;                                                                      \
      rsimd_tick_ += len;                                                              \
      if (rsimd_tick_ >= rsimd_stride) {                                               \
        rsimd_tick_ = 0;                                                               \
        rsimd_check_interrupt();                                                       \
      }                                                                                \
    }                                                                                  \
  } while (0)

/* ---- Reductions --------------------------------------------------------- */

/* Reductions with a scalar result. simd_range() finishes the min and the
   max of one min/max fold and combines them; scans (cumsum ...) and
   which_na return vectors and do not use rsimd_reduce_finish(). */
typedef enum {
  RSIMD_RED_SUM,
  RSIMD_RED_PROD,
  RSIMD_RED_MEAN,
  RSIMD_RED_MIN,
  RSIMD_RED_MAX,
  RSIMD_RED_WHICH_MIN,
  RSIMD_RED_WHICH_MAX,
  RSIMD_RED_ANY,
  RSIMD_RED_ALL,
  RSIMD_RED_ANY_NA,
  RSIMD_RED_COUNT_NA,
  RSIMD_RED_SUM_SQ,
  RSIMD_RED_SUM_ABS,
  RSIMD_RED_DOT,
  RSIMD_RED_NORM,
  RSIMD_RED_DIST,
  RSIMD_RED_COSINE,
  RSIMD_RED_VAR,
  RSIMD_RED_SD,
  RSIMD_RED_OP_COUNT
} rsimd_reduce_op;

/* "sum", "prod", ..., indexed by rsimd_reduce_op. */
extern const char *const rsimd_reduce_op_names[RSIMD_RED_OP_COUNT];

/* The identity of `op`: zero accumulators, idx = -1, no flags, except
   prod (f64 = 1), min and which_min (f64 = Inf, i64 = INT64_MAX), max and
   which_max (f64 = -Inf, i64 = INT64_MIN). For every op the maximum of
   the min/max kernels starts at f64_hi = -Inf, i64_hi = INT64_MIN. */
void rsimd_reduce_result_init(rsimd_reduce_result *r, int op);

/* The R value of a finished reduction over n elements of type `type`, with
   base R's result types:
     sum        I32/LGL -> integer, or double when the 64-bit accumulator
                exceeds INT_MAX in magnitude (or overflowed into f64);
                F64 -> double; I64 -> integer64
     prod, mean, var, sd, sum_sq, sum_abs, dot, norm, dist, cosine
                -> double
     min, max   I32/LGL -> integer; F64 -> double; I64 -> integer64
     which_min, which_max
                1-based idx + 1 (integer(0) when idx < 0); integer, or
                double when n > INT_MAX; also for U8
     any, all, any_na -> logical
     count_na   i64, as a double
   The double value is rsimd_reduce_value() for o->precision (f64, the
   compensated pair or the pairwise leaf tree). Errors with
   "invalid 'type' (<type>) of argument" for a type the op does not take.

   Missing values and empty input (na.h describes the rules):
     - with na.rm = FALSE, saw_na gives NA of the result type, else saw_nan
       gives NaN, for every op but which_*, any, all, any_na and count_na;
     - any/all use three-valued logic from any_true, any_false, saw_na;
     - count is the number of surviving elements: mean of none is NaN,
       var and sd of fewer than two are NA, min and max of none are Inf and
       -Inf (double, whatever the type) with base R's warning "no
       non-missing arguments to min; returning Inf". integer64 min and max
       are not covered by the empty rule yet. */
SEXP rsimd_reduce_finish(int op, rsimd_etype type, R_xlen_t n, const rsimd_reduce_result *r,
                         const rsimd_opts *o);

/* c(lo, hi) for simd_range(), from the finished min and max (scalars):
   integer when both are integer, else double. */
SEXP rsimd_range_pair(SEXP lo, SEXP hi);

/* Base R's warning for checked integer arithmetic, "NAs produced by integer
   overflow"; entry points call it once per call when a kernel reported
   overflow. */
void rsimd_warn_int_overflow(void);

/* ---- Results ------------------------------------------------------------ */

/* A bare vector of n elements of type t (REALSXP for F64 and I64), not
   PROTECTed. */
SEXP rsimd_alloc_like(rsimd_etype t, R_xlen_t n);
/* Writable data of a result allocated by rsimd_alloc_like(). */
void *rsimd_out_ptr(SEXP out);
/* The attribute policy: results are bare, except that a class attribute of
   exactly "integer64" (copied only to a double result) or "simd_vec" is
   copied from x to out. For "simd_vec" the "rsimd_impl" attribute is copied
   too, and "rsimd_no_na" only when no_na is non-zero (the op guarantees an
   NA-free result); otherwise it is left unset. */
void rsimd_copy_class(SEXP x, SEXP out, int no_na);

/* ---- Options and arguments ---------------------------------------------- */

/* Fills o from the entry point's arguments; pass R_NilValue for an argument
   the op does not take (na.rm -> FALSE, na_check -> TRUE, precision ->
   fast). na_check is cleared when no_na_hint is set: an input known to be
   NA-free needs no check even when one was asked for. precision is the
   integer code RSIMD_PREC_*. */
void rsimd_opts_init(rsimd_opts *o, SEXP na_rm, SEXP na_check, SEXP precision, int no_na_hint);

/* Scalar argument readers; each errors naming the argument:
   "'<name>' must be TRUE or FALSE", "... a single integer",
   "... a single number", "... a single string". */
int rsimd_arg_lgl1(SEXP x, const char *name);
int rsimd_arg_int1(SEXP x, const char *name);
double rsimd_arg_dbl1(SEXP x, const char *name);
const char *rsimd_arg_str(SEXP x, const char *name);

/* ---- Complex ------------------------------------------------------------ */

/* The n complex values as 2n doubles (re, im, re, im ...), without a copy,
   for ops that treat the real and imaginary parts alike. */
static inline const double *rsimd_c128_as_f64(const Rcomplex *x, R_xlen_t n) {
  (void) n;
  return (const double *) x;
}

#endif /* RSIMD_RVEC_H */
