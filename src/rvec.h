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

/* 1 for the types with 32-bit int storage (integer and logical). */
static inline int rsimd_is_int_like(rsimd_etype t) { return t == RSIMD_I32 || t == RSIMD_LGL; }

/* The element type of x, or RSIMD_BAD for anything that is not a plain
   atomic vector of a supported type (classed vectors such as factors and
   dates, data frames, lists, S4 objects, character ...). The classes
   taken are integer64, recognised by its class so that bit64 need not be
   loaded, and simd_vec. */
rsimd_etype rsimd_etype_of(SEXP x);
/* As rsimd_etype_of(), but errors for RSIMD_BAD, naming the argument:
   "'x' must be an atomic vector (double, integer, logical, raw, complex or
   integer64), not <class>". */
rsimd_etype rsimd_check_atomic(SEXP x, const char *arg);

/* Rf_error() with the user's call into rsimd (.user_call()) rather than
   the call of the closure that made the .Call(): from a simd_vec method,
   the generic's call (sum(x), not simd_sum(x, na.rm = na.rm)). */
void rsimd_error(const char *fmt, ...)
#ifdef __GNUC__
  __attribute__((noreturn, format(printf, 1, 2)))
#endif
  ;

/* Errors for an operand of type t that the op does not take. From an
   exported simd_* function the message is rsimd's "<fun>() does not
   support '<arg>' of type <type>", naming the argument that has type t in
   the function's frame (`arg` when it does, or none is found; NULL if
   unknown), with the user's call; from anywhere else (a simd_vec method
   of a base generic) it is base R's `msg` (NULL for "invalid 'type'
   (<type>) of argument"), with the generic's call. Evaluates
   .type_error() in the namespace. */
void rsimd_type_error(rsimd_etype t, const char *arg, const char *msg)
#ifdef __GNUC__
  __attribute__((noreturn))
#endif
  ;
/* The same for operands of types t1 and t2 that the op does not take
   together: "<fun>() cannot combine <t1> and <t2> operands", or base R's
   `msg`. */
void rsimd_mix_error(rsimd_etype t1, rsimd_etype t2, const char *msg)
#ifdef __GNUC__
  __attribute__((noreturn))
#endif
  ;
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
   its data pointer without materialising ALTREP objects; records a simd_vec
   operand (see rsimd_entry()). Returns 0. */
int rsimd_in_init(rsimd_in *v, SEXP x, const char *arg);

/* A pointer to elements [i, i + *len). On the contiguous path that is the
   data itself and *len = n - i; on the ALTREP path the elements are copied
   into `buf` (caller-owned, RSIMD_CHUNK elements of the type) and
   *len = min(RSIMD_CHUNK, n - i). Requires 0 <= i < n. Inputs are never
   materialised: the data pointer comes from DATAPTR_OR_NULL. */
const void *rsimd_in_region(const rsimd_in *v, R_xlen_t i, R_xlen_t *len, void *buf);

/* The endpoints of a compact sequence of base R (1:n, seq_len(n) and the
   like, integer or double, step 1 or -1) that is not expanded: its
   elements are first, first + step, ..., last, integers of magnitude
   below 2^53. Returns 1 and fills *first and *last for such an input of
   type RSIMD_I32 or RSIMD_F64, 0 for any other (always 0 before R 4.6.0,
   which has no API naming an ALTREP class). */
int rsimd_in_seq(const rsimd_in *v, double *first, double *last);

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
#define RSIMD_FOREACH_CHUNK(in, T, px, len, off, ...) \
  RSIMD_FOREACH_CHUNK_FROM(in, T, 0, px, len, off, __VA_ARGS__)

/* As RSIMD_FOREACH_CHUNK over the elements from index `start` on. */
#define RSIMD_FOREACH_CHUNK_FROM(in, T, start, px, len, off, ...)                      \
  do {                                                                                 \
    T rsimd_buf_[RSIMD_CHUNK];                                                         \
    R_xlen_t rsimd_tick_ = 0;                                                          \
    for (R_xlen_t off = (start); off < (in)->n;) {                                     \
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

/* Length-1 broadcast only: equal lengths give n = nx; a zero-length
   operand gives n = 0 (as in base R, so (0, 3) is not an error); ny == 1
   sets y_scalar (n = nx); otherwise nx == 1 sets x_scalar (n = ny). Any
   other pair errors: "lengths of 'x' (<nx>) and 'y' (<ny>) must be equal
   or one of them must be 1". Both operands must already have the common type (the R side
   promotes them). Returns 0. */
int rsimd_bin_init(rsimd_bin *b, SEXP x, SEXP y);

/* As RSIMD_FOREACH_CHUNK for both operands of `b` (an rsimd_bin pointer):
   px and py point at the chunk's elements of x and y. A scalar operand's
   pointer addresses its single element on every chunk (kernels read x[0]
   when x_scalar is set). */
#define RSIMD_FOREACH_CHUNK2(b, T, px, py, len, off, ...) \
  RSIMD_FOREACH_CHUNK2T(b, T, T, px, py, len, off, __VA_ARGS__)

/* As RSIMD_FOREACH_CHUNK2 for operands of different element types: TX is
   the C type of x's elements, TY that of y's. */
#define RSIMD_FOREACH_CHUNK2T(b, TX, TY, px, py, len, off, ...)                        \
  do {                                                                                 \
    TX rsimd_bufx_[RSIMD_CHUNK];                                                       \
    TY rsimd_bufy_[RSIMD_CHUNK];                                                       \
    R_xlen_t rsimd_tick_ = 0, rsimd_l_;                                                \
    const TX *rsimd_sx_ =                                                              \
        (b)->x_scalar ? (const TX *) rsimd_in_region(&(b)->x, 0, &rsimd_l_, rsimd_bufx_) \
                      : NULL;                                                          \
    const TY *rsimd_sy_ =                                                              \
        (b)->y_scalar ? (const TY *) rsimd_in_region(&(b)->y, 0, &rsimd_l_, rsimd_bufy_) \
                      : NULL;                                                          \
    for (R_xlen_t off = 0; off < (b)->n;) {                                            \
      R_xlen_t len = (b)->n - off;                                                     \
      const TX *px = rsimd_sx_;                                                        \
      const TY *py = rsimd_sy_;                                                        \
      if (len > rsimd_stride) len = rsimd_stride;                                      \
      if (px == NULL) {                                                                \
        px = (const TX *) rsimd_in_region(&(b)->x, off, &rsimd_l_, rsimd_bufx_);       \
        if (rsimd_l_ < len) len = rsimd_l_;                                            \
      }                                                                                \
      if (py == NULL) {                                                                \
        py = (const TY *) rsimd_in_region(&(b)->y, off, &rsimd_l_, rsimd_bufy_);       \
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

/* ---- Elementwise operands ----------------------------------------------- */

#define RSIMD_EW_MAX_ARGS 3

/* The operands of an elementwise op (one to three), each of any element
   type: entry points pick the kernel and its RSIMD_EW_I32(k) flags from
   in[k].type. */
typedef struct {
  int k;                       /* number of operands */
  rsimd_in in[RSIMD_EW_MAX_ARGS];
  int flags;                   /* RSIMD_EW_SCALAR(k) bits: operand k is broadcast */
  R_xlen_t n;                  /* result length */
  int no_na_hint;              /* every operand is known NA-free */
} rsimd_ew;

/* Classifies the k operands args[0 .. k - 1], named names[0 .. k - 1] in
   messages, and applies the length-1 broadcast rule: the operands that do
   not have length 1 must all have the same length, which is the result
   length (1 if every operand has length 1), and the operands of length 1
   are broadcast when it is not 1. A zero-length operand makes the result
   zero-length whatever the other lengths, as in base R. Otherwise errors, for two operands as
   rsimd_bin_init() does, for three with "lengths of 'x' (<nx>), 'y' (<ny>)
   and 'z' (<nz>) must be equal or 1". Returns 0. */
int rsimd_ew_init(rsimd_ew *e, int k, const SEXP *args, const char *const *names);

/* Errors (an internal error) if an operand of e is complex. */
void rsimd_ew_no_c128(const rsimd_ew *e);

/* Runs the statements in `...` once per chunk of the operands of `e` (an
   rsimd_ew pointer) with
     const void *p[RSIMD_EW_MAX_ARGS]  the chunk's elements of each operand
                                       (a broadcast operand's single
                                       element on every chunk),
     R_xlen_t len, off                 as in RSIMD_FOREACH_CHUNK.
   The body may `break` out of the loop (early exit) but must not use
   `continue`. Interrupts are checked after every rsimd_stride elements. The chunk buffers hold RSIMD_CHUNK elements
   of up to 8 bytes, so a complex operand is an internal error: callers
   divert complex operands before the loop. */
#define RSIMD_FOREACH_CHUNK_EW(e, p, len, off, ...)                                    \
  do {                                                                                 \
    _Alignas(Rcomplex) unsigned char rsimd_ebuf_[RSIMD_EW_MAX_ARGS]                    \
                                                [RSIMD_CHUNK * sizeof(double)];        \
    R_xlen_t rsimd_tick_ = 0, rsimd_l_;                                                \
    const void *p[RSIMD_EW_MAX_ARGS] = {NULL, NULL, NULL};                             \
    int rsimd_k_;                                                                      \
    rsimd_ew_no_c128(e);                                                               \
    for (rsimd_k_ = 0; rsimd_k_ < (e)->k; rsimd_k_++) {                                \
      if ((e)->flags & RSIMD_EW_SCALAR(rsimd_k_)) {                                    \
        p[rsimd_k_] = rsimd_in_region(&(e)->in[rsimd_k_], 0, &rsimd_l_, rsimd_ebuf_[rsimd_k_]); \
      }                                                                                \
    }                                                                                  \
    for (R_xlen_t off = 0; off < (e)->n;) {                                            \
      R_xlen_t len = (e)->n - off;                                                     \
      if (len > rsimd_stride) len = rsimd_stride;                                      \
      for (rsimd_k_ = 0; rsimd_k_ < (e)->k; rsimd_k_++) {                              \
        if (!((e)->flags & RSIMD_EW_SCALAR(rsimd_k_))) {                               \
          p[rsimd_k_] = rsimd_in_region(&(e)->in[rsimd_k_], off, &rsimd_l_, rsimd_ebuf_[rsimd_k_]); \
          if (rsimd_l_ < len) len = rsimd_l_;                                          \
        }                                                                              \
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
                F64 -> double; I64 -> integer64 (exact, or NA)
     prod, mean, var, sd, sum_sq, sum_abs, dot, norm, dist, cosine
                -> double
     min, max   I32/LGL -> integer; F64 -> double; I64 -> integer64
     which_min, which_max
                1-based idx + 1 (integer(0) when idx < 0); integer, or
                double when n > INT_MAX; also for U8
     any, all, any_na -> logical
     count_na   i64, as a double
   The double value is rsimd_reduce_value() for o->precision (f64, the
   compensated pair or the pairwise leaf tree). Errors (rsimd_type_error())
   for a type the op does not take.

   Missing values and empty input (na.h describes the rules):
     - with na.rm = FALSE, saw_na gives NA of the result type, else saw_nan
       gives NaN, for every op but which_*, any, all, any_na and count_na;
       var and sd give NA for either (base R);
     - any/all use three-valued logic from any_true, any_false, saw_na;
     - count is the number of surviving elements: mean of none is NaN,
       var and sd of fewer than two are NA, min and max of none are Inf and
       -Inf (double) with base R's warning "no non-missing arguments to
       min; returning Inf", except for integer64, where they are +INT64_MAX
       and -INT64_MAX with bit64's warning "no non-NA value, returning the
       highest possible integer64 value +9223372036854775807" (lowest ...
       -9223372036854775807);
     - an integer64 sum whose exact total (i64 + carry * 2^64) is outside
       int64 or is INT64_MIN gives NA with bit64's overflow warning. */
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
/* Sets class "integer64" on out (a double vector holding int64_t bits,
   PROTECTed by the caller). Loads bit64's namespace first, so that its S3
   methods are registered, unless it is known to be loaded (this evaluates
   R code); errors, with an install hint, if bit64 is not installed. */
void rsimd_set_i64_class(SEXP out);
/* A length-1 integer64 vector holding v. */
SEXP rsimd_scalar_i64(int64_t v);
/* bit64's warning for integer64 overflow, "NAs produced by integer64
   overflow". */
void rsimd_warn_i64_overflow(void);

/* ---- simd_vec operands -------------------------------------------------- */

/* Every .Call entry point that runs kernels calls rsimd_entry() first. It
   resets rsimd_active to the selected table and forgets the simd_vec
   operands of the previous call (which may have ended in an error with a
   pinned table active).

   rsimd_in_init() (and so rsimd_bin_init() and rsimd_ew_init()) then
   records each simd_vec operand (one whose class contains "simd_vec"):
     - a valid NA-free flag (rsimd_sv_flag() == 1) sets the operand's
       no_na_hint, so na_check is skipped as for an ALTREP that promises
       no NA;
     - attribute rsimd_impl (a tier name) pins the call: rsimd_active
       becomes that tier's resolved table. It errors for an invalid value,
       for a tier that is not available on this machine and for operands
       pinned to different tiers.
   A plain operand of length 1 is checked for NA directly, so that a scalar
   does not stop an NA-free call from skipping its checks. */
void rsimd_entry(void);
/* rsimd_entry() first honours a change to option rsimd.impl made without
   simd_use() (by evaluating .sync_option() when the option differs from the
   value it last synced with); a selection (C_simd_select) calls
   rsimd_impl_selected() so that the next call checks again. */
void rsimd_impl_selected(void);

/* Warnings are deferred to the end of the call. Rf_warning() can run R
   code (a calling handler), which could call rsimd again and overwrite
   the per-call state above while this call still needs it. rsimd_warn()
   (printf-style; the message is truncated to 255 bytes, and only the first
   8 of a call are kept) records a warning, and every entry point returns
   through rsimd_exit(out), which issues the recorded warnings in order
   and returns out, with the user's call (as rsimd_error()). A call that
   records a warning and then errors must
   call rsimd_warn_flush() before the error. rsimd_entry() discards any
   warnings left by a previous call. */
void rsimd_warn(const char *fmt, ...)
#ifdef __GNUC__
  __attribute__((format(printf, 1, 2)))
#endif
  ;
void rsimd_warn_flush(void);
SEXP rsimd_exit(SEXP out);

/* The NA-free flag of a simd_vec lives in attribute rsimd_na_token, an
   external pointer to the object itself whose protected value is the flag
   (TRUE or FALSE), so that R code cannot forge it. Base R copies
   attributes onto new data (pmin(), storage.mode<-, Re(), ...), so the
   flag counts only while the token points at x and is not shared: a copy
   shares the token, which raises its reference count for good (GC never
   lowers it), so the copy and the original both read as unknown, and a
   later object at a reused address cannot pick the flag up either. The
   token's tag is the class attribute x had when it was stamped, and the
   flag counts only while x still has that very vector: class(x) <- NULL
   on an unshared x does not copy it, so a primitive [<- could then change
   the data in place, but putting the class back installs another vector:
   the stamp gives x a class vector no other object has. A token is also
   never valid after unserialize() (the pointer comes back NULL). Tokens
   and classes are read without Rf_getAttrib(), which marks the value it
   returns as shared.

   rsimd_sv_flag() returns 1 (known NA-free), 0 (known to contain a
   missing value) or -1 (unknown) for x; raw data is always 1.
   rsimd_sv_stamp() sets the flag of x in place to flag (1, 0, or -1 to
   remove it) with a fresh token and returns x; x must not be shared. */
int rsimd_sv_flag(SEXP x);
SEXP rsimd_sv_stamp(SEXP x, int flag);
/* Removes in place, and returns, a token of x that points at another
   object, together with the flag; see C_simd_sv_release(). */
SEXP rsimd_sv_release(SEXP x);
/* A copy of the data of atomic x without attributes (memcpy), or
   R_NilValue for a type it does not handle or an ALTREP x without a data
   pointer. Copying x itself would share its token. */
SEXP rsimd_sv_bare(SEXP x);

/* For entry points whose result is a value vector (not a reduction, mask or
   index): when an operand of the call was a simd_vec, out (a fresh result,
   or an input, which is then shallow-copied) becomes a simd_vec with class
   "simd_vec", or c("simd_vec", "integer64") for an integer64 result,
   attribute rsimd_impl set to the call's pin (absent when unpinned) and
   the NA-free flag set (rsimd_sv_stamp()) when the op cannot make a missing value from
   non-missing input (keeps_na_free) and every operand was known NA-free,
   or when out is raw; otherwise the flag is absent. Other results are
   returned unchanged. out need not be PROTECTed; the value returned is not
   PROTECTed. */
SEXP rsimd_sv_result(SEXP out, int keeps_na_free);

/* ---- Options and arguments ---------------------------------------------- */

/* Fills o from the entry point's arguments, with precision fast (set
   o->precision from rsimd_arg_precision() for an op that takes it). Pass
   R_NilValue for na.rm when the op does not take it (FALSE), and for
   na_check to read option rsimd.na_check (TRUE when unset); an op whose
   na.rm or na_check is fixed calls rsimd_opts_init_fixed() instead.
   na_check is cleared when no_na_hint is set: an input known to be NA-free
   needs no check even when one was asked for. */
void rsimd_opts_init(rsimd_opts *o, SEXP na_rm, SEXP na_check, int no_na_hint);
/* As rsimd_opts_init() with na.rm and na_check given as 0 or 1. */
void rsimd_opts_init_fixed(rsimd_opts *o, int na_rm, int na_check, int no_na_hint);
/* The na.rm flag from an argument holding one, FALSE for R_NilValue. */
int rsimd_arg_na_rm(SEXP na_rm);
/* The precision code RSIMD_PREC_* from an argument holding one, or from
   option rsimd.precision for R_NilValue (fast when unset). */
int rsimd_arg_precision(SEXP precision);
/* The math accuracy code (0 accurate, 1 fast) from an argument holding
   one, or from option rsimd.math_accuracy for R_NilValue (accurate when
   unset). */
int rsimd_arg_accuracy(SEXP accuracy);

/* Scalar argument readers; each errors naming the argument:
   "'<name>' must be TRUE or FALSE", "... a single integer",
   "... a single number", "... a single string". */
int rsimd_arg_lgl1(SEXP x, const char *name);
int rsimd_arg_int1(SEXP x, const char *name);
double rsimd_arg_dbl1(SEXP x, const char *name);
/* As rsimd_arg_dbl1(), but a missing value (NA, NaN or a logical NA) is
   allowed and returned as NA_real_ or NaN. */
double rsimd_arg_num1(SEXP x, const char *name);
const char *rsimd_arg_str(SEXP x, const char *name);
/* 1 if x is a single string equal to one of names (NULL-terminated). */
int rsimd_str_in(SEXP x, const char *const *names);
/* The index of name in table (count entries; NULL entries are skipped).
   Errors "unknown <what> '<name>'" when it is not there. */
int rsimd_lookup(const char *name, const char *const *table, int count, const char *what);
/* The index of the op name op (a single string from the R side) in names
   (count entries); an unknown name is an internal error. */
int rsimd_arg_choice(SEXP op, const char *const *names, int count);

/* ---- Complex ------------------------------------------------------------ */

/* The n complex values as 2n doubles (re, im, re, im ...), without a copy,
   for ops that treat the real and imaginary parts alike. */
static inline const double *rsimd_c128_as_f64(const Rcomplex *x, R_xlen_t n) {
  (void) n;
  return (const double *) x;
}

#endif /* RSIMD_RVEC_H */
