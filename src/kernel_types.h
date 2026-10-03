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

/* The running state of a reduction. Entry points initialise it with
   rsimd_reduce_result_init(), call the kernel once per chunk to fold that
   chunk in, and convert it to the R result with rsimd_reduce_finish(). */
typedef struct {
  double f64;     /* accumulated value (fast/pairwise/compensated per opts) */
  double comp;    /* compensation term, compensated mode only */
  int64_t i64;    /* integer accumulator (sum/count), int64 min/max */
  R_xlen_t idx;   /* 0-based index of the extremum, -1 if none */
  R_xlen_t count; /* non-NA elements seen (for mean/var, na.rm) */
  unsigned saw_na : 1, saw_nan : 1, overflow : 1, any_true : 1, any_false : 1;
} rsimd_reduce_result;

/* Signature of every elementwise binary kernel, shown for double; other
   element types substitute their C type. A scalar operand (x_scalar or
   y_scalar set) is read from x[0] or y[0] and broadcast. */
typedef void (*rsimd_binop_f64)(const double *x, const double *y, R_xlen_t n, int x_scalar,
                                int y_scalar, double *out, const rsimd_opts *o);

#endif /* RSIMD_KERNEL_TYPES_H */
