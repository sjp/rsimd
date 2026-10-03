/* Reduction kernels. Each folds one chunk of the input into the running
 * rsimd_reduce_result with the rules of na.h: precision mode, NA and NaN
 * tracking, na.rm. The none tier uses the scalar folds and no SIMD code at
 * all, so it is the reference the vector tiers are tested against; every
 * other tier uses the vector folds of its own width.
 */

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  rsimd_fold_f64(x, NULL, n, RSIMD_TERM_X, r, o);
}

void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_sum_i32((const int32_t *) x, n, r, o);
}

#else /* vector tiers */

#ifndef RSIMD_NO_F64_SIMD
void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  rsimd_vfold_f64(x, NULL, n, RSIMD_TERM_X, r, o);
}
#else
#define RSIMD_SKIP_sum_f64 1
#endif

void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_vfold_sum_i32((const int32_t *) x, n, r, o);
}

#endif /* vector tiers */
