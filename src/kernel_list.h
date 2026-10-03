/* The kernel slots of struct rsimd_kernels, one line per op/type pair:
 *
 *   RSIMD_OP(slot_name, return_type, (argument list))
 *
 * This file has no include guard: it is included repeatedly with different
 * definitions of RSIMD_OP to generate the function pointer types, the struct
 * members, the slot names and every tier's table (src/dispatch.h,
 * src/kernels/table.inc.c). Adding an op is one line here plus the kernel
 * body, RSIMD_KERNEL(slot_name), in a kernels/<family>.inc.c file.
 *
 * Every tier must define each slot's kernel, except that a kernel source may
 * write `#define RSIMD_SKIP_<slot_name> 1` to leave the slot empty in that
 * tier (an op with no useful SIMD version, or one not written yet); the
 * dispatcher then uses the next tier down. The none tier cannot skip slots.
 */

/* Internal slots. tier_name returns the id of the tier the kernel was
   compiled for. fill_probe exists only in the none tier, so in every other
   tier it is always filled from below. Both are used by the tests of the
   dispatcher. */
RSIMD_OP(tier_name, const char *, (void))
RSIMD_OP(fill_probe, const char *, (void))

/* Internal self-test slots: they apply the NA, precision and overflow
   helpers of na.h to whole chunks so that the tests can compare every
   tier with none through unexported .debug_* functions. op codes are in
   selftest.h. */
RSIMD_OP(selftest_fold_f64, void,
         (const double *x, const double *y, R_xlen_t n, int term, rsimd_reduce_result *r,
          const rsimd_opts *o))
RSIMD_OP(selftest_sum_i32, void,
         (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(selftest_lgl, void,
         (const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(selftest_arith_i32, int,
         (int op, const int *x, const int *y, R_xlen_t n, int x_scalar, int y_scalar, int *out,
          const rsimd_opts *o))
RSIMD_OP(selftest_arith_f64, void,
         (int op, const double *x, const double *y, R_xlen_t n, int x_scalar, int y_scalar,
          double *out, const rsimd_opts *o))

/* Reductions. Each call folds one chunk into the running result (see
   rsimd_reduce_result in kernel_types.h); the entry point merges chunks
   and finishes. Logical vectors use the integer kernels. */
RSIMD_OP(sum_f64, void, (const double *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(sum_i32, void, (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
