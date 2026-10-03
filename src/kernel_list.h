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
/* The sum of x - c, in the precision mode: the refinement pass of mean. */
RSIMD_OP(sum_dev_f64, void,
         (const double *x, R_xlen_t n, double c, rsimd_reduce_result *r, const rsimd_opts *o))
/* Products, in double; the precision mode does not apply. */
RSIMD_OP(prod_f64, void, (const double *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(prod_i32, void, (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
/* Minimum into f64 (i64) and maximum into f64_hi (i64_hi), over the
   elements that are not missing; count is the number of those. A zero
   extremum keeps the sign of the first zero in the input, as base R. */
RSIMD_OP(minmax_f64, void,
         (const double *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(minmax_i32, void, (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
/* The 0-based index of the first element equal to v (numerically, so 0
   matches -0), or -1. Used by which_min and which_max after min/max. */
RSIMD_OP(find_f64, R_xlen_t, (const double *x, R_xlen_t n, double v))
RSIMD_OP(find_i32, R_xlen_t, (const int *x, R_xlen_t n, int v))
/* which_min (max = 0) or which_max (max = 1) of raw bytes: the extremum in
   i64 and its index in idx, off being the chunk's offset. */
RSIMD_OP(which_u8, void,
         (const Rbyte *x, R_xlen_t n, R_xlen_t off, int max, rsimd_reduce_result *r))
/* any/all flags (any_true, any_false, saw_na) of the elements as logical
   values, stopping early as `stop` (RSIMD_STOP_*) says. A double is TRUE
   when non-zero and NA when NaN; an integer TRUE when non-zero; a raw byte
   TRUE when non-zero. */
RSIMD_OP(anyall_f64, void,
         (const double *x, R_xlen_t n, int stop, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(anyall_i32, void,
         (const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(anyall_u8, void, (const Rbyte *x, R_xlen_t n, int stop, rsimd_reduce_result *r))
/* Missing values (NA or NaN; for complex, either part), by `mode`
   (RSIMD_NAMODE_*): ANY sets saw_na at the first one and returns; COUNT
   adds their number to i64; WHICH_I32 and WHICH_F64 write their 1-based
   indices (off + i + 1) to out[i64], out[i64 + 1] ... as int or double,
   advancing i64. */
RSIMD_OP(na_f64, void,
         (const double *x, R_xlen_t n, int mode, R_xlen_t off, void *out, rsimd_reduce_result *r))
RSIMD_OP(na_i32, void,
         (const int *x, R_xlen_t n, int mode, R_xlen_t off, void *out, rsimd_reduce_result *r))
RSIMD_OP(na_c128, void,
         (const Rcomplex *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
          rsimd_reduce_result *r))

/* Sums of x^2 and |x| in the precision mode; the _i32 kernels read
   integer and logical elements as doubles (an NA as NA_real_), except
   sumabs_i32, which sums exactly in i64 as sum_i32 does. */
RSIMD_OP(sumsq_f64, void,
         (const double *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(sumsq_i32, void, (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(sumabs_f64, void,
         (const double *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(sumabs_i32, void, (const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o))
/* Sums of x*y (dot) and (x - y)^2 (dist) over pairs, in the precision
   mode, for operands of the element types `types` (RSIMD_PAIR_*); a pair
   is missing when either element is, and na.rm removes the pair. */
RSIMD_OP(dot_f64, void,
         (const void *x, const void *y, R_xlen_t n, int types, rsimd_reduce_result *r,
          const rsimd_opts *o))
RSIMD_OP(dist_f64, void,
         (const void *x, const void *y, R_xlen_t n, int types, rsimd_reduce_result *r,
          const rsimd_opts *o))
/* cosine's three sums in one pass over the chunk: x*y into r[0], x^2 into
   r[1] and y^2 into r[2]. */
RSIMD_OP(cosine_f64, void,
         (const void *x, const void *y, R_xlen_t n, int types, rsimd_reduce_result *r,
          const rsimd_opts *o))
/* The sum of (x - c)^2 in the precision mode: the second pass of var. */
RSIMD_OP(var_pass2_f64, void,
         (const double *x, R_xlen_t n, double c, rsimd_reduce_result *r, const rsimd_opts *o))
RSIMD_OP(var_pass2_i32, void,
         (const int *x, R_xlen_t n, double c, rsimd_reduce_result *r, const rsimd_opts *o))

/* Prefix scans. Each writes out[i] for the chunk's elements, continuing
   from the state s of the previous chunk, and returns the index of the
   first element it did not write, or -1 when it wrote all n: the first
   NaN (NA included) of a double input, the first NA of an integer input,
   or the element where an integer cumsum leaves the int32 range (which
   also sets s->overflow). The entry point fills the rest. cumsum_f64 uses
   a sequential compensated sum in compensated mode and a vector scan
   otherwise; cumminmax takes the minimum (max = 0) or maximum (max = 1),
   the later element winning ties as in base R. */
RSIMD_OP(cumsum_f64, R_xlen_t,
         (const double *x, R_xlen_t n, double *out, rsimd_scan_state *s, const rsimd_opts *o))
RSIMD_OP(cumsum_i32, R_xlen_t, (const int *x, R_xlen_t n, int *out, rsimd_scan_state *s))
RSIMD_OP(cumprod_f64, R_xlen_t, (const double *x, R_xlen_t n, double *out, rsimd_scan_state *s))
RSIMD_OP(cumprod_i32, R_xlen_t, (const int *x, R_xlen_t n, double *out, rsimd_scan_state *s))
RSIMD_OP(cumminmax_f64, R_xlen_t,
         (const double *x, R_xlen_t n, int max, double *out, rsimd_scan_state *s))
RSIMD_OP(cumminmax_i32, R_xlen_t,
         (const int *x, R_xlen_t n, int max, int *out, rsimd_scan_state *s))
