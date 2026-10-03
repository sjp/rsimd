#ifndef RSIMD_SELFTEST_H
#define RSIMD_SELFTEST_H

/* Operation codes of the internal self-test kernel slots (kernel_list.h). */

/* selftest_arith_i32 */
enum {
  RSIMD_ST_ADD,
  RSIMD_ST_SUB,
  RSIMD_ST_MUL,
  RSIMD_ST_NEG,
  RSIMD_ST_ABS,
  RSIMD_ST_ADD_WRAP,
  RSIMD_ST_SUB_WRAP,
  RSIMD_ST_MUL_WRAP,
  RSIMD_ST_NEG_WRAP,
  RSIMD_ST_ABS_WRAP,
  RSIMD_ST_IDIV,
  RSIMD_ST_MOD,
  RSIMD_ST_I32_COUNT
};

/* selftest_arith_f64: pmin and pmax propagate NaN, like base R. */
enum {
  RSIMD_ST_F_ADD,
  RSIMD_ST_F_SUB,
  RSIMD_ST_F_MUL,
  RSIMD_ST_F_DIV,
  RSIMD_ST_F_PMIN,
  RSIMD_ST_F_PMAX,
  RSIMD_ST_F64_COUNT
};

#endif /* RSIMD_SELFTEST_H */
