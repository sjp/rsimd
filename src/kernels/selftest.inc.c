/* Kernels of the internal self-test slots (kernel_list.h). They run the
 * helpers of na.h over whole chunks: the none tier through the scalar
 * forms, every other tier through the vector forms, so comparing a tier
 * with none tests the vector form of each rule. They are not user-facing
 * operations, but they follow the kernel contract (no R calls, any n,
 * scalar operands broadcast from x[0] or y[0]).
 */

#include "selftest.h"

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(selftest_fold_f64)(const double *x, const double *y, R_xlen_t n, int term,
                                     rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(selftest_fold_f64)(const double *x, const double *y, R_xlen_t n, int term,
                                     rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_f64(x, y, n, term, r, o);
}

void RSIMD_KERNEL(selftest_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                    const rsimd_opts *o);
void RSIMD_KERNEL(selftest_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                    const rsimd_opts *o) {
  rsimd_fold_sum_i32((const int32_t *) x, n, 0, r, o);
}

void RSIMD_KERNEL(selftest_lgl)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                                const rsimd_opts *o);
void RSIMD_KERNEL(selftest_lgl)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                                const rsimd_opts *o) {
  rsimd_fold_lgl((const int32_t *) x, n, stop, r, o);
}

int RSIMD_KERNEL(selftest_arith_i32)(int op, const int *x, const int *y, R_xlen_t n,
                                     int x_scalar, int y_scalar, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(selftest_arith_i32)(int op, const int *x, const int *y, R_xlen_t n,
                                     int x_scalar, int y_scalar, int *out, const rsimd_opts *o) {
  const int check = o->na_check;
  int ovf = 0;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int32_t a = x[x_scalar ? 0 : i], b = y[y_scalar ? 0 : i], r;
    switch (op) {
    case RSIMD_ST_ADD: r = rsimd_add_i32(a, b, check, &ovf); break;
    case RSIMD_ST_SUB: r = rsimd_sub_i32(a, b, check, &ovf); break;
    case RSIMD_ST_MUL: r = rsimd_mul_i32(a, b, check, &ovf); break;
    case RSIMD_ST_NEG:
    case RSIMD_ST_NEG_WRAP: r = rsimd_neg_i32(a); break;
    case RSIMD_ST_ABS:
    case RSIMD_ST_ABS_WRAP: r = rsimd_abs_i32(a); break;
    case RSIMD_ST_ADD_WRAP: r = rsimd_add_wrap_i32(a, b, check); break;
    case RSIMD_ST_SUB_WRAP: r = rsimd_sub_wrap_i32(a, b, check); break;
    case RSIMD_ST_MUL_WRAP: r = rsimd_mul_wrap_i32(a, b, check); break;
    case RSIMD_ST_IDIV: r = rsimd_intdiv_i32(a, b, 0, check); break;
    case RSIMD_ST_MOD: r = rsimd_intdiv_i32(a, b, 1, check); break;
    default: r = RSIMD_NA_I32; break;
    }
    out[i] = r;
  }
  return ovf;
}

void RSIMD_KERNEL(selftest_arith_f64)(int op, const double *x, const double *y, R_xlen_t n,
                                      int x_scalar, int y_scalar, double *out,
                                      const rsimd_opts *o);
void RSIMD_KERNEL(selftest_arith_f64)(int op, const double *x, const double *y, R_xlen_t n,
                                      int x_scalar, int y_scalar, double *out,
                                      const rsimd_opts *o) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    double a = x[x_scalar ? 0 : i], b = y[y_scalar ? 0 : i], r;
    switch (op) {
    case RSIMD_ST_F_ADD: r = a + b; break;
    case RSIMD_ST_F_SUB: r = a - b; break;
    case RSIMD_ST_F_MUL: r = a * b; break;
    case RSIMD_ST_F_DIV: r = a / b; break;
    case RSIMD_ST_F_PMIN: r = isnan(a) || isnan(b) ? a + b : (a < b ? a : b); break;
    case RSIMD_ST_F_PMAX: r = isnan(a) || isnan(b) ? a + b : (a > b ? a : b); break;
    default: r = 0.0; break;
    }
    out[i] = o->na_check ? rsimd_na_merge_f64(r, a, b) : r;
  }
}

#else /* vector tiers */

#ifndef RSIMD_NO_F64_SIMD
void RSIMD_KERNEL(selftest_fold_f64)(const double *x, const double *y, R_xlen_t n, int term,
                                     rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(selftest_fold_f64)(const double *x, const double *y, R_xlen_t n, int term,
                                     rsimd_reduce_result *r, const rsimd_opts *o) {
  switch (term) {
  case RSIMD_TERM_SQ: rsimd_vfold_f64(x, y, n, RSIMD_TERM_SQ, r, o); break;
  case RSIMD_TERM_ABS: rsimd_vfold_f64(x, y, n, RSIMD_TERM_ABS, r, o); break;
  case RSIMD_TERM_XY: rsimd_vfold_f64(x, y, n, RSIMD_TERM_XY, r, o); break;
  default: rsimd_vfold_f64(x, y, n, RSIMD_TERM_X, r, o); break;
  }
}

void RSIMD_KERNEL(selftest_arith_f64)(int op, const double *x, const double *y, R_xlen_t n,
                                      int x_scalar, int y_scalar, double *out,
                                      const rsimd_opts *o);
void RSIMD_KERNEL(selftest_arith_f64)(int op, const double *x, const double *y, R_xlen_t n,
                                      int x_scalar, int y_scalar, double *out,
                                      const rsimd_opts *o) {
  const rsimd_vf64 bx = rsimd_vf64_set1(x[0]), by = rsimd_vf64_set1(y[0]);
  ptrdiff_t i;
  for (i = 0; i < n; i += RSIMD_LANES_64) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64 a = x_scalar ? bx : rsimd_vf64_loadu_p(pg, x + i, 0.0);
    rsimd_vf64 b = y_scalar ? by : rsimd_vf64_loadu_p(pg, y + i, 0.0);
    rsimd_vf64 r;
    switch (op) {
    case RSIMD_ST_F_ADD: r = rsimd_vf64_add(a, b); break;
    case RSIMD_ST_F_SUB: r = rsimd_vf64_sub(a, b); break;
    case RSIMD_ST_F_MUL: r = rsimd_vf64_mul(a, b); break;
    case RSIMD_ST_F_DIV: r = rsimd_vf64_div(a, b); break;
    case RSIMD_ST_F_PMIN:
    case RSIMD_ST_F_PMAX:
      /* The layer's min/max return b for a NaN operand; propagate it. */
      r = op == RSIMD_ST_F_PMIN ? rsimd_vf64_min(a, b) : rsimd_vf64_max(a, b);
      r = rsimd_vf64_blend(r, rsimd_vf64_add(a, b),
                           rsimd_mf64_or(rsimd_vf64_is_nan(a), rsimd_vf64_is_nan(b)));
      break;
    default: r = rsimd_vf64_zero(); break;
    }
    if (o->na_check) r = rsimd_vf64_na_merge(r, a, b);
    rsimd_vf64_storeu_p(pg, out + i, r);
  }
}
#else
#define RSIMD_SKIP_selftest_fold_f64 1
#define RSIMD_SKIP_selftest_arith_f64 1
#endif

void RSIMD_KERNEL(selftest_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                    const rsimd_opts *o);
void RSIMD_KERNEL(selftest_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                    const rsimd_opts *o) {
  rsimd_vfold_sum_i32((const int32_t *) x, n, 0, r, o);
}

void RSIMD_KERNEL(selftest_lgl)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                                const rsimd_opts *o);
void RSIMD_KERNEL(selftest_lgl)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                                const rsimd_opts *o) {
  rsimd_vfold_lgl((const int32_t *) x, n, stop, r, o);
}

int RSIMD_KERNEL(selftest_arith_i32)(int op, const int *x, const int *y, R_xlen_t n,
                                     int x_scalar, int y_scalar, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(selftest_arith_i32)(int op, const int *x, const int *y, R_xlen_t n,
                                     int x_scalar, int y_scalar, int *out, const rsimd_opts *o) {
  const int32_t *px = (const int32_t *) x, *py = (const int32_t *) y;
  int32_t *pout = (int32_t *) out;
  const int check = o->na_check;
  ptrdiff_t i;

  if (op == RSIMD_ST_IDIV || op == RSIMD_ST_MOD) {
    const int mod = op == RSIMD_ST_MOD;
#ifdef RSIMD_NO_F64_SIMD
    for (i = 0; i < n; i++) {
      pout[i] = rsimd_intdiv_i32(px[x_scalar ? 0 : i], py[y_scalar ? 0 : i], mod, check);
    }
#else
    const rsimd_vf64 bx = rsimd_vf64_set1(px[0]), by = rsimd_vf64_set1(py[0]);
    for (i = 0; i < n; i += RSIMD_LANES_64) {
      rsimd_p64 pg = rsimd_p64_while(i, n);
      rsimd_vf64 a = x_scalar ? bx : rsimd_vf64_loadu_i32_p(pg, px + i, 0);
      rsimd_vf64 b = y_scalar ? by : rsimd_vf64_loadu_i32_p(pg, py + i, 1);
      rsimd_vf64_storeu_i32_p(pg, pout + i, rsimd_vf64_intdiv(a, b, mod, check));
    }
#endif
    return 0;
  }

  {
    const rsimd_vi32 bx = rsimd_vi32_set1(px[0]), by = rsimd_vi32_set1(py[0]);
    rsimd_mi32 ovf = rsimd_mi32_none();
    for (i = 0; i < n; i += RSIMD_LANES_32) {
      rsimd_p32 pg = rsimd_p32_while(i, n);
      rsimd_vi32 a = x_scalar ? bx : rsimd_vi32_loadu_p(pg, px + i, 0);
      rsimd_vi32 b = y_scalar ? by : rsimd_vi32_loadu_p(pg, py + i, 0);
      rsimd_vi32 r;
      rsimd_mi32 bad = rsimd_mi32_none();
      /* 0: NA maps to itself (neg, abs); 1: wrapping; 2: checked. */
      int kind = op >= RSIMD_ST_ADD_WRAP ? 1 : 2;
      switch (op) {
      case RSIMD_ST_ADD:
      case RSIMD_ST_ADD_WRAP: r = rsimd_vi32_add(a, b); break;
      case RSIMD_ST_SUB:
      case RSIMD_ST_SUB_WRAP: r = rsimd_vi32_sub(a, b); break;
      case RSIMD_ST_MUL:
      case RSIMD_ST_MUL_WRAP: r = rsimd_vi32_mul(a, b); break;
      case RSIMD_ST_NEG:
      case RSIMD_ST_NEG_WRAP:
        r = rsimd_vi32_neg_wrap(a);
        kind = 0;
        break;
      case RSIMD_ST_ABS:
      case RSIMD_ST_ABS_WRAP:
        r = rsimd_vi32_abs_wrap(a);
        kind = 0;
        break;
      default:
        r = rsimd_vi32_set1(RSIMD_NA_I32);
        kind = 0;
        break;
      }
      if (kind == 2) {
        if (op == RSIMD_ST_ADD) bad = rsimd_vi32_add_ovf(a, b, r);
        else if (op == RSIMD_ST_SUB) bad = rsimd_vi32_sub_ovf(a, b, r);
        else bad = rsimd_vi32_mul_ovf(a, b, r);
      }
      if (kind != 0) {
        if (check) {
          /* NA operands give NA and never count as overflow; the NA mask
             is applied last. */
          rsimd_mi32 na = rsimd_vi32_na2(a, b);
          ovf = rsimd_mi32_or(ovf, rsimd_mi32_andnot(na, bad));
          bad = rsimd_mi32_or(bad, na);
        } else {
          ovf = rsimd_mi32_or(ovf, bad);
        }
        r = rsimd_vi32_set_na(r, bad);
      }
      rsimd_vi32_storeu_p(pg, pout + i, r);
    }
    return rsimd_mi32_any(ovf);
  }
}

#endif /* vector tiers */
