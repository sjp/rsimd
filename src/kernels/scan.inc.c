/* Prefix scans: cumsum, cumprod, cummin, cummax.
 *
 * Every scan runs in double: integer elements convert exactly, an int32
 * cumsum stays exact because it stops as soon as a value leaves the int32
 * range, and int32 minima and maxima are exact too. A kernel stops at the
 * first missing input (or integer overflow) and returns its index; the
 * entry point fills the rest of the result (NaN, then NA from the first NA
 * onward, for doubles; NA for integers).
 *
 * The none tier runs the sequential loop: the running value combined with
 * each element in turn, as base R does (min and max keep the running value
 * only when it is strictly smaller or larger, so the later of equal values
 * wins, which decides the sign of a zero result). The fixed-width vector
 * tiers scan W lanes at a time in log2(W) shift-and-combine steps (Hillis
 * and Steele), in blocks of four vectors: each vector is scanned, chained
 * to the one before it in the block, combined with the running value
 * broadcast to every lane and stored, and the last lane is the next
 * running value. Each step combines an earlier partial with a later one in
 * that order, so the min/max tie rule and the integer results are exactly
 * the sequential ones; double sums and products are reassociated within
 * each block of 4 W elements and differ from the sequential loop in the
 * last bits. The tail, and a vector that
 * holds a missing input or would overflow, go to the sequential loop.
 *
 * Compensated mode makes cumsum_f64 run a sequential Neumaier sum on every
 * tier (identical results); pairwise mode is treated as fast. cumprod
 * ignores the precision mode, as prod does. The SVE tiers have no scan
 * kernels: the dispatcher uses the neon ones.
 */

enum { RSIMD_SCAN_SUM = 0, RSIMD_SCAN_PROD = 1, RSIMD_SCAN_MIN = 2, RSIMD_SCAN_MAX = 3 };

/* The sequential scan of elements i .. n - 1 (int32 elements when in_i32,
   int32 results when out_i32), continuing from s->f64. Returns the index of
   the first element not written, or -1. */
static inline R_xlen_t RSIMD_KERNEL(scan_seq_)(const void *x, const int in_i32, R_xlen_t i,
                                               R_xlen_t n, void *out, const int out_i32,
                                               const int op, rsimd_scan_state *s) {
  double acc = s->f64;
  for (; i < n; i++) {
    double v;
    if (in_i32) {
      int e = ((const int *) x)[i];
      if (e == RSIMD_NA_I32) break;
      v = (double) e;
    } else {
      v = ((const double *) x)[i];
      if (isnan(v)) break;
    }
    switch (op) {
    case RSIMD_SCAN_SUM: v = acc + v; break;
    case RSIMD_SCAN_PROD: v = acc * v; break;
    case RSIMD_SCAN_MIN: v = acc < v ? acc : v; break;
    default: v = acc > v ? acc : v; break;
    }
    if (out_i32) {
      if (op == RSIMD_SCAN_SUM && !(v >= -(double) INT32_MAX && v <= (double) INT32_MAX)) {
        s->overflow = 1;
        break;
      }
      ((int *) out)[i] = (int) v;
    } else {
      ((double *) out)[i] = v;
    }
    acc = v;
  }
  s->f64 = acc;
  return i < n ? i : -1;
}

/* Compensated cumsum: the Neumaier sum of the elements so far, every tier
   alike. */
static inline R_xlen_t RSIMD_KERNEL(cumsum_comp_)(const double *x, R_xlen_t n, double *out,
                                                  rsimd_scan_state *s) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (isnan(x[i])) return i;
    rsimd_neumaier_add(&s->f64, &s->comp, x[i]);
    out[i] = rsimd_neumaier_value(s->f64, s->comp);
  }
  return -1;
}

#if RSIMD_TIER_IS(none) || RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2) || \
  defined(RSIMD_NO_F64_SIMD)
#define RSIMD_SCAN_(x, in_i32, n, out, out_i32, op, s) \
  RSIMD_KERNEL(scan_seq_)((x), (in_i32), 0, (n), (out), (out_i32), (op), (s))
#else

/* b combined with the earlier partial a. */
RSIMD_ALWAYS_INLINE rsimd_vf64 RSIMD_KERNEL(scan_op_)(rsimd_vf64 a, rsimd_vf64 b, const int op) {
  switch (op) {
  case RSIMD_SCAN_SUM: return rsimd_vf64_add(a, b);
  case RSIMD_SCAN_PROD: return rsimd_vf64_mul(a, b);
  case RSIMD_SCAN_MIN: return rsimd_vf64_min(a, b); /* a < b ? a : b */
  default: return rsimd_vf64_max(a, b);             /* a > b ? a : b */
  }
}

/* The vector scan; see the top of the file. */
RSIMD_ALWAYS_INLINE R_xlen_t RSIMD_KERNEL(scan_vec_)(const void *x, const int in_i32, R_xlen_t n,
                                                    void *out, const int out_i32, const int op,
                                                    rsimd_scan_state *s) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 id = rsimd_vf64_set1(op == RSIMD_SCAN_SUM    ? 0.0
                                        : op == RSIMD_SCAN_PROD ? 1.0
                                        : op == RSIMD_SCAN_MIN  ? HUGE_VAL
                                                                : -HUGE_VAL);
  const rsimd_vf64 na = rsimd_vf64_set1(RSIMD_NA_I32_AS_F64);
  const rsimd_vf64 hi = rsimd_vf64_set1((double) INT32_MAX);
  const rsimd_vf64 lo = rsimd_vf64_set1(-(double) INT32_MAX);
  rsimd_vf64 carry = rsimd_vf64_set1(s->f64);
  R_xlen_t i = 0;
  /* Blocks of four vectors: each vector is scanned on its own, the four are
     chained, and the carry is combined into each once, so the loop-carried
     chain is one op and one broadcast per 4 W elements; one missing-value
     (and overflow) test per block. A block that fails a test is left to
     the vector loop below, which finds the element. */
#define RSIMD_SCAN_LD_(j)                                                                \
  (in_i32 ? rsimd_vf64_loadu_i32((const int32_t *) x + i + (j) * W)                      \
          : rsimd_vf64_loadu((const double *) x + i + (j) * W))
#define RSIMD_SCAN_MISSING_(v) (in_i32 ? rsimd_vf64_cmp_eq((v), na) : rsimd_vf64_is_nan(v))
#define RSIMD_SCAN_OUT_(v) rsimd_mf64_or(rsimd_vf64_cmp_gt((v), hi), rsimd_vf64_cmp_lt((v), lo))
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_vf64 v0 = RSIMD_SCAN_LD_(0), v1 = RSIMD_SCAN_LD_(1);
    rsimd_vf64 v2 = RSIMD_SCAN_LD_(2), v3 = RSIMD_SCAN_LD_(3);
    ptrdiff_t k;
    if (rsimd_mf64_any(rsimd_mf64_or(rsimd_mf64_or(RSIMD_SCAN_MISSING_(v0), RSIMD_SCAN_MISSING_(v1)),
                                     rsimd_mf64_or(RSIMD_SCAN_MISSING_(v2), RSIMD_SCAN_MISSING_(v3))))) {
      break;
    }
    for (k = 1; k < W; k *= 2) {
      v0 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_shift_up(v0, (int) k, id), v0, op);
      v1 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_shift_up(v1, (int) k, id), v1, op);
      v2 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_shift_up(v2, (int) k, id), v2, op);
      v3 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_shift_up(v3, (int) k, id), v3, op);
    }
    v1 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_bcast_last(v0), v1, op);
    v2 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_bcast_last(v1), v2, op);
    v3 = RSIMD_KERNEL(scan_op_)(rsimd_vf64_bcast_last(v2), v3, op);
    v0 = RSIMD_KERNEL(scan_op_)(carry, v0, op);
    v1 = RSIMD_KERNEL(scan_op_)(carry, v1, op);
    v2 = RSIMD_KERNEL(scan_op_)(carry, v2, op);
    v3 = RSIMD_KERNEL(scan_op_)(carry, v3, op);
    if (out_i32) {
      if (op == RSIMD_SCAN_SUM &&
          rsimd_mf64_any(rsimd_mf64_or(rsimd_mf64_or(RSIMD_SCAN_OUT_(v0), RSIMD_SCAN_OUT_(v1)),
                                       rsimd_mf64_or(RSIMD_SCAN_OUT_(v2), RSIMD_SCAN_OUT_(v3))))) {
        break;
      }
      rsimd_vf64_storeu_i32((int32_t *) out + i, v0);
      rsimd_vf64_storeu_i32((int32_t *) out + i + W, v1);
      rsimd_vf64_storeu_i32((int32_t *) out + i + 2 * W, v2);
      rsimd_vf64_storeu_i32((int32_t *) out + i + 3 * W, v3);
    } else {
      rsimd_vf64_storeu((double *) out + i, v0);
      rsimd_vf64_storeu((double *) out + i + W, v1);
      rsimd_vf64_storeu((double *) out + i + 2 * W, v2);
      rsimd_vf64_storeu((double *) out + i + 3 * W, v3);
    }
    carry = rsimd_vf64_bcast_last(v3);
  }
#undef RSIMD_SCAN_LD_
#undef RSIMD_SCAN_MISSING_
#undef RSIMD_SCAN_OUT_
  for (; i + W <= n; i += W) {
    rsimd_vf64 v = in_i32 ? rsimd_vf64_loadu_i32((const int32_t *) x + i)
                          : rsimd_vf64_loadu((const double *) x + i);
    ptrdiff_t k;
    if (rsimd_mf64_any(in_i32 ? rsimd_vf64_cmp_eq(v, na) : rsimd_vf64_is_nan(v))) break;
    for (k = 1; k < W; k *= 2) {
      v = RSIMD_KERNEL(scan_op_)(rsimd_vf64_shift_up(v, (int) k, id), v, op);
    }
    v = RSIMD_KERNEL(scan_op_)(carry, v, op);
    if (out_i32) {
      if (op == RSIMD_SCAN_SUM &&
          rsimd_mf64_any(rsimd_mf64_or(rsimd_vf64_cmp_gt(v, hi), rsimd_vf64_cmp_lt(v, lo)))) {
        break;
      }
      rsimd_vf64_storeu_i32((int32_t *) out + i, v);
    } else {
      rsimd_vf64_storeu((double *) out + i, v);
    }
    carry = rsimd_vf64_bcast_last(v);
  }
  s->f64 = rsimd_vf64_first(carry);
  return RSIMD_KERNEL(scan_seq_)(x, in_i32, i, n, out, out_i32, op, s);
}
#define RSIMD_SCAN_(x, in_i32, n, out, out_i32, op, s) \
  RSIMD_KERNEL(scan_vec_)((x), (in_i32), (n), (out), (out_i32), (op), (s))
#endif

#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2) || defined(RSIMD_NO_F64_SIMD)
/* The neon (or none) kernels are used. */
#define RSIMD_SKIP_cumsum_f64 1
#define RSIMD_SKIP_cumsum_i32 1
#define RSIMD_SKIP_cumprod_f64 1
#define RSIMD_SKIP_cumprod_i32 1
#define RSIMD_SKIP_cumminmax_f64 1
#define RSIMD_SKIP_cumminmax_i32 1
#else

R_xlen_t RSIMD_KERNEL(cumsum_f64)(const double *x, R_xlen_t n, double *out, rsimd_scan_state *s,
                                  const rsimd_opts *o);
R_xlen_t RSIMD_KERNEL(cumsum_f64)(const double *x, R_xlen_t n, double *out, rsimd_scan_state *s,
                                  const rsimd_opts *o) {
  if (o->precision == RSIMD_PREC_COMPENSATED) return RSIMD_KERNEL(cumsum_comp_)(x, n, out, s);
  return RSIMD_SCAN_(x, 0, n, out, 0, RSIMD_SCAN_SUM, s);
}

R_xlen_t RSIMD_KERNEL(cumsum_i32)(const int *x, R_xlen_t n, int *out, rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(cumsum_i32)(const int *x, R_xlen_t n, int *out, rsimd_scan_state *s) {
  return RSIMD_SCAN_(x, 1, n, out, 1, RSIMD_SCAN_SUM, s);
}

R_xlen_t RSIMD_KERNEL(cumprod_f64)(const double *x, R_xlen_t n, double *out,
                                   rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(cumprod_f64)(const double *x, R_xlen_t n, double *out,
                                   rsimd_scan_state *s) {
  return RSIMD_SCAN_(x, 0, n, out, 0, RSIMD_SCAN_PROD, s);
}

R_xlen_t RSIMD_KERNEL(cumprod_i32)(const int *x, R_xlen_t n, double *out, rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(cumprod_i32)(const int *x, R_xlen_t n, double *out, rsimd_scan_state *s) {
  return RSIMD_SCAN_(x, 1, n, out, 0, RSIMD_SCAN_PROD, s);
}

R_xlen_t RSIMD_KERNEL(cumminmax_f64)(const double *x, R_xlen_t n, int max, double *out,
                                     rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(cumminmax_f64)(const double *x, R_xlen_t n, int max, double *out,
                                     rsimd_scan_state *s) {
  return max ? RSIMD_SCAN_(x, 0, n, out, 0, RSIMD_SCAN_MAX, s)
             : RSIMD_SCAN_(x, 0, n, out, 0, RSIMD_SCAN_MIN, s);
}

R_xlen_t RSIMD_KERNEL(cumminmax_i32)(const int *x, R_xlen_t n, int max, int *out,
                                     rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(cumminmax_i32)(const int *x, R_xlen_t n, int max, int *out,
                                     rsimd_scan_state *s) {
  return max ? RSIMD_SCAN_(x, 1, n, out, 1, RSIMD_SCAN_MAX, s)
             : RSIMD_SCAN_(x, 1, n, out, 1, RSIMD_SCAN_MIN, s);
}

#endif

#undef RSIMD_SCAN_
