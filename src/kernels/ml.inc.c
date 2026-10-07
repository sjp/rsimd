/* The passes of softmax and log-softmax (softmax_f64, op codes
 * RSIMD_SOFTMAX_* of kernel_types.h). The entry points materialise the
 * input as doubles, take its maximum m and then run these passes in place,
 * so x and out may be the same buffer: every element is read before it is
 * written, a vector at a time.
 *
 * The exponentials are summed in blocks of RSIMD_PAIRWISE_LEAF elements:
 * each block is computed (into out, or a buffer on the stack for SUM) and
 * then folded with the helpers of sum_f64, so a pairwise sum has the same
 * leaves as simd_sum() of the exponentials. The input holds no NaN and
 * x[i] - m <= 0, so no missing-value checks are needed.
 *
 * The none tier uses libm's exp; the SIMD tiers use SLEEF and leave the
 * slot empty when built without it.
 */

#if RSIMD_TIER_IS(none) || defined(RSIMD_HAVE_SLEEF)

#if RSIMD_TIER_IS(none)
#define RSIMD_ML_FOLD rsimd_fold_f64
#else
#define RSIMD_ML_FOLD rsimd_vfold_f64
#endif

/* dst[i] = exp(x[i] - m) for the n <= RSIMD_PAIRWISE_LEAF elements of a
   block. */
static inline void RSIMD_KERNEL(ml_exp_block_)(const double *x, ptrdiff_t n, double m,
                                               double *dst) {
#if RSIMD_TIER_IS(none)
  ptrdiff_t i;
  for (i = 0; i < n; i++) dst[i] = exp(x[i] - m);
#else
  const rsimd_vf64 vm = rsimd_vf64_set1(m);
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vf64_storeu(dst + i, rsimd_sleef_exp(rsimd_vf64_sub(rsimd_vf64_loadu(x + i), vm)));
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64 a = rsimd_vf64_sub(rsimd_vf64_loadu_p(pg, x + i, m), vm);
    rsimd_vf64_storeu_p(pg, dst + i, rsimd_sleef_exp(a));
  }
#endif
}

void RSIMD_KERNEL(softmax_f64)(int op, const double *x, R_xlen_t n, double m, double c,
                               double *out, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(softmax_f64)(int op, const double *x, R_xlen_t n, double m, double c,
                               double *out, rsimd_reduce_result *r, const rsimd_opts *o) {
  /* The precision mode only: the input has no missing values. */
  const rsimd_opts fo = {0, 0, o->precision, RSIMD_EXT_BOTH};
  double buf[RSIMD_PAIRWISE_LEAF];
  ptrdiff_t i, len;
  switch (op) {
  case RSIMD_SOFTMAX_EXP_SUM:
  case RSIMD_SOFTMAX_SUM:
    for (i = 0; i < n; i += RSIMD_PAIRWISE_LEAF) {
      double *dst = op == RSIMD_SOFTMAX_SUM ? buf : out + i;
      len = n - i < RSIMD_PAIRWISE_LEAF ? n - i : RSIMD_PAIRWISE_LEAF;
      RSIMD_KERNEL(ml_exp_block_)(x + i, len, m, dst);
      RSIMD_ML_FOLD(dst, NULL, len, RSIMD_TERM_X, r, &fo);
    }
    break;
  case RSIMD_SOFTMAX_DIV:
  case RSIMD_SOFTMAX_LOG: {
#if RSIMD_TIER_IS(none)
    if (op == RSIMD_SOFTMAX_DIV) {
      for (i = 0; i < n; i++) out[i] = x[i] / c;
    } else {
      for (i = 0; i < n; i++) out[i] = (x[i] - m) - c;
    }
#else
    const rsimd_vf64 vm = rsimd_vf64_set1(m), vc = rsimd_vf64_set1(c);
    const int div = op == RSIMD_SOFTMAX_DIV;
#define RSIMD_ML_ELT(a) \
  (div ? rsimd_vf64_div((a), vc) : rsimd_vf64_sub(rsimd_vf64_sub((a), vm), vc))
    for (i = 0; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
      rsimd_vf64_storeu(out + i, RSIMD_ML_ELT(rsimd_vf64_loadu(x + i)));
    }
    if (i < n) {
      rsimd_p64 pg = rsimd_p64_while(i, n);
      rsimd_vf64_storeu_p(pg, out + i, RSIMD_ML_ELT(rsimd_vf64_loadu_p(pg, x + i, 1.0)));
    }
#undef RSIMD_ML_ELT
#endif
    break;
  }
  default: break;
  }
}

#undef RSIMD_ML_FOLD

#else /* a SIMD tier without SLEEF: the slot is filled from below */
#define RSIMD_SKIP_softmax_f64 1
#endif
