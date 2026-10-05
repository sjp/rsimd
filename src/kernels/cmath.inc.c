/* Elementwise math on complex numbers: Mod and Arg (math1_c128).
 *
 * Base R computes Mod(z) with libm's cabs and Arg(z) with carg, which are
 * hypot(re, im) and atan2(im, re). The none tier calls those, so it is
 * base R on every platform, special values included (glibc's hypot, for
 * one, returns NA rather than Inf for NA and an infinite part, because R's
 * NA is a signalling NaN). The SIMD tiers use SLEEF's hypot (0.5 ULP) and
 * atan2 (1 ULP; 3.5 ULP for both with RSIMD_MATH_FAST) for elements with
 * both parts finite and libm for the others.
 *
 * Included after complex.inc.c (rsimd_c128_load_) and math.inc.c.
 */

/* Mod or Arg of one element with libm. */
static inline double rsimd_cmath_1(int op, Rcomplex z) {
  return op == RSIMD_CMATH_MOD ? hypot(z.r, z.i) : atan2(z.i, z.r);
}

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(math1_c128)(int op, const Rcomplex *x, R_xlen_t n, double *out);
void RSIMD_KERNEL(math1_c128)(int op, const Rcomplex *x, R_xlen_t n, double *out) {
  R_xlen_t i;
  op &= ~RSIMD_MATH_FAST;
  for (i = 0; i < n; i++) out[i] = rsimd_cmath_1(op, x[i]);
}

#elif !defined(RSIMD_NO_F64_SIMD) && defined(RSIMD_HAVE_SLEEF)

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(math1_c128_)(const int op, const int fast, const Rcomplex *x,
                                                   R_xlen_t n, double *out) {
  const ptrdiff_t W = RSIMD_LANES_64;
  ptrdiff_t i = 0, j;
  for (; i + W <= n; i += W) {
    rsimd_vf64 re, im, r;
    rsimd_c128_load_(x, i, &re, &im);
    if (op == RSIMD_CMATH_MOD) r = fast ? rsimd_sleef_hypot_fast(re, im) : rsimd_sleef_hypot(re, im);
    else r = fast ? rsimd_sleef_atan2_fast(im, re) : rsimd_sleef_atan2(im, re);
    rsimd_vf64_storeu(out + i, r);
    /* x - x is NaN exactly when x is infinite or NaN. */
    if (rsimd_mf64_any(rsimd_mf64_or(rsimd_vf64_is_nan(rsimd_vf64_sub(re, re)),
                                     rsimd_vf64_is_nan(rsimd_vf64_sub(im, im))))) {
      for (j = i; j < i + W; j++) {
        if (!isfinite(x[j].r) || !isfinite(x[j].i)) out[j] = rsimd_cmath_1(op, x[j]);
      }
    }
  }
  for (; i < n; i++) out[i] = rsimd_cmath_1(op, x[i]);
}

void RSIMD_KERNEL(math1_c128)(int op, const Rcomplex *x, R_xlen_t n, double *out);
void RSIMD_KERNEL(math1_c128)(int op, const Rcomplex *x, R_xlen_t n, double *out) {
  switch (op) {
  case RSIMD_CMATH_MOD: RSIMD_KERNEL(math1_c128_)(RSIMD_CMATH_MOD, 0, x, n, out); break;
  case RSIMD_CMATH_MOD | RSIMD_MATH_FAST:
    RSIMD_KERNEL(math1_c128_)(RSIMD_CMATH_MOD, 1, x, n, out);
    break;
  case RSIMD_CMATH_ARG: RSIMD_KERNEL(math1_c128_)(RSIMD_CMATH_ARG, 0, x, n, out); break;
  default: RSIMD_KERNEL(math1_c128_)(RSIMD_CMATH_ARG, 1, x, n, out); break;
  }
}

#else /* 32-bit ARM, or a SIMD tier without SLEEF: filled from below */
#define RSIMD_SKIP_math1_c128 1
#endif
