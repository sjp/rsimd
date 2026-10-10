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

/* ---- Elementary functions (cmath1_c128, cmath2_c128) ---------------------- */

/* One element of cmath1_c128 as base R: NA in both parts for an NA part,
   else b's function. Returns 1 for a NaN part from an element without
   one. */
static inline int rsimd_cmath1_1(const rsimd_cmath_base *b, int op, const Rcomplex *x,
                                 Rcomplex *out) {
  if (rsimd_is_na_f64(x->r) || rsimd_is_na_f64(x->i)) {
    out->r = out->i = rsimd_na_real();
    return 0;
  }
  b->f1[op](x, out);
  return (isnan(out->r) || isnan(out->i)) && !isnan(x->r) && !isnan(x->i);
}

/* One element of cmath2_c128 as base R: LOGB and ATAN2 give NA in both
   parts when all four parts are NA, and report a NaN part from operands
   without one; POW calls b's function alone. */
static inline int rsimd_cmath2_1(const rsimd_cmath_base *b, int op, const Rcomplex *x,
                                 const Rcomplex *y, Rcomplex *out) {
  if (op == RSIMD_CM2_POW) {
    b->f2[op](x, y, out);
    return 0;
  }
  if (rsimd_is_na_f64(x->r) && rsimd_is_na_f64(x->i) && rsimd_is_na_f64(y->r) &&
      rsimd_is_na_f64(y->i)) {
    out->r = out->i = rsimd_na_real();
    return 0;
  }
  b->f2[op](x, y, out);
  return (isnan(out->r) || isnan(out->i)) &&
         !(isnan(x->r) || isnan(x->i) || isnan(y->r) || isnan(y->i));
}

/* Is y a whole real number with |y| <= 65536 (base R's repeated-squaring
   powers)? */
static inline int rsimd_cm_whole_(Rcomplex y) {
  return y.i == 0 && y.r == trunc(y.r) && fabs(y.r) <= 65536;
}

/* x^k by base R's R_cpow_n, with the multiply (in the variants vr and vi,
   constants where specialised) and divide of `a`, each exactly base R's
   operator. */
RSIMD_ALWAYS_INLINE Rcomplex rsimd_cm_ipow_1_(const int vr, const int vi,
                                              const rsimd_c128_arith *a, Rcomplex x, int k) {
  Rcomplex z, one;
  int m = k < 0 ? -k : k;
  one.r = 1.0;
  one.i = 0.0;
  if (m == 0) return one;
  if (m == 1) {
    z = x;
  } else {
    z = one;
    while (m > 0) {
      if (m & 1) z = rsimd_cmul_in_(vr, vi, a, z, x);
      if (m == 1) break;
      m >>= 1;
      x = rsimd_cmul_in_(vr, vi, a, x, x);
    }
  }
  return k < 0 ? rsimd_cdiv_1_(a, one, z) : z;
}

/* cmath2_c128 for elements [i, n) one by one (flags RSIMD_EW_SCALAR(k)). */
static inline int rsimd_cmath2_from(const rsimd_cmath_base *b, int op, const Rcomplex *x,
                                    const Rcomplex *y, R_xlen_t i, R_xlen_t n, int flags,
                                    Rcomplex *out) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  int nan = 0;
  for (; i < n; i++) nan |= rsimd_cmath2_1(b, op, x + i * sx, y + i * sy, out + i);
  return nan;
}

#if RSIMD_TIER_IS(none)

int RSIMD_KERNEL(cmath1_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out,
                              const rsimd_cmath_base *b);
int RSIMD_KERNEL(cmath1_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out,
                              const rsimd_cmath_base *b) {
  R_xlen_t i;
  int nan = 0;
  op &= ~RSIMD_MATH_FAST;
  for (i = 0; i < n; i++) nan |= rsimd_cmath1_1(b, op, x + i, out + i);
  return nan;
}

/* x^y for every element: whole powers of non-zero bases multiply inline
   in the variants vr and vi (base R's operators exactly); the rest is base
   R's power. */
RSIMD_ALWAYS_INLINE void rsimd_cm_pow_from_(const int vr, const int vi, const Rcomplex *x,
                                            const Rcomplex *y, R_xlen_t n, int flags,
                                            Rcomplex *out, const rsimd_cmath_base *b,
                                            const rsimd_c128_arith *a) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    const Rcomplex *xi = x + i * sx, *yi = y + i * sy;
    if (rsimd_cm_whole_(*yi) && !(xi->r == 0 && xi->i == 0)) {
      out[i] = rsimd_cm_ipow_1_(vr, vi, a, *xi, (int) yi->r);
    } else {
      b->f2[RSIMD_CM2_POW](xi, yi, out + i);
    }
  }
}

int RSIMD_KERNEL(cmath2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                              Rcomplex *out, const rsimd_cmath_base *b,
                              const rsimd_c128_arith *a);
int RSIMD_KERNEL(cmath2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                              Rcomplex *out, const rsimd_cmath_base *b,
                              const rsimd_c128_arith *a) {
  op &= ~RSIMD_MATH_FAST;
  if (op != RSIMD_CM2_POW) return rsimd_cmath2_from(b, op, x, y, 0, n, flags, out);
  /* The multiply variants of GCC builds are specialised, as in ew2_c128. */
  if (a->mul_re == RSIMD_CMUL_FMA1 && a->mul_im == RSIMD_CMUL_FMA1) {
    rsimd_cm_pow_from_(RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1, x, y, n, flags, out, b, a);
  } else if (a->mul_re == RSIMD_CMUL_UNFUSED && a->mul_im == RSIMD_CMUL_UNFUSED) {
    rsimd_cm_pow_from_(RSIMD_CMUL_UNFUSED, RSIMD_CMUL_UNFUSED, x, y, n, flags, out, b, a);
  } else {
    rsimd_cm_pow_from_(a->mul_re, a->mul_im, x, y, n, flags, out, b, a);
  }
  return 0;
}

#elif !defined(RSIMD_NO_F64_SIMD) && defined(RSIMD_HAVE_SLEEF)

/* The vector tiers compute finite elements with the formulas below, built
   from SLEEF's real functions (their fast variants with `fast`), and send
   every element with a non-finite part, a formula result with a NaN part
   or an input in a region the formula cannot handle to b. Each formula
   returns the mask of the lanes it cannot handle. */

#define RSIMD_CM_V(c) rsimd_vf64_set1(c)
#define RSIMD_CM_F1(f, a) (fast ? rsimd_sleef_##f##_fast(a) : rsimd_sleef_##f(a))
#define RSIMD_CM_F2(f, a, b) (fast ? rsimd_sleef_##f##_fast(a, b) : rsimd_sleef_##f(a, b))

RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_cm_copysign_(rsimd_vf64 a, rsimd_vf64 s) {
  const rsimd_vf64 m = RSIMD_CM_V(-0.0);
  return rsimd_vf64_or(rsimd_vf64_andnot(m, a), rsimd_vf64_and(m, s));
}

RSIMD_ALWAYS_INLINE void rsimd_cm_sincos_(const int fast, rsimd_vf64 a, rsimd_vf64 *s,
                                          rsimd_vf64 *c) {
  if (fast) rsimd_sleef_sincos_fast(a, s, c);
  else rsimd_sleef_sincos(a, s, c);
}

/* a + b = *s + *e exactly (Knuth's two-sum). */
RSIMD_ALWAYS_INLINE void rsimd_cm_two_sum_(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 *s,
                                           rsimd_vf64 *e) {
  rsimd_vf64 t = rsimd_vf64_add(a, b), bb = rsimd_vf64_sub(t, a);
  *s = t;
  *e = rsimd_vf64_add(rsimd_vf64_sub(a, rsimd_vf64_sub(t, bb)), rsimd_vf64_sub(b, bb));
}

/* x^2 + y^2 - 1 with the squares and the sum carried in double-double, so
   that it is accurate relative to its value where it cancels (|z| near
   1). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_cm_x2y2m1_(rsimd_vf64 x, rsimd_vf64 y) {
  const rsimd_vf64 p1 = rsimd_vf64_mul(x, x), p2 = rsimd_vf64_mul(y, y);
  const rsimd_vf64 e1 = rsimd_vf64_fma(x, x, rsimd_vf64_neg(p1)),
                   e2 = rsimd_vf64_fma(y, y, rsimd_vf64_neg(p2));
  rsimd_vf64 s1, r1, s2, r2;
  rsimd_cm_two_sum_(rsimd_vf64_max(p1, p2), RSIMD_CM_V(-1.0), &s1, &r1);
  rsimd_cm_two_sum_(s1, rsimd_vf64_min(p1, p2), &s2, &r2);
  return rsimd_vf64_add(s2, rsimd_vf64_add(rsimd_vf64_add(r1, r2), rsimd_vf64_add(e1, e2)));
}

RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_is_zero_(rsimd_vf64 x, rsimd_vf64 y) {
  return rsimd_mf64_and(rsimd_vf64_cmp_eq(x, rsimd_vf64_zero()),
                        rsimd_vf64_cmp_eq(y, rsimd_vf64_zero()));
}

RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_cm_maxabs_(rsimd_vf64 x, rsimd_vf64 y) {
  return rsimd_vf64_max(rsimd_vf64_abs(x), rsimd_vf64_abs(y));
}

/* exp(x + iy) = e^x (cos y + i sin y); |x| > 708 is left out (overflow and
   subnormal results). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_exp_(const int fast, rsimd_vf64 x, rsimd_vf64 y,
                                             rsimd_vf64 *re, rsimd_vf64 *im) {
  rsimd_vf64 e = rsimd_sleef_exp(x), s, c;
  rsimd_cm_sincos_(fast, y, &s, &c);
  *re = rsimd_vf64_mul(e, c);
  *im = rsimd_vf64_mul(e, s);
  return rsimd_vf64_cmp_gt(rsimd_vf64_abs(x), RSIMD_CM_V(708.0));
}

/* log(x + iy) = log|z| + i atan2(y, x), with log|z| as log1p(x^2 + y^2 -
   1) / 2 for 1/4 <= |z|^2 <= 4, where log(hypot(x, y)) loses accuracy: the
   error of hypot is magnified by 1 / |log|z||, so the fast 3.5-ULP hypot
   alone would give up to 6 ULP between |z| = 1/2 and 2.
   Parts both below 2^-1000 (z = 0 included) or either beyond 2^1022
   (where |z| can overflow) are left out. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_log_(const int fast, rsimd_vf64 x, rsimd_vf64 y,
                                             rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 p = rsimd_vf64_add(rsimd_vf64_mul(x, x), rsimd_vf64_mul(y, y));
  const rsimd_mf64 near = rsimd_mf64_and(rsimd_vf64_cmp_ge(p, RSIMD_CM_V(0.25)),
                                         rsimd_vf64_cmp_le(p, RSIMD_CM_V(4.0)));
  const rsimd_vf64 m = rsimd_cm_maxabs_(x, y);
  rsimd_vf64 r = rsimd_vf64_zero();
  if (!rsimd_mf64_all(near)) r = RSIMD_CM_F1(log, RSIMD_CM_F2(hypot, x, y));
  if (rsimd_mf64_any(near)) {
    rsimd_vf64 l = rsimd_vf64_mul(RSIMD_CM_V(0.5), rsimd_sleef_log1p(rsimd_cm_x2y2m1_(x, y)));
    r = rsimd_vf64_blend(r, l, near);
  }
  *re = r;
  *im = RSIMD_CM_F2(atan2, y, x);
  /* SLEEF's hypot loses accuracy on subnormal parts, and overflows where
     |z| > DBL_MAX. */
  return rsimd_mf64_or(rsimd_vf64_cmp_lt(m, RSIMD_CM_V(0x1p-1000)),
                       rsimd_vf64_cmp_gt(m, RSIMD_CM_V(0x1p1022)));
}

/* sqrt(x + iy) with t = sqrt((|x| + |z|) / 2): t + i y / 2t for x >= 0,
   |y| / 2t + i copysign(t, y) for x < 0. Parts beyond 2^1020 (where |x| +
   |z| overflows) or both below 2^-1020 (where the halving loses bits) are
   left out; z = 0 gives NaN and is redone. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_sqrt_(const int fast, rsimd_vf64 x, rsimd_vf64 y,
                                              rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 m = rsimd_cm_maxabs_(x, y), ax = rsimd_vf64_abs(x);
  const rsimd_vf64 t = rsimd_vf64_sqrt(
    rsimd_vf64_mul(RSIMD_CM_V(0.5), rsimd_vf64_add(ax, RSIMD_CM_F2(hypot, x, y))));
  const rsimd_vf64 t2 = rsimd_vf64_add(t, t);
  const rsimd_mf64 neg = rsimd_vf64_cmp_lt(x, rsimd_vf64_zero());
  *re = rsimd_vf64_blend(t, rsimd_vf64_div(rsimd_vf64_abs(y), t2), neg);
  *im = rsimd_vf64_blend(rsimd_vf64_div(y, t2), rsimd_cm_copysign_(t, y), neg);
  return rsimd_mf64_or(rsimd_vf64_cmp_gt(m, RSIMD_CM_V(0x1p1020)),
                       rsimd_vf64_cmp_lt(m, RSIMD_CM_V(0x1p-1020)));
}

/* sinh(h) and cosh(h) from one expm1: with e = expm1(|h|) and E = e + 1,
   sinh = (e + e / E) / 2 (no cancellation for small |h|) and cosh = (E +
   1 / E) / 2. SLEEF has no faster expm1, so fast mode is the same. */
RSIMD_ALWAYS_INLINE void rsimd_cm_sinhcosh_(rsimd_vf64 h, rsimd_vf64 *sh, rsimd_vf64 *ch) {
  const rsimd_vf64 half = RSIMD_CM_V(0.5), e = rsimd_sleef_expm1(rsimd_vf64_abs(h));
  const rsimd_vf64 E = rsimd_vf64_add(e, RSIMD_CM_V(1.0));
  *sh = rsimd_cm_copysign_(rsimd_vf64_mul(half, rsimd_vf64_add(e, rsimd_vf64_div(e, E))), h);
  *ch = rsimd_vf64_mul(half, rsimd_vf64_add(E, rsimd_vf64_div(RSIMD_CM_V(1.0), E)));
}

/* sin, cos, sinh and cosh from the sine and cosine of one part and the
   hyperbolic sine and cosine of the other (h, left out beyond 709). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_trig_(const int op, const int fast, rsimd_vf64 x,
                                              rsimd_vf64 y, rsimd_vf64 *re, rsimd_vf64 *im) {
  const int hyp = op == RSIMD_CM_SINH || op == RSIMD_CM_COSH;
  const rsimd_vf64 t = hyp ? y : x, h = hyp ? x : y;
  rsimd_vf64 s, c, sh, ch;
  rsimd_cm_sinhcosh_(h, &sh, &ch);
  rsimd_cm_sincos_(fast, t, &s, &c);
  switch (op) {
  case RSIMD_CM_SIN: /* sin x cosh y + i cos x sinh y */
    *re = rsimd_vf64_mul(s, ch);
    *im = rsimd_vf64_mul(c, sh);
    break;
  case RSIMD_CM_COS: /* cos x cosh y - i sin x sinh y */
    *re = rsimd_vf64_mul(c, ch);
    *im = rsimd_vf64_neg(rsimd_vf64_mul(s, sh));
    break;
  case RSIMD_CM_SINH: /* sinh x cos y + i cosh x sin y */
    *re = rsimd_vf64_mul(sh, c);
    *im = rsimd_vf64_mul(ch, s);
    break;
  default: /* cosh x cos y + i sinh x sin y */
    *re = rsimd_vf64_mul(ch, c);
    *im = rsimd_vf64_mul(sh, s);
    break;
  }
  return rsimd_vf64_cmp_gt(rsimd_vf64_abs(h), RSIMD_CM_V(709.0));
}

/* tanh(x + iy) by Kahan's formula: with t = tan y, s = sinh x, b = 1 + t^2
   and d = 1 + b s^2, (b sqrt(1 + s^2) s + i t) / d. tan(z) is -i tanh(iz).
   |x| > 20 (the hyperbolic part) is left out. Below 2^-1000 t and s are y
   and x themselves, as SLEEF loses the sign of the smallest subnormals on
   some tiers. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_tanh_(const int tan, const int fast, rsimd_vf64 x,
                                              rsimd_vf64 y, rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 h = tan ? y : x, a = tan ? x : y;
  const rsimd_vf64 t = rsimd_math_tiny(RSIMD_CM_F1(tan, a), a);
  const rsimd_vf64 s = rsimd_math_tiny(RSIMD_CM_F1(sinh, h), h);
  const rsimd_vf64 one = RSIMD_CM_V(1.0), ss = rsimd_vf64_mul(s, s);
  const rsimd_vf64 b = rsimd_vf64_fma(t, t, one), d = rsimd_vf64_fma(b, ss, one);
  const rsimd_vf64 p = rsimd_vf64_div(rsimd_vf64_mul(rsimd_vf64_mul(b, rsimd_vf64_sqrt(rsimd_vf64_add(one, ss))), s), d);
  const rsimd_vf64 q = rsimd_vf64_div(t, d);
  /* tanh: p + iq; tan(x + iy) = -i tanh(-y + ix): q + ip with s = sinh y. */
  *re = tan ? q : p;
  *im = tan ? p : q;
  return rsimd_vf64_cmp_gt(rsimd_vf64_abs(h), RSIMD_CM_V(20.0));
}

/* The real part of asin (acos = 0) or acos (acos = 1) of |x| + i|y| and the
   magnitude of the imaginary part, by Hull, Fairgrieve and Tang's
   algorithm (ACM TOMS 23, 1997) without its regions for very large and
   very small parts: parts beyond 2^500 and |y| below 2^-500 but not zero
   are left out. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_asin_core_(const int acos, const int fast, rsimd_vf64 ax,
                                                   rsimd_vf64 ay, rsimd_vf64 *rr, rsimd_vf64 *ri) {
  const rsimd_vf64 one = RSIMD_CM_V(1.0), half = RSIMD_CM_V(0.5), zero = rsimd_vf64_zero();
  const rsimd_vf64 xp1 = rsimd_vf64_add(ax, one), xm1 = rsimd_vf64_sub(ax, one),
                   onemx = rsimd_vf64_sub(one, ax);
  /* Within the parts handled, the squares neither overflow nor lose bits
     to underflow, so sqrt(a^2 + b^2) replaces hypot. */
  const rsimd_vf64 yy = rsimd_vf64_mul(ay, ay);
  const rsimd_vf64 r = rsimd_vf64_sqrt(rsimd_vf64_fma(xp1, xp1, yy)),
                   s = rsimd_vf64_sqrt(rsimd_vf64_fma(xm1, xm1, yy));
  const rsimd_vf64 A = rsimd_vf64_mul(half, rsimd_vf64_add(r, s)), B = rsimd_vf64_div(ax, A);
  const rsimd_vf64 rx = rsimd_vf64_add(r, xp1);
  const rsimd_mf64 xle1 = rsimd_vf64_cmp_le(ax, one), small_b = rsimd_vf64_cmp_le(B, RSIMD_CM_V(0.6417));
  const rsimd_mf64 small_a = rsimd_vf64_cmp_le(A, RSIMD_CM_V(1.5));
  rsimd_vf64 re = zero, im = zero;
  /* Real part. */
  if (rsimd_mf64_any(small_b)) re = acos ? RSIMD_CM_F1(acos, B) : RSIMD_CM_F1(asin, B);
  if (!rsimd_mf64_all(small_b)) {
    const rsimd_vf64 apx = rsimd_vf64_add(A, ax);
    /* |x| <= 1: sqrt((A + x) / 2 (y^2 / (r + x + 1) + s + 1 - x));
       |x| > 1: y sqrt((A + x) / 2 (1 / (r + x + 1) + 1 / (s + x - 1))). */
    const rsimd_vf64 d1 = rsimd_vf64_sqrt(rsimd_vf64_mul(
      rsimd_vf64_mul(half, apx), rsimd_vf64_add(rsimd_vf64_div(yy, rx), rsimd_vf64_add(s, onemx))));
    const rsimd_vf64 d2 = rsimd_vf64_mul(
      ay, rsimd_vf64_sqrt(rsimd_vf64_mul(half, rsimd_vf64_add(rsimd_vf64_div(apx, rx),
                                                              rsimd_vf64_div(apx, rsimd_vf64_add(s, xm1))))));
    const rsimd_vf64 d = rsimd_vf64_blend(d2, d1, xle1);
    const rsimd_vf64 v = acos ? RSIMD_CM_F1(atan, rsimd_vf64_div(d, ax))
                              : RSIMD_CM_F1(atan, rsimd_vf64_div(ax, d));
    re = rsimd_vf64_blend(v, re, small_b);
  }
  /* Imaginary part. */
  if (rsimd_mf64_any(small_a)) {
    /* A - 1 = (y^2 / (r + x + 1) + y^2 / (s + 1 - x)) / 2 for |x| < 1, else
       (y^2 / (r + x + 1) + s + x - 1) / 2. */
    const rsimd_vf64 t1 = rsimd_vf64_div(yy, rx);
    const rsimd_vf64 am1 = rsimd_vf64_mul(
      half, rsimd_vf64_add(t1, rsimd_vf64_blend(rsimd_vf64_add(s, xm1),
                                                rsimd_vf64_div(yy, rsimd_vf64_add(s, onemx)),
                                                rsimd_vf64_cmp_lt(ax, one))));
    im = rsimd_sleef_log1p(
      rsimd_vf64_add(am1, rsimd_vf64_sqrt(rsimd_vf64_mul(am1, rsimd_vf64_add(A, one)))));
  }
  if (!rsimd_mf64_all(small_a)) {
    const rsimd_vf64 v = RSIMD_CM_F1(
      log, rsimd_vf64_add(A, rsimd_vf64_sqrt(rsimd_vf64_fma(A, A, RSIMD_CM_V(-1.0)))));
    im = rsimd_vf64_blend(v, im, small_a);
  }
  *rr = re;
  *ri = im;
  return rsimd_mf64_or(
    rsimd_vf64_cmp_gt(rsimd_vf64_max(ax, ay), RSIMD_CM_V(0x1p500)),
    rsimd_mf64_and(rsimd_vf64_cmp_gt(ay, zero), rsimd_vf64_cmp_lt(ay, RSIMD_CM_V(0x1p-500))));
}

/* The lanes on base R's branch cut for asin and acos: b = 0 with |a| > 1
   (and for atan with the parts the other way round). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_on_cut_(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_mf64_and(rsimd_vf64_cmp_eq(b, rsimd_vf64_zero()),
                        rsimd_vf64_cmp_gt(rsimd_vf64_abs(a), RSIMD_CM_V(1.0)));
}

/* Base R's z_asin on its branch cut (y = 0, |x| > 1): with t1 = |x + 1| / 2,
   t2 = |x - 1| / 2 and a = t1 + t2, asin(t1 - t2) + i log(a + sqrt(a^2 -
   1)), the imaginary part negated for x > 1, the operations as base R's
   (only asin and log are SLEEF's). t1 - t2 is +-1 unless x + 1 or x - 1
   rounds; beyond +-1 it gives base R's NaN, which the caller redoes. a^2 -
   1 is rounded as base R's code rounds it (`how`, an RSIMD_CUT_* value).
   Returns the lanes it leaves to base R: where `how` is unknown, a < 1.5,
   where a^2 - 1 cancels and the two roundings give different results. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_asin_cut_(const int fast, const int how, rsimd_vf64 x,
                                                  rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 one = RSIMD_CM_V(1.0), half = RSIMD_CM_V(0.5);
  const rsimd_vf64 t1 = rsimd_vf64_mul(half, rsimd_vf64_abs(rsimd_vf64_add(x, one))),
                   t2 = rsimd_vf64_mul(half, rsimd_vf64_abs(rsimd_vf64_sub(x, one)));
  const rsimd_vf64 a = rsimd_vf64_add(t1, t2), d = rsimd_vf64_sub(t1, t2);
  const rsimd_vf64 a21 = how == RSIMD_CUT_FUSED ? rsimd_vf64_fma(a, a, rsimd_vf64_neg(one))
                                                : rsimd_vf64_sub(rsimd_vf64_mul(a, a), one);
  const rsimd_vf64 l = RSIMD_CM_F1(log, rsimd_vf64_add(a, rsimd_vf64_sqrt(a21)));
  const rsimd_mf64 unit = rsimd_vf64_cmp_eq(rsimd_vf64_abs(d), one);
  /* asin(+-1) is +-pi/2 correctly rounded, as the C library gives it. */
  rsimd_vf64 as = rsimd_vf64_mul(d, RSIMD_CM_V(M_PI_2));
  if (!rsimd_mf64_all(unit)) as = rsimd_vf64_blend(RSIMD_CM_F1(asin, d), as, unit);
  *im = rsimd_vf64_blend(l, rsimd_vf64_neg(l), rsimd_vf64_cmp_gt(x, one));
  /* Base R adds the real part of im * I, a zero, or NaN where im is
     infinite (a^2 overflows). */
  *re = rsimd_vf64_add(as, rsimd_vf64_mul(*im, rsimd_vf64_zero()));
  if (how != RSIMD_CUT_UNKNOWN) return rsimd_mf64_none();
  return rsimd_vf64_cmp_lt(a, RSIMD_CM_V(1.5));
}

/* asin (acos = 0) or acos (acos = 1) of x + iy: base R's branch-cut code
   in the lanes on its cut (base R's acos there is pi/2 - asin), Hull et
   al. elsewhere, each computed only if some lane needs it; `how` as in
   rsimd_cm_asin_cut_(). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_asin_(const int acos, const int fast, const int how,
                                              rsimd_vf64 x, rsimd_vf64 y, rsimd_vf64 *re,
                                              rsimd_vf64 *im) {
  const rsimd_mf64 cut = rsimd_cm_on_cut_(x, y);
  rsimd_vf64 r = rsimd_vf64_zero(), i = r;
  rsimd_mf64 m = rsimd_mf64_none();
  if (!rsimd_mf64_all(cut)) {
    m = rsimd_cm_asin_core_(acos, fast, rsimd_vf64_abs(x), rsimd_vf64_abs(y), &r, &i);
    if (acos) { /* (R or pi - R for x < 0, -copysign(I, y)) */
      r = rsimd_vf64_blend(r, rsimd_vf64_sub(RSIMD_CM_V(M_PI), r),
                           rsimd_vf64_cmp_lt(x, rsimd_vf64_zero()));
      i = rsimd_vf64_neg(rsimd_cm_copysign_(i, y));
    } else { /* (copysign(R, x), copysign(I, y)) */
      r = rsimd_cm_copysign_(r, x);
      i = rsimd_cm_copysign_(i, y);
    }
  }
  if (rsimd_mf64_any(cut)) {
    rsimd_vf64 cr, ci;
    const rsimd_mf64 cm = rsimd_cm_asin_cut_(fast, how, x, &cr, &ci);
    if (acos) {
      cr = rsimd_vf64_sub(RSIMD_CM_V(M_PI_2), cr);
      ci = rsimd_vf64_neg(ci);
    }
    r = rsimd_vf64_blend(r, cr, cut);
    i = rsimd_vf64_blend(i, ci, cut);
    m = rsimd_mf64_or(rsimd_mf64_andnot(cut, m), rsimd_mf64_and(cut, cm));
  }
  *re = r;
  *im = i;
  return m;
}

/* atan(x + iy) as C99's catan: (atan2(2x, 1 - x^2 - y^2) + i log1p(4|y| /
   (x^2 + (1 - |y|)^2)) sign(y)) / (2, 4), with 1 - x^2 - y^2 in
   double-double. Parts beyond 2^500 are left out, and so are x = 0 with
   |y| >= 1 (the poles, and base R's branch cut, which rsimd_cm_atan_()
   computes) and x^2 + (1 - |y|)^2 below 2^-960 (near the poles, where it
   underflows). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_atan_core_(const int fast, rsimd_vf64 x, rsimd_vf64 y,
                                                   rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 ay = rsimd_vf64_abs(y), one = RSIMD_CM_V(1.0), oy = rsimd_vf64_sub(one, ay);
  const rsimd_vf64 d = rsimd_vf64_neg(rsimd_cm_x2y2m1_(x, y));
  const rsimd_vf64 den = rsimd_vf64_fma(x, x, rsimd_vf64_mul(oy, oy));
  const rsimd_vf64 q = rsimd_vf64_div(rsimd_vf64_mul(RSIMD_CM_V(4.0), ay), den);
  *re = rsimd_vf64_mul(RSIMD_CM_V(0.5), RSIMD_CM_F2(atan2, rsimd_vf64_add(x, x), d));
  *im = rsimd_cm_copysign_(rsimd_vf64_mul(RSIMD_CM_V(0.25), rsimd_sleef_log1p(q)), y);
  return rsimd_mf64_or(rsimd_mf64_or(rsimd_vf64_cmp_gt(rsimd_cm_maxabs_(x, y), RSIMD_CM_V(0x1p500)),
                                     rsimd_vf64_cmp_lt(den, RSIMD_CM_V(0x1p-960))),
                       rsimd_mf64_and(rsimd_vf64_cmp_eq(x, rsimd_vf64_zero()),
                                      rsimd_vf64_cmp_ge(ay, one)));
}

/* atan(x + iy): base R's branch-cut code in the lanes on its cut (x = 0,
   |y| > 1), +-pi/2 + i log((y + 1)^2 / (y - 1)^2) / 4 with base R's
   operations (only the logarithm is SLEEF's), C99's formula elsewhere,
   each computed only if some lane needs it. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_atan_(const int fast, rsimd_vf64 x, rsimd_vf64 y,
                                              rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_mf64 cut = rsimd_cm_on_cut_(y, x);
  rsimd_vf64 r = rsimd_vf64_zero(), i = r;
  rsimd_mf64 m = rsimd_mf64_none();
  if (!rsimd_mf64_all(cut)) m = rsimd_cm_atan_core_(fast, x, y, &r, &i);
  if (rsimd_mf64_any(cut)) {
    const rsimd_vf64 one = RSIMD_CM_V(1.0);
    const rsimd_vf64 yp = rsimd_vf64_add(y, one), ym = rsimd_vf64_sub(y, one);
    const rsimd_vf64 q = rsimd_vf64_div(rsimd_vf64_mul(yp, yp), rsimd_vf64_mul(ym, ym));
    r = rsimd_vf64_blend(r, rsimd_cm_copysign_(RSIMD_CM_V(M_PI_2), y), cut);
    i = rsimd_vf64_blend(i, rsimd_vf64_mul(RSIMD_CM_V(0.25), RSIMD_CM_F1(log, q)), cut);
    m = rsimd_mf64_andnot(cut, m);
  }
  *re = r;
  *im = i;
  return m;
}

/* z * i and -i * z as base R's z_asinh and z_atanh compute them: full
   complex products with 0 + 1i and -0 - 1i (exact, but they decide the
   signs of zero parts). */
RSIMD_ALWAYS_INLINE void rsimd_cm_mul_i_(rsimd_vf64 x, rsimd_vf64 y, rsimd_vf64 *re,
                                         rsimd_vf64 *im) {
  const rsimd_vf64 zero = rsimd_vf64_zero();
  *re = rsimd_vf64_sub(rsimd_vf64_mul(x, zero), y);
  *im = rsimd_vf64_add(x, rsimd_vf64_mul(y, zero));
}
RSIMD_ALWAYS_INLINE void rsimd_cm_mul_negi_(rsimd_vf64 x, rsimd_vf64 y, rsimd_vf64 *re,
                                            rsimd_vf64 *im) {
  const rsimd_vf64 nz = RSIMD_CM_V(-0.0);
  *re = rsimd_vf64_sub(rsimd_vf64_mul(nz, x), rsimd_vf64_neg(y));
  *im = rsimd_vf64_add(rsimd_vf64_mul(nz, y), rsimd_vf64_neg(x));
}

/* Function op of x + iy, finite; `how` is b->asin_cut. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm1_v_(const int op, const int fast, const int how,
                                            rsimd_vf64 x, rsimd_vf64 y, rsimd_vf64 *re,
                                            rsimd_vf64 *im) {
  rsimd_vf64 r, i;
  rsimd_mf64 m;
  switch (op) {
  case RSIMD_CM_SQRT: return rsimd_cm_sqrt_(fast, x, y, re, im);
  case RSIMD_CM_EXP: return rsimd_cm_exp_(fast, x, y, re, im);
  case RSIMD_CM_LOG: return rsimd_cm_log_(fast, x, y, re, im);
  case RSIMD_CM_SIN:
  case RSIMD_CM_COS:
  case RSIMD_CM_SINH:
  case RSIMD_CM_COSH: return rsimd_cm_trig_(op, fast, x, y, re, im);
  case RSIMD_CM_TAN: return rsimd_cm_tanh_(1, fast, x, y, re, im);
  case RSIMD_CM_TANH: return rsimd_cm_tanh_(0, fast, x, y, re, im);
  case RSIMD_CM_ASIN: return rsimd_cm_asin_(0, fast, how, x, y, re, im);
  case RSIMD_CM_ACOS: return rsimd_cm_asin_(1, fast, how, x, y, re, im);
  case RSIMD_CM_ACOSH: /* C99's cacosh: (I, copysign(acos's real part, y)) */
    m = rsimd_cm_asin_core_(1, fast, rsimd_vf64_abs(x), rsimd_vf64_abs(y), &r, &i);
    r = rsimd_vf64_blend(r, rsimd_vf64_sub(RSIMD_CM_V(M_PI), r),
                         rsimd_vf64_cmp_lt(x, rsimd_vf64_zero()));
    *re = i;
    *im = rsimd_cm_copysign_(r, y);
    return m;
  case RSIMD_CM_ASINH: /* -i asin(iz), multiplied as base R does */
    rsimd_cm_mul_i_(x, y, &x, &y);
    m = rsimd_cm_asin_(0, fast, how, x, y, &r, &i);
    rsimd_cm_mul_negi_(r, i, re, im);
    return m;
  case RSIMD_CM_ATAN: return rsimd_cm_atan_(fast, x, y, re, im);
  default: /* ATANH: -i atan(iz) */
    rsimd_cm_mul_i_(x, y, &x, &y);
    m = rsimd_cm_atan_(fast, x, y, &r, &i);
    rsimd_cm_mul_negi_(r, i, re, im);
    return m;
  }
}

/* Loads the W complex numbers from x[i] (x[0] for all with `bcast`), or
   the m < W left from it padded with 1 + 0i. */
RSIMD_ALWAYS_INLINE void rsimd_cm_load_(const Rcomplex *x, ptrdiff_t i, ptrdiff_t m, int bcast,
                                        rsimd_vf64 *re, rsimd_vf64 *im) {
  const ptrdiff_t W = RSIMD_LANES_64;
  if (bcast) {
    *re = rsimd_vf64_set1(x[0].r);
    *im = rsimd_vf64_set1(x[0].i);
  } else if (m == W) {
    rsimd_c128_load_(x, i, re, im);
  } else {
    Rcomplex t[RSIMD_MAX_LANES_64];
    ptrdiff_t j;
    for (j = 0; j < W; j++) {
      t[j].r = j < m ? x[i + j].r : 1.0;
      t[j].i = j < m ? x[i + j].i : 0.0;
    }
    rsimd_c128_load_(t, 0, re, im);
  }
}

/* Stores m <= W complex numbers to out[i]. */
RSIMD_ALWAYS_INLINE void rsimd_cm_store_(Rcomplex *out, ptrdiff_t i, ptrdiff_t m, rsimd_vf64 re,
                                         rsimd_vf64 im) {
  if (m == RSIMD_LANES_64) {
    rsimd_c128_store_(out, i, re, im);
  } else {
    Rcomplex t[RSIMD_MAX_LANES_64];
    rsimd_c128_store_(t, 0, re, im);
    memcpy(out + i, t, (size_t) m * sizeof(Rcomplex));
  }
}

/* The lanes of m as doubles in d (non-zero for a set lane); 1 if any. */
RSIMD_ALWAYS_INLINE int rsimd_cm_lanes_(rsimd_mf64 m, double *d) {
  if (!rsimd_mf64_any(m)) return 0;
  rsimd_vf64_storeu(d, rsimd_vf64_blend(rsimd_vf64_zero(), RSIMD_CM_V(1.0), m));
  return 1;
}

RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm_nonfinite_(rsimd_vf64 x) {
  return rsimd_vf64_is_nan(rsimd_vf64_sub(x, x));
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(cmath1_)(const int op, const int fast, const Rcomplex *x,
                                              R_xlen_t n, Rcomplex *out,
                                              const rsimd_cmath_base *b) {
  const ptrdiff_t W = RSIMD_LANES_64;
  double lanes[RSIMD_MAX_LANES_64];
  ptrdiff_t i, j;
  int nan = 0;
  for (i = 0; i < n; i += W) {
    const ptrdiff_t m = n - i < W ? (ptrdiff_t) n - i : W;
    rsimd_vf64 xr, xi, re, im;
    rsimd_mf64 redo;
    rsimd_cm_load_(x, i, m, 0, &xr, &xi);
    redo = rsimd_cm1_v_(op, fast, b->asin_cut, xr, xi, &re, &im);
    redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_cm_nonfinite_(xr), rsimd_cm_nonfinite_(xi)));
    redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_vf64_is_nan(re), rsimd_vf64_is_nan(im)));
    /* On an axis, exp and the (hyperbolic) sines and cosines give a zero
       part whose sign C libraries choose differently (macOS's libm from
       glibc), so those elements are base R's. */
    if (op == RSIMD_CM_EXP || op == RSIMD_CM_SIN || op == RSIMD_CM_COS || op == RSIMD_CM_SINH ||
        op == RSIMD_CM_COSH) {
      redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_vf64_cmp_eq(xr, rsimd_vf64_zero()),
                                               rsimd_vf64_cmp_eq(xi, rsimd_vf64_zero())));
    }
    rsimd_cm_store_(out, i, m, re, im);
    if (rsimd_cm_lanes_(redo, lanes)) {
      for (j = 0; j < m; j++) {
        if (lanes[j] != 0) nan |= rsimd_cmath1_1(b, op, x + i + j, out + i + j);
      }
    }
  }
  return nan;
}

int RSIMD_KERNEL(cmath1_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out,
                              const rsimd_cmath_base *b);
int RSIMD_KERNEL(cmath1_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out,
                              const rsimd_cmath_base *b) {
  const int fast = (op & RSIMD_MATH_FAST) != 0;
#define RSIMD_CM_CASE(o)                                                                     \
  case o: return fast ? RSIMD_KERNEL(cmath1_)(o, 1, x, n, out, b)                          \
                      : RSIMD_KERNEL(cmath1_)(o, 0, x, n, out, b);
  switch (op & ~RSIMD_MATH_FAST) {
    RSIMD_CM_CASE(RSIMD_CM_SQRT)
    RSIMD_CM_CASE(RSIMD_CM_EXP)
    RSIMD_CM_CASE(RSIMD_CM_LOG)
    RSIMD_CM_CASE(RSIMD_CM_SIN)
    RSIMD_CM_CASE(RSIMD_CM_COS)
    RSIMD_CM_CASE(RSIMD_CM_TAN)
    RSIMD_CM_CASE(RSIMD_CM_SINH)
    RSIMD_CM_CASE(RSIMD_CM_COSH)
    RSIMD_CM_CASE(RSIMD_CM_TANH)
    RSIMD_CM_CASE(RSIMD_CM_ASIN)
    RSIMD_CM_CASE(RSIMD_CM_ACOS)
    RSIMD_CM_CASE(RSIMD_CM_ATAN)
    RSIMD_CM_CASE(RSIMD_CM_ASINH)
    RSIMD_CM_CASE(RSIMD_CM_ACOSH)
    RSIMD_CM_CASE(RSIMD_CM_ATANH)
  default: return 0;
  }
#undef RSIMD_CM_CASE
}

/* x^k for whole k (|k| <= 65536) as base R: binary powering by the
   multiply of variants vr and vi (z = 1, then z = z * x for each set bit
   and x = x * x), x itself for |k| = 1, 1 + 0i for 0, and the reciprocal
   1 / x^|k| by the division (`fused`) for k < 0. Lanes whose result has
   a NaN part are redone by the caller. */
RSIMD_ALWAYS_INLINE void rsimd_cm_ipow_(int vr, int vi, int fused, rsimd_vf64 xr, rsimd_vf64 xi,
                                        rsimd_vf64 k, rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 one = RSIMD_CM_V(1.0), zero = rsimd_vf64_zero();
  rsimd_vf64 K = rsimd_vf64_abs(k), zr = one, zi = zero, pr = xr, pi = xi;
  rsimd_mf64 act = rsimd_vf64_cmp_ge(K, one), m;
  while (rsimd_mf64_any(act)) {
    const rsimd_vf64 h = rsimd_vf64_floor(rsimd_vf64_mul(K, RSIMD_CM_V(0.5)));
    const rsimd_mf64 odd = rsimd_mf64_and(act, rsimd_vf64_cmp_ne(K, rsimd_vf64_add(h, h)));
    rsimd_vf64 tr, ti;
    if (rsimd_mf64_any(odd)) {
      rsimd_cmul_v_(vr, vi, zr, zi, pr, pi, &tr, &ti);
      zr = rsimd_vf64_blend(zr, tr, odd);
      zi = rsimd_vf64_blend(zi, ti, odd);
    }
    act = rsimd_mf64_and(act, rsimd_vf64_cmp_gt(K, one));
    K = h;
    if (rsimd_mf64_any(act)) {
      rsimd_cmul_v_(vr, vi, pr, pi, pr, pi, &tr, &ti);
      pr = tr;
      pi = ti;
    }
  }
  m = rsimd_vf64_cmp_eq(rsimd_vf64_abs(k), one);
  zr = rsimd_vf64_blend(zr, xr, m);
  zi = rsimd_vf64_blend(zi, xi, m);
  m = rsimd_vf64_cmp_lt(k, zero);
  if (rsimd_mf64_any(m)) {
    rsimd_vf64 qr, qi;
    rsimd_cdiv_v_(fused, one, zero, zr, zi, &qr, &qi);
    zr = rsimd_vf64_blend(zr, qr, m);
    zi = rsimd_vf64_blend(zi, qi, m);
  }
  *re = zr;
  *im = zi;
}

/* op (POW, LOGB or ATAN2) of finite x and y. */
/* (a + bi) / (c + di) as base R divides: the vector formula of the
   load-time variant, or base R's operator lane by lane where no variant
   matched (RSIMD_CDIV_SCALAR). */
RSIMD_ALWAYS_INLINE void rsimd_cm_div_(const rsimd_c128_arith *ar, rsimd_vf64 a, rsimd_vf64 b,
                                       rsimd_vf64 c, rsimd_vf64 d, rsimd_vf64 *re,
                                       rsimd_vf64 *im) {
  if (ar->div != RSIMD_CDIV_SCALAR) {
    rsimd_cdiv_v_(ar->div == RSIMD_CDIV_FMA, a, b, c, d, re, im);
  } else {
    Rcomplex x[RSIMD_MAX_LANES_64], y[RSIMD_MAX_LANES_64], q[RSIMD_MAX_LANES_64];
    ptrdiff_t j;
    rsimd_c128_store_(x, 0, a, b);
    rsimd_c128_store_(y, 0, c, d);
    for (j = 0; j < RSIMD_LANES_64; j++) ar->div1(x + j, y + j, q + j);
    rsimd_c128_load_(q, 0, re, im);
  }
}

/* `ylog` (LOGB with a broadcast base): yr + i yi is already log(base). */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_cm2_v_(const int op, const int fast, const int ylog,
                                            const rsimd_c128_arith *a, rsimd_vf64 xr,
                                            rsimd_vf64 xi, rsimd_vf64 yr, rsimd_vf64 yi,
                                            rsimd_vf64 *re, rsimd_vf64 *im) {
  const int fused = a->div != RSIMD_CDIV_UNFUSED;
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_vf64 lr, li, mr, mi;
  rsimd_mf64 redo;
  if (op == RSIMD_CM2_LOGB) {
    redo = rsimd_cm_log_(fast, xr, xi, &lr, &li);
    if (ylog) {
      mr = yr;
      mi = yi;
    } else {
      redo = rsimd_mf64_or(redo, rsimd_cm_log_(fast, yr, yi, &mr, &mi));
    }
    rsimd_cm_div_(a, lr, li, mr, mi, re, im);
    return redo;
  }
  if (op == RSIMD_CM2_ATAN2) {
    /* catan(x / y), plus pi for Re y < 0, minus 2 pi above pi; y = 0 is
       base R's special case. */
    rsimd_cm_div_(a, xr, xi, yr, yi, &lr, &li);
    redo = rsimd_mf64_or(rsimd_cm_atan_core_(fast, lr, li, &mr, &mi), rsimd_cm_is_zero_(yr, yi));
    mr = rsimd_vf64_blend(mr, rsimd_vf64_add(mr, RSIMD_CM_V(M_PI)), rsimd_vf64_cmp_lt(yr, zero));
    mr = rsimd_vf64_blend(mr, rsimd_vf64_sub(mr, RSIMD_CM_V(2 * M_PI)),
                          rsimd_vf64_cmp_gt(mr, RSIMD_CM_V(M_PI)));
    *re = mr;
    *im = mi;
    return redo;
  }
  /* POW: 0^y is base R's special case; whole real y with |y| <= 65536
     powers by multiplication, other y use exp(y log x). */
  {
    const rsimd_mf64 whole = rsimd_mf64_and(
      rsimd_mf64_and(rsimd_vf64_cmp_eq(yi, zero), rsimd_vf64_cmp_eq(yr, rsimd_vf64_trunc(yr))),
      rsimd_vf64_cmp_le(rsimd_vf64_abs(yr), RSIMD_CM_V(65536.0)));
    rsimd_vf64 gr = zero, gi = zero, wr, wi;
    redo = rsimd_cm_is_zero_(xr, xi);
    if (!rsimd_mf64_all(whole)) {
      redo = rsimd_mf64_or(redo, rsimd_cm_log_(fast, xr, xi, &lr, &li));
      wr = rsimd_vf64_sub(rsimd_vf64_mul(yr, lr), rsimd_vf64_mul(yi, li));
      wi = rsimd_vf64_add(rsimd_vf64_mul(yr, li), rsimd_vf64_mul(yi, lr));
      redo = rsimd_mf64_or(redo, rsimd_mf64_and(rsimd_mf64_not(whole),
                                                rsimd_cm_exp_(fast, wr, wi, &gr, &gi)));
      redo = rsimd_mf64_or(redo, rsimd_mf64_and(rsimd_mf64_not(whole),
                                                rsimd_mf64_or(rsimd_cm_nonfinite_(wr),
                                                              rsimd_cm_nonfinite_(wi))));
    }
    if (rsimd_mf64_any(whole)) {
      if (a->mul_re == RSIMD_CMUL_SCALAR || a->mul_im == RSIMD_CMUL_SCALAR ||
          (a->div == RSIMD_CDIV_SCALAR && rsimd_mf64_any(rsimd_vf64_cmp_lt(yr, zero)))) {
        redo = rsimd_mf64_or(redo, whole);
      } else {
        rsimd_cm_ipow_(a->mul_re, a->mul_im, fused, xr, xi,
                       rsimd_vf64_blend(zero, yr, whole), &lr, &li);
        gr = rsimd_vf64_blend(gr, lr, whole);
        gi = rsimd_vf64_blend(gi, li, whole);
      }
    }
    *re = gr;
    *im = gi;
    return redo;
  }
}

/* x^k for one whole k (|k| <= 65536) and every x: the repeated squaring
   of rsimd_cm_ipow_1_() on whole vectors, in the multiply variant vr, vi
   (constants where specialised), for the elements of x with finite parts
   and a result without a NaN part; the others are redone by b. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cm_ipow_k_)(const int vr, const int vi, const int fused,
                                                  const Rcomplex *x, R_xlen_t n, int k,
                                                  Rcomplex *out, const rsimd_cmath_base *b) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int m0 = k < 0 ? -k : k;
  double lanes[RSIMD_MAX_LANES_64];
  ptrdiff_t i, j;
  Rcomplex y;
  y.r = (double) k;
  y.i = 0.0;
  for (i = 0; i < n; i += W) {
    const ptrdiff_t m = n - i < W ? (ptrdiff_t) n - i : W;
    rsimd_vf64 xr, xi, zr, zi, pr, pi, tr, ti;
    rsimd_mf64 redo;
    int mm = m0;
    rsimd_cm_load_(x, i, m, 0, &xr, &xi);
    if (mm == 0) {
      zr = RSIMD_CM_V(1.0);
      zi = rsimd_vf64_zero();
    } else if (mm == 1) {
      zr = xr;
      zi = xi;
    } else {
      zr = RSIMD_CM_V(1.0);
      zi = rsimd_vf64_zero();
      pr = xr;
      pi = xi;
      for (;;) {
        if (mm & 1) {
          rsimd_cmul_v_(vr, vi, zr, zi, pr, pi, &tr, &ti);
          zr = tr;
          zi = ti;
        }
        if (mm == 1) break;
        mm >>= 1;
        rsimd_cmul_v_(vr, vi, pr, pi, pr, pi, &tr, &ti);
        pr = tr;
        pi = ti;
      }
    }
    if (k < 0) {
      rsimd_cdiv_v_(fused, RSIMD_CM_V(1.0), rsimd_vf64_zero(), zr, zi, &tr, &ti);
      zr = tr;
      zi = ti;
    }
    redo = rsimd_mf64_or(rsimd_cm_is_zero_(xr, xi),
                         rsimd_mf64_or(rsimd_cm_nonfinite_(xr), rsimd_cm_nonfinite_(xi)));
    redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_vf64_is_nan(zr), rsimd_vf64_is_nan(zi)));
    rsimd_cm_store_(out, i, m, zr, zi);
    if (rsimd_cm_lanes_(redo, lanes)) {
      for (j = 0; j < m; j++) {
        if (lanes[j] != 0) b->f2[RSIMD_CM2_POW](x + i + j, &y, out + i + j);
      }
    }
  }
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(cmath2_)(const int op, const int fast, const Rcomplex *x,
                                              const Rcomplex *y, R_xlen_t n, int flags,
                                              Rcomplex *out, const rsimd_cmath_base *b,
                                              const rsimd_c128_arith *a) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  double lanes[RSIMD_MAX_LANES_64];
  ptrdiff_t i, j;
  int nan = 0, ylog = 0;
  rsimd_vf64 lyr = rsimd_vf64_zero(), lyi = rsimd_vf64_zero();
  if (op == RSIMD_CM2_POW && sy && !sx && rsimd_cm_whole_(y[0]) && a->mul_re != RSIMD_CMUL_SCALAR &&
      a->mul_im != RSIMD_CMUL_SCALAR && (y[0].r >= 0 || a->div != RSIMD_CDIV_SCALAR)) {
    /* One whole exponent: no masks, and the common variants specialised. */
    const int k = (int) y[0].r, fused = a->div != RSIMD_CDIV_UNFUSED;
    if (a->mul_re == RSIMD_CMUL_FMA1 && a->mul_im == RSIMD_CMUL_FMA1) {
      RSIMD_KERNEL(cm_ipow_k_)(RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1, fused, x, n, k, out, b);
    } else if (a->mul_re == RSIMD_CMUL_UNFUSED && a->mul_im == RSIMD_CMUL_UNFUSED) {
      RSIMD_KERNEL(cm_ipow_k_)(RSIMD_CMUL_UNFUSED, RSIMD_CMUL_UNFUSED, fused, x, n, k, out, b);
    } else {
      RSIMD_KERNEL(cm_ipow_k_)(a->mul_re, a->mul_im, fused, x, n, k, out, b);
    }
    return 0;
  }
  if (op == RSIMD_CM2_LOGB && sy && isfinite(y[0].r) && isfinite(y[0].i)) {
    /* A broadcast base: its log once. */
    rsimd_vf64 br = rsimd_vf64_set1(y[0].r), bi = rsimd_vf64_set1(y[0].i);
    double t[RSIMD_MAX_LANES_64];
    rsimd_mf64 bad = rsimd_cm_log_(fast, br, bi, &lyr, &lyi);
    bad = rsimd_mf64_or(bad, rsimd_mf64_or(rsimd_cm_nonfinite_(lyr), rsimd_cm_nonfinite_(lyi)));
    ylog = !rsimd_cm_lanes_(bad, t);
  }
  for (i = 0; i < n; i += W) {
    const ptrdiff_t m = n - i < W ? (ptrdiff_t) n - i : W;
    rsimd_vf64 xr, xi, yr, yi, re, im;
    rsimd_mf64 redo;
    rsimd_cm_load_(x, i, m, sx, &xr, &xi);
    if (ylog) {
      yr = lyr;
      yi = lyi;
    } else {
      rsimd_cm_load_(y, i, m, sy, &yr, &yi);
    }
    redo = rsimd_cm2_v_(op, fast, ylog, a, xr, xi, yr, yi, &re, &im);
    redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_cm_nonfinite_(xr), rsimd_cm_nonfinite_(xi)));
    if (!ylog) {
      redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_cm_nonfinite_(yr), rsimd_cm_nonfinite_(yi)));
    }
    redo = rsimd_mf64_or(redo, rsimd_mf64_or(rsimd_vf64_is_nan(re), rsimd_vf64_is_nan(im)));
    rsimd_cm_store_(out, i, m, re, im);
    if (rsimd_cm_lanes_(redo, lanes)) {
      for (j = 0; j < m; j++) {
        if (lanes[j] != 0) {
          nan |= rsimd_cmath2_1(b, op, x + (sx ? 0 : i + j), y + (sy ? 0 : i + j), out + i + j);
        }
      }
    }
  }
  return nan;
}

int RSIMD_KERNEL(cmath2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                              Rcomplex *out, const rsimd_cmath_base *b,
                              const rsimd_c128_arith *a);
int RSIMD_KERNEL(cmath2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                              Rcomplex *out, const rsimd_cmath_base *b,
                              const rsimd_c128_arith *a) {
  const int fast = (op & RSIMD_MATH_FAST) != 0;
  switch (op & ~RSIMD_MATH_FAST) {
  case RSIMD_CM2_POW:
    return fast ? RSIMD_KERNEL(cmath2_)(RSIMD_CM2_POW, 1, x, y, n, flags, out, b, a)
                : RSIMD_KERNEL(cmath2_)(RSIMD_CM2_POW, 0, x, y, n, flags, out, b, a);
  case RSIMD_CM2_LOGB:
    return fast ? RSIMD_KERNEL(cmath2_)(RSIMD_CM2_LOGB, 1, x, y, n, flags, out, b, a)
                : RSIMD_KERNEL(cmath2_)(RSIMD_CM2_LOGB, 0, x, y, n, flags, out, b, a);
  default:
    return fast ? RSIMD_KERNEL(cmath2_)(RSIMD_CM2_ATAN2, 1, x, y, n, flags, out, b, a)
                : RSIMD_KERNEL(cmath2_)(RSIMD_CM2_ATAN2, 0, x, y, n, flags, out, b, a);
  }
}

#undef RSIMD_CM_V
#undef RSIMD_CM_F1
#undef RSIMD_CM_F2

#else /* 32-bit ARM, or a SIMD tier without SLEEF: filled from below */
#define RSIMD_SKIP_cmath1_c128 1
#define RSIMD_SKIP_cmath2_c128 1
#endif
