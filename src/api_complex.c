/* Complex numbers: the .Call entry points of Conj, Re and Im, Mod and Arg
   and the load-time choice of the arithmetic variants, and the complex
   paths of the arithmetic, reduction, scan, comparison and predicate
   entry points.

   Addition, subtraction and negation reuse the double kernels on the 2n
   doubles of the interleaved Rcomplex layout, so each part follows the
   double rules (NA payloads, na_check). A complex scalar operand is two
   doubles, which the double kernels' single-scalar broadcast cannot
   express: it is repeated into a small block that the kernel reads as an
   ordinary operand. Results are bare vectors. */

#include <complex.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "soft_fma.h"
#include "api_complex.h"

/* Complex numbers per broadcast block: 4 KiB on the stack. */
#define RSIMD_C128_BCAST 256

/* x op y (RSIMD_EW_ADD or RSIMD_EW_SUB) for one chunk of len elements;
   xs or ys is set when that operand is a scalar. */
static void c128_add_chunk(int op, const Rcomplex *px, const Rcomplex *py, R_xlen_t len, int xs,
                           int ys, Rcomplex *out, const rsimd_opts *o) {
  Rcomplex blk[RSIMD_C128_BCAST];
  R_xlen_t i;
  if (!xs && !ys) {
    rsimd_active->ew2_f64(op, rsimd_c128_as_f64(px, len), rsimd_c128_as_f64(py, len), 2 * len, 0,
                          (double *) out, o);
    return;
  }
  for (i = 0; i < RSIMD_C128_BCAST; i++) blk[i] = xs ? px[0] : py[0];
  for (i = 0; i < len; i += RSIMD_C128_BCAST) {
    R_xlen_t l = len - i < RSIMD_C128_BCAST ? len - i : RSIMD_C128_BCAST;
    rsimd_active->ew2_f64(op, (const double *) (xs ? blk : px + i),
                          (const double *) (ys ? blk : py + i), 2 * l, 0, (double *) (out + i), o);
  }
}

SEXP rsimd_c128_add(int op, SEXP x, SEXP y, const rsimd_opts *o) {
  rsimd_bin b;
  SEXP out;
  Rcomplex *po;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
    Rf_error("internal error: complex arithmetic on %s and %s", rsimd_etype_names[b.x.type],
             rsimd_etype_names[b.y.type]);
  }
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, b.n));
  po = (Rcomplex *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off, {
    c128_add_chunk(op, px, py, len, b.x_scalar, b.y_scalar, po + off, o);
  });
  UNPROTECT(1);
  return out;
}

/* ---- Multiply, divide and cumprod as base R computes them --------------

   Base R multiplies and divides complex numbers with C99's double complex
   operators and steps cumprod with plain double arithmetic, so the result
   depends on the compiler that built R: whether it contracts a*b + c into
   a fused multiply-add, in which operand order (which decides whether NA
   or NaN wins where both operands have one), and which __divdc3 its
   runtime has. The functions below are those operators compiled with the
   package's flags, which are R's (the tier units turn contraction off).
   When the package loads, R computes products, quotients and a cumprod of
   inputs from c128_probe_inputs() with base R's operators, and
   C_simd_c128_probe() picks the kernel variants and operators that
   reproduce all of them (rsimd_c128_arith in kernel_types.h). */

static double complex to_c99(const Rcomplex *x) {
  double complex z = 0;
  __real__ z = x->r;
  __imag__ z = x->i;
  return z;
}

static void from_c99(double complex z, Rcomplex *out) {
  out->r = creal(z);
  out->i = cimag(z);
}

/* x * y, and the same product with the operands the other way round. */
static void c99_mul_xy(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  from_c99(to_c99(x) * to_c99(y), out);
}
static void c99_mul_yx(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  from_c99(to_c99(y) * to_c99(x), out);
}
static void c99_div(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  from_c99(to_c99(x) / to_c99(y), out);
}
/* Base R's cumprod step, x the element and y the product so far. */
static void cumprod_step(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  Rcomplex t = *y;
  out->r = x->r * t.r - x->i * t.i;
  out->i = x->r * t.i + x->i * t.r;
}

/* ---- Base R's operators, independent of this file's compiler ------------

   The operators above are base R's when this file's compiler rounds as
   R's did, which it need not: Ubuntu's amd64v3 GCC targets x86-64-v3 and
   fuses a*b + c (and so does its libgcc's __muldc3 and __divdc3), while R
   from CRAN or Posit is built for baseline x86-64 and fuses nothing. So the
   probe also tries the operators below: C11 Annex G's multiply (G.5.1,
   with its recovery of infinities) and Smith's division with the scaling
   of GCC 12 and later, then Annex G's recovery, fused as libgcc's GCC
   build fuses them or not at all.

   Every rounding and every NaN is explicit, so that this file's compiler
   cannot change them: rd() rounds on its own (no contraction can fuse it),
   rsimd_fma fuses, and where an operand is NaN the result is the NaN that
   the hardware gives for the operands in the order libgcc's build has them
   (x86 the first operand's; Arm a signalling NaN's first, as R's NA is
   stored, then the first operand's), quieted. That order is the one GCC
   gives libgcc's source, which the operations below follow (ratio * b
   rather than b * ratio, say). */

static double rd(double x) {
  volatile double v = x;
  return v;
}

static int is_snan(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return isnan(x) && !((u >> 51) & 1);
}

static double quiet(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  u |= (uint64_t) 1 << 51;
  memcpy(&x, &u, sizeof x);
  return x;
}

/* The NaN of p op q, p or q NaN, p the first operand. */
static double nan2(double p, double q) {
#if defined(__aarch64__) || defined(__arm__)
  if (is_snan(p)) return quiet(p);
  if (is_snan(q)) return quiet(q);
#endif
  return quiet(isnan(p) ? p : q);
}

/* The NaN of fma(p, q, r), one of them NaN. Arm takes the addend first,
   and gives its default NaN for inf * 0 + a quiet NaN. */
static double nan3(double p, double q, double r) {
#if defined(__aarch64__) || defined(__arm__)
  if (is_snan(r)) return quiet(r);
  if (isnan(r) && ((isinf(p) && q == 0) || (p == 0 && isinf(q)))) return NAN;
  if (is_snan(p)) return quiet(p);
  if (is_snan(q)) return quiet(q);
  return quiet(isnan(r) ? r : isnan(p) ? p : q);
#else
  return quiet(isnan(p) ? p : isnan(q) ? q : r);
#endif
}

static double mul_(double p, double q) { return isnan(p) || isnan(q) ? nan2(p, q) : rd(p * q); }
static double add_(double p, double q) { return isnan(p) || isnan(q) ? nan2(p, q) : rd(p + q); }
static double sub_(double p, double q) { return isnan(p) || isnan(q) ? nan2(p, q) : rd(p - q); }
static double div_(double p, double q) { return isnan(p) || isnan(q) ? nan2(p, q) : rd(p / q); }
static double fma_(double p, double q, double r) {
  return isnan(p) || isnan(q) || isnan(r) ? nan3(p, q, r) : rsimd_fma(p, q, r);
}
/* p * q + r, p * q - r and r - p * q, fused or not. */
static double madd_(int fused, double p, double q, double r) {
  return fused ? fma_(p, q, r) : add_(mul_(p, q), r);
}
static double msub_(int fused, double p, double q, double r) {
  return fused ? fma_(p, q, isnan(r) ? r : -r) : sub_(mul_(p, q), r);
}
static double rsub_(int fused, double r, double p, double q) {
  return fused ? fma_(isnan(p) ? p : -p, q, r) : sub_(r, mul_(p, q));
}

/* (a + bi)(c + di): the parts ac - bd and ad + bc with every product
   rounded (each is also needed on its own), and where both are NaN Annex
   G's recovery: an infinite factor is boxed to +-1 and 0 parts and the
   NaNs of the other become zeros, as do all NaNs where a partial product
   overflowed, and the parts are recomputed times infinity, fused (as
   fma(a, c, -bd) and fma(a, d, bc)) or not. */
static void ag_mul(int fused, double a, double b, double c, double d, Rcomplex *out) {
  const double ac = mul_(a, c), bd = mul_(b, d), ad = mul_(a, d), cb = mul_(c, b);
  double x = sub_(ac, bd), y = add_(ad, cb);
  if (isnan(x) && isnan(y)) {
    int recalc = 0;
    if (isinf(a) || isinf(b)) {
      a = copysign(isinf(a) ? 1.0 : 0.0, a);
      b = copysign(isinf(b) ? 1.0 : 0.0, b);
      if (isnan(c)) c = copysign(0.0, c);
      if (isnan(d)) d = copysign(0.0, d);
      recalc = 1;
    }
    if (isinf(c) || isinf(d)) {
      c = copysign(isinf(c) ? 1.0 : 0.0, c);
      d = copysign(isinf(d) ? 1.0 : 0.0, d);
      if (isnan(a)) a = copysign(0.0, a);
      if (isnan(b)) b = copysign(0.0, b);
      recalc = 1;
    }
    if (!recalc && (isinf(ac) || isinf(bd) || isinf(ad) || isinf(cb))) {
      if (isnan(a)) a = copysign(0.0, a);
      if (isnan(b)) b = copysign(0.0, b);
      if (isnan(c)) c = copysign(0.0, c);
      if (isnan(d)) d = copysign(0.0, d);
      recalc = 1;
    }
    if (recalc) {
      x = INFINITY * (fused ? rsimd_fma(a, c, -rd(b * d)) : rd(a * c) - rd(b * d));
      y = INFINITY * (fused ? rsimd_fma(a, d, rd(b * c)) : rd(a * d) + rd(b * c));
    }
  }
  out->r = x;
  out->i = y;
}

/* (a + bi) / (c + di) by Smith's method: the operands halved when the
   divisor's larger part is near overflow, or scaled up by 2^52 when it is
   below DBL_EPSILON or a dividend part is subnormal against a small
   divisor; the ratio of the divisor's parts, applied to the quotient of a
   dividend part instead where it is subnormal or zero. Where both parts
   are NaN, Annex G's recovery (of the scaled operands): a zero divisor
   gives infinities, an infinite dividend over a finite divisor an
   infinity, a finite dividend over an infinite divisor zeros. */
static void ag_div(int fused, double a, double b, double c, double d, Rcomplex *out) {
  const double rbig = DBL_MAX / 2, rmax2 = rbig * DBL_EPSILON;
  const int sw = fabs(c) < fabs(d);
  double big = sw ? fabs(d) : fabs(c), ratio, denom, x, y;
  if (big >= rbig) {
    a = div_(a, 2);
    b = div_(b, 2);
    c = div_(c, 2);
    d = div_(d, 2);
    big /= 2;
  }
  if (big < DBL_EPSILON || (fabs(a) < DBL_MIN && fabs(b) < rmax2 && big < rmax2) ||
      (fabs(b) < DBL_MIN && fabs(a) < rmax2 && big < rmax2)) {
    a = mul_(a, 1 / DBL_EPSILON);
    b = mul_(b, 1 / DBL_EPSILON);
    c = mul_(c, 1 / DBL_EPSILON);
    d = mul_(d, 1 / DBL_EPSILON);
  }
  if (sw) {
    ratio = div_(c, d);
    denom = madd_(fused, c, ratio, d);
    if (fabs(ratio) > DBL_MIN) {
      x = madd_(fused, a, ratio, b);
      y = msub_(fused, ratio, b, a);
    } else {
      x = madd_(fused, div_(a, d), c, b);
      y = msub_(fused, div_(b, d), c, a);
    }
  } else {
    ratio = div_(d, c);
    denom = madd_(fused, d, ratio, c);
    if (fabs(ratio) > DBL_MIN) {
      x = madd_(fused, b, ratio, a);
      y = rsub_(fused, b, ratio, a);
    } else {
      x = madd_(fused, div_(b, c), d, a);
      y = rsub_(fused, b, div_(a, c), d);
    }
  }
  x = div_(x, denom);
  y = div_(y, denom);
  if (isnan(x) && isnan(y)) {
    if (c == 0.0 && d == 0.0 && (!isnan(a) || !isnan(b))) {
      x = mul_(copysign(INFINITY, c), a);
      y = mul_(copysign(INFINITY, c), b);
    } else if ((isinf(a) || isinf(b)) && isfinite(c) && isfinite(d)) {
      a = copysign(isinf(a) ? 1.0 : 0.0, a);
      b = copysign(isinf(b) ? 1.0 : 0.0, b);
      x = INFINITY * (rd(a * c) + rd(b * d));
      y = INFINITY * (rd(b * c) - rd(a * d));
    } else if ((isinf(c) || isinf(d)) && isfinite(a) && isfinite(b)) {
      c = copysign(isinf(c) ? 1.0 : 0.0, c);
      d = copysign(isinf(d) ? 1.0 : 0.0, d);
      x = 0.0 * (rd(a * c) + rd(b * d));
      y = 0.0 * (rd(b * c) - rd(a * d));
    }
  }
  out->r = x;
  out->i = y;
}

/* The operators as rsimd_c128_fn: unfused (u) or fused (f), x * y or with
   the operands the other way round (yx). */
static void ag_mul_u(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_mul(0, x->r, x->i, y->r, y->i, out);
}
static void ag_mul_u_yx(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_mul(0, y->r, y->i, x->r, x->i, out);
}
static void ag_mul_f(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_mul(1, x->r, x->i, y->r, y->i, out);
}
static void ag_mul_f_yx(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_mul(1, y->r, y->i, x->r, x->i, out);
}
static void ag_div_u(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_div(0, x->r, x->i, y->r, y->i, out);
}
static void ag_div_f(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  ag_div(1, x->r, x->i, y->r, y->i, out);
}

/* Until the probe has run: base R's operators for every element. */
static rsimd_c128_arith c128_arith = {RSIMD_CMUL_SCALAR, RSIMD_CMUL_SCALAR, RSIMD_CDIV_SCALAR,
                                      RSIMD_CMUL_SCALAR, RSIMD_CMUL_SCALAR, c99_mul_xy, c99_div,
                                      cumprod_step};
const rsimd_c128_arith *rsimd_c128_arith_get(void) { return &c128_arith; }

/* What the probe chose, for simd_current(): mul, div and cumprod. */
static char c128_names[3][96] = {"scalar", "scalar", "scalar"};

/* A deterministic sequence for the probe's inputs (splitmix64), so that
   loading the package does not touch R's random number generator. */
static uint64_t probe_next(uint64_t *s) {
  uint64_t z = (*s += 0x9E3779B97F4A7C15u);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

/* A random double with exponent in [lo, hi] (the full significand random,
   a random sign); lo below -1022 gives subnormals. */
static double probe_double(uint64_t *s, int lo, int hi) {
  uint64_t r = probe_next(s);
  double m = 1.0 + (double) (r >> 12) / 4503599627370496.0;
  int e = lo + (int) (probe_next(s) % (uint64_t) (hi - lo + 1));
  return (r & 1 ? -1.0 : 1.0) * ldexp(m, e);
}

#define PROBE_RANDOM 1024
#define PROBE_SPECIALS 8
#define PROBE_CUMPROD 512
#define PROBE_ASIN 8

/* The probe's operands: list(z, w, x, a) with z and w for x * y and x / y
   (random pairs of moderate and of extreme magnitudes, then every
   combination of special values for the four parts), x for cumprod
   (random numbers of modulus near 1, so the product stays finite) and a
   for asin (real numbers just above 1, where the roundings of alpha^2 - 1
   differ). */
static SEXP simd_c128_probe_inputs_impl(void) {
  const double sp[PROBE_SPECIALS] = {1.0, 0.0, -0.0, -2.5, R_PosInf, R_NegInf, R_NaN, NA_REAL};
  const R_xlen_t ns = PROBE_SPECIALS * PROBE_SPECIALS * PROBE_SPECIALS * PROBE_SPECIALS;
  const R_xlen_t n = 2 * PROBE_RANDOM + ns;
  uint64_t seed = 20261005u;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4)), names;
  Rcomplex *z, *w, *x, *a;
  R_xlen_t i;
  SET_VECTOR_ELT(out, 0, rsimd_alloc_like(RSIMD_C128, n));
  SET_VECTOR_ELT(out, 1, rsimd_alloc_like(RSIMD_C128, n));
  SET_VECTOR_ELT(out, 2, rsimd_alloc_like(RSIMD_C128, PROBE_CUMPROD));
  SET_VECTOR_ELT(out, 3, rsimd_alloc_like(RSIMD_C128, PROBE_ASIN));
  z = (Rcomplex *) rsimd_out_ptr(VECTOR_ELT(out, 0));
  w = (Rcomplex *) rsimd_out_ptr(VECTOR_ELT(out, 1));
  x = (Rcomplex *) rsimd_out_ptr(VECTOR_ELT(out, 2));
  a = (Rcomplex *) rsimd_out_ptr(VECTOR_ELT(out, 3));
  for (i = 0; i < 2 * PROBE_RANDOM; i++) {
    /* The first half has moderate magnitudes, where only the rounding
       differs between the variants; the second spans the exponent range
       (overflow, subnormals and every scaling case of the division), with
       some zero parts. */
    int lo = i < PROBE_RANDOM ? -8 : -1074, hi = i < PROBE_RANDOM ? 8 : 1023;
    z[i].r = probe_double(&seed, lo, hi);
    z[i].i = probe_double(&seed, lo, hi);
    w[i].r = probe_double(&seed, lo, hi);
    w[i].i = probe_double(&seed, lo, hi);
    if (i >= PROBE_RANDOM && i % 5 == 0) z[i].i = 0.0;
    if (i >= PROBE_RANDOM && i % 7 == 0) w[i].i = 0.0;
    if (i >= PROBE_RANDOM && i % 11 == 0) w[i].r = 0.0;
  }
  for (i = 0; i < ns; i++) {
    R_xlen_t k = i;
    z[2 * PROBE_RANDOM + i].r = sp[k % PROBE_SPECIALS];
    k /= PROBE_SPECIALS;
    z[2 * PROBE_RANDOM + i].i = sp[k % PROBE_SPECIALS];
    k /= PROBE_SPECIALS;
    w[2 * PROBE_RANDOM + i].r = sp[k % PROBE_SPECIALS];
    k /= PROBE_SPECIALS;
    w[2 * PROBE_RANDOM + i].i = sp[k % PROBE_SPECIALS];
  }
  for (i = 0; i < PROBE_CUMPROD; i++) {
    double t = probe_double(&seed, -2, 1);
    x[i].r = cos(t) * (1.0 + probe_double(&seed, -12, -10));
    x[i].i = sin(t) * (1.0 + probe_double(&seed, -12, -10));
  }
  for (i = 0; i < PROBE_ASIN; i++) {
    a[i].r = 1 + ldexp((double) (2 * i + 3), -30);
    a[i].i = 0.0;
  }
  names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, Rf_mkChar("z"));
  SET_STRING_ELT(names, 1, Rf_mkChar("w"));
  SET_STRING_ELT(names, 2, Rf_mkChar("x"));
  SET_STRING_ELT(names, 3, Rf_mkChar("a"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP C_simd_c128_probe_inputs(void) {
  rsimd_entry();
  return rsimd_exit(simd_c128_probe_inputs_impl());
}

/* One part as identical() compares it: NA and NaN are told apart (other
   NaN payloads and signs are not), and a zero's sign counts. With
   `nan_blind`, any two NaNs are equal. */
static int same_part(double a, double b, int nan_blind) {
  if (isnan(a) || isnan(b)) {
    return isnan(a) && isnan(b) && (nan_blind || R_IsNA(a) == R_IsNA(b));
  }
  return memcmp(&a, &b, sizeof a) == 0;
}

/* 2 if out matches base in every element, 1 if only NA versus NaN
   differs, 0 otherwise. */
static int match_level(const Rcomplex *out, const Rcomplex *base, R_xlen_t n) {
  int level = 2;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (same_part(out[i].r, base[i].r, 0) && same_part(out[i].i, base[i].i, 0)) continue;
    if (!same_part(out[i].r, base[i].r, 1) || !same_part(out[i].i, base[i].i, 1)) return 0;
    level = 1;
  }
  return level;
}

/* What a kernel computes: the formula (op, v1, v2 as formula_c128), and
   f where that gives a NaN part; with v1 SCALAR, f alone. */
static void probe_run(int op, int v1, int v2, rsimd_c128_fn f, const Rcomplex *x,
                      const Rcomplex *y, R_xlen_t n, Rcomplex *out) {
  R_xlen_t i;
  if (v1 != RSIMD_CMUL_SCALAR) rsimd_active->formula_c128(op, v1, v2, x, y, n, out);
  for (i = 0; i < n; i++) {
    if (v1 == RSIMD_CMUL_SCALAR || isnan(out[i].r) || isnan(out[i].i)) f(x + i, y + i, out + i);
  }
}

static const char *const mul_re_names[] = {"ac - bd", "fma(a, c, -bd)", "fma(-b, d, ac)"};
static const char *const mul_im_names[] = {"ad + bc", "fma(b, c, ad)", "fma(a, d, bc)"};

/* The name of choice k; `own` when its NaN results come from the Annex G
   operators above rather than this compiler's. */
static void set_name(int k, const char *s, int own, int level) {
  snprintf(c128_names[k], sizeof c128_names[k], "%s%s%s", s, own ? "; Annex G" : "",
           level == 1 ? " (NA/NaN may differ)" : "");
}

static void set_mul_name(int k, int vr, int vi, int own, int level) {
  char buf[64];
  if (vr == RSIMD_CMUL_SCALAR) {
    set_name(k, "scalar", own, level);
    return;
  }
  snprintf(buf, sizeof buf, "%s, %s", mul_re_names[vr], mul_im_names[vi]);
  set_name(k, buf, own, level);
}

/* The variants, in order of preference: those of GCC builds first. */
static const int variants[9][2] = {
  {RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1},     {RSIMD_CMUL_UNFUSED, RSIMD_CMUL_UNFUSED},
  {RSIMD_CMUL_FMA2, RSIMD_CMUL_FMA1},     {RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA2},
  {RSIMD_CMUL_FMA2, RSIMD_CMUL_FMA2},     {RSIMD_CMUL_FMA1, RSIMD_CMUL_UNFUSED},
  {RSIMD_CMUL_UNFUSED, RSIMD_CMUL_FMA1},  {RSIMD_CMUL_FMA2, RSIMD_CMUL_UNFUSED},
  {RSIMD_CMUL_UNFUSED, RSIMD_CMUL_FMA2}};

/* Picks the variants from base R's results: prod = z * w, quot = z / w,
   cp = cumprod(x). Each op takes the first variant (with the first
   operator) that matches base R exactly, else the first that matches up
   to NA versus NaN; else its operator alone (SCALAR). Returns the names
   of the choices, as c128_variants(). */
/* The elements of a probe vector: complex, of length n (or any length
   with n < 0, which is then returned in *len), and contiguous. */
static const Rcomplex *probe_data(SEXP v, R_xlen_t n, R_xlen_t *len) {
  rsimd_in in;
  rsimd_in_init(&in, v, "probe");
  if (in.type != RSIMD_C128 || in.ptr == NULL || in.n < 1 || (n >= 0 && in.n != n)) {
    Rf_error("internal error: invalid complex probe");
  }
  if (len != NULL) *len = in.n;
  return (const Rcomplex *) in.ptr;
}

static SEXP simd_c128_probe_impl(SEXP z, SEXP w, SEXP prod, SEXP quot, SEXP x, SEXP cp,
                                 SEXP as, SEXP asin_base) {
  /* For each variant, the operators for its NaN results in order: this
     compiler's, then the Annex G ones (unfused first where the variant is,
     as an unfused build's are; fused first otherwise). */
  const rsimd_c128_fn mul_fns[2][6] = {
    {c99_mul_xy, c99_mul_yx, ag_mul_u, ag_mul_u_yx, ag_mul_f, ag_mul_f_yx},
    {c99_mul_xy, c99_mul_yx, ag_mul_f, ag_mul_f_yx, ag_mul_u, ag_mul_u_yx}};
  /* For testing, RSIMD_DEBUG_OWN_COMPLEX_OPS=true leaves out this
     compiler's operators, so that the Annex G ones are checked against
     base R wherever the tests run. */
  const char *own_env = getenv("RSIMD_DEBUG_OWN_COMPLEX_OPS");
  const int f0 = own_env != NULL && strcmp(own_env, "true") == 0 ? 1 : 0;
  const Rcomplex *pz, *pw, *pp, *pq, *px, *pc;
  Rcomplex *out, *prev, one;
  rsimd_c128_arith a = c128_arith;
  int best = 0, k, f, level;
  R_xlen_t i, n, m;

  pz = probe_data(z, -1, &n);
  pw = probe_data(w, n, NULL);
  pp = probe_data(prod, n, NULL);
  pq = probe_data(quot, n, NULL);
  px = probe_data(x, -1, &m);
  pc = probe_data(cp, m, NULL);
  {
    R_xlen_t na;
    const Rcomplex *pa = probe_data(as, -1, &na);
    rsimd_cmath_probe_asin(pa, probe_data(asin_base, na, NULL), na);
  }
  out = (Rcomplex *) R_alloc((size_t) (n > m ? n : m), sizeof(Rcomplex));

  /* x * y: the variant and the operator for its NaN results (mul_fns). */
  a.mul_re = a.mul_im = RSIMD_CMUL_SCALAR;
  a.mul1 = f0 ? ag_mul_u : c99_mul_xy;
  for (k = 0; k <= 9 && best < 2; k++) {
    int vr = k < 9 ? variants[k][0] : RSIMD_CMUL_SCALAR, vi = k < 9 ? variants[k][1] : vr;
    const rsimd_c128_fn *fns = mul_fns[vr != RSIMD_CMUL_UNFUSED && vr != RSIMD_CMUL_SCALAR];
    for (f = 2 * f0; f < 6; f++) {
      probe_run(RSIMD_EW_MUL, vr, vi, fns[f], pz, pw, n, out);
      level = match_level(out, pp, n);
      if (level > best) {
        best = level;
        a.mul_re = vr;
        a.mul_im = vi;
        a.mul1 = fns[f];
      }
      if (best == 2) break;
    }
  }
  set_mul_name(0, a.mul_re, a.mul_im, a.mul1 != c99_mul_xy && a.mul1 != c99_mul_yx, best);

  /* x / y: libgcc's algorithm fused or not, else the operator alone; for
     the NaN results this compiler's operator, then the Annex G one rounded
     as the variant (both ways for SCALAR). */
  best = 0;
  a.div = RSIMD_CDIV_SCALAR;
  a.div1 = f0 ? ag_div_u : c99_div;
  for (k = 0; k < 3 && best < 2; k++) {
    int v = k == 0 ? RSIMD_CDIV_FMA : k == 1 ? RSIMD_CDIV_UNFUSED : RSIMD_CDIV_SCALAR;
    const rsimd_c128_fn div_fns[3] = {c99_div, v == RSIMD_CDIV_FMA ? ag_div_f : ag_div_u, ag_div_f};
    for (f = f0; f < (v == RSIMD_CDIV_SCALAR ? 3 : 2) && best < 2; f++) {
      probe_run(RSIMD_EW_DIV, v == RSIMD_CDIV_SCALAR ? RSIMD_CMUL_SCALAR : v, 0, div_fns[f], pz, pw,
                n, out);
      level = match_level(out, pq, n);
      if (level > best) {
        best = level;
        a.div = v;
        a.div1 = div_fns[f];
      }
    }
  }
  set_name(1, a.div == RSIMD_CDIV_FMA        ? "libgcc, fused"
              : a.div == RSIMD_CDIV_UNFUSED ? "libgcc, unfused"
                                            : "scalar",
           a.div1 != c99_div, best);

  /* cumprod: each step is the element times the previous product. */
  prev = (Rcomplex *) R_alloc((size_t) m, sizeof(Rcomplex));
  one.r = 1.0;
  one.i = 0.0;
  for (i = 0; i < m; i++) prev[i] = i == 0 ? one : pc[i - 1];
  best = 0;
  a.cp_re = a.cp_im = RSIMD_CMUL_SCALAR;
  for (k = 0; k <= 9 && best < 2; k++) {
    int vr = k < 9 ? variants[k][0] : RSIMD_CMUL_SCALAR, vi = k < 9 ? variants[k][1] : vr;
    if (vr == RSIMD_CMUL_SCALAR) {
      for (i = 0; i < m; i++) cumprod_step(px + i, prev + i, out + i);
    } else {
      rsimd_active->formula_c128(RSIMD_EW_MUL, vr, vi, px, prev, m, out);
    }
    level = match_level(out, pc, m);
    if (level > best) {
      best = level;
      a.cp_re = vr;
      a.cp_im = vi;
    }
  }
  set_mul_name(2, a.cp_re, a.cp_im, 0, best);

  c128_arith = a;
  return c128_variants();
}

SEXP C_simd_c128_probe(SEXP z, SEXP w, SEXP prod, SEXP quot, SEXP x, SEXP cp, SEXP as,
                       SEXP asin_base) {
  rsimd_entry();
  return rsimd_exit(simd_c128_probe_impl(z, w, prod, quot, x, cp, as, asin_base));
}

SEXP c128_variants(void) {
  static const char *const names[] = {"mul", "div", "cumprod"};
  SEXP out = PROTECT(Rf_allocVector(STRSXP, 3)), nm = PROTECT(Rf_allocVector(STRSXP, 3));
  int k;
  for (k = 0; k < 3; k++) {
    SET_STRING_ELT(out, k, Rf_mkChar(c128_names[k]));
    SET_STRING_ELT(nm, k, Rf_mkChar(names[k]));
  }
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

SEXP C_simd_c128_variants(void) {
  rsimd_entry();
  return rsimd_exit(c128_variants());
}

/* Internal, for the tests: sets the variants by code (mul_re, mul_im,
   div, cp_re, cp_im; -1 SCALAR) with the default operators, so that the
   kernels of a variant this R build does not use can be checked against
   formula_c128. Returns the previous codes. */
static SEXP simd_c128_set_variants_impl(SEXP codes) {
  SEXP old;
  const int *pc;
  int *po;
  rsimd_in in;
  rsimd_in_init(&in, codes, "codes");
  if (in.type != RSIMD_I32 || in.n != 5 || in.ptr == NULL) {
    Rf_error("internal error: invalid variant codes");
  }
  pc = (const int *) in.ptr;
  old = PROTECT(rsimd_alloc_like(RSIMD_I32, 5));
  po = (int *) rsimd_out_ptr(old);
  po[0] = c128_arith.mul_re;
  po[1] = c128_arith.mul_im;
  po[2] = c128_arith.div;
  po[3] = c128_arith.cp_re;
  po[4] = c128_arith.cp_im;
  c128_arith.mul_re = pc[0];
  c128_arith.mul_im = pc[1];
  c128_arith.div = pc[2];
  c128_arith.cp_re = pc[3];
  c128_arith.cp_im = pc[4];
  UNPROTECT(1);
  return old;
}

SEXP C_simd_c128_set_variants(SEXP codes) {
  rsimd_entry();
  return rsimd_exit(simd_c128_set_variants_impl(codes));
}

/* x * y or x / y (op RSIMD_EW_MUL or RSIMD_EW_DIV) of two complex
   vectors, bit-identical to base R. */
SEXP rsimd_c128_muldiv(int op, SEXP x, SEXP y) {
  rsimd_bin b;
  SEXP out;
  Rcomplex *po;
  int flags;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
    Rf_error("internal error: complex arithmetic on %s and %s", rsimd_etype_names[b.x.type],
             rsimd_etype_names[b.y.type]);
  }
  flags = (b.x_scalar ? RSIMD_EW_SCALAR(0) : 0) | (b.y_scalar ? RSIMD_EW_SCALAR(1) : 0);
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, b.n));
  po = (Rcomplex *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off, {
    rsimd_active->ew2_c128(op, px, py, len, flags, po + off, &c128_arith);
  });
  UNPROTECT(1);
  return out;
}

/* One complex number as a length-1 complex vector. */
static SEXP c128_scalar(double re, double im) {
  SEXP out = PROTECT(rsimd_alloc_like(RSIMD_C128, 1));
  Rcomplex *po = (Rcomplex *) rsimd_out_ptr(out);
  po->r = re;
  po->i = im;
  UNPROTECT(1);
  return out;
}

/* The R value of the product s. */
static SEXP cprod_finish(rsimd_cprod_state *s, const rsimd_opts *o);

SEXP rsimd_c128_prod(const rsimd_in *in, const rsimd_opts *o) {
  rsimd_cprod_state s;

  memset(&s, 0, sizeof s);
  s.p.r = 1.0;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off,
                      { rsimd_active->prod_c128(px, len, &s, &c128_arith, o); });
  return cprod_finish(&s, o);
}

SEXP rsimd_c128_prod2(int op, const rsimd_bin *b, const rsimd_opts *o) {
  rsimd_cprod_state s;
  Rcomplex blk[RSIMD_C128_BCAST];

  memset(&s, 0, sizeof s);
  s.p.r = 1.0;
  RSIMD_FOREACH_CHUNK2(b, Rcomplex, px, py, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i += RSIMD_C128_BCAST) {
      R_xlen_t l = len - i < RSIMD_C128_BCAST ? len - i : RSIMD_C128_BCAST;
      c128_add_chunk(op, b->x_scalar ? px : px + i, b->y_scalar ? py : py + i, l, b->x_scalar,
                     b->y_scalar, blk, o);
      rsimd_active->prod_c128(blk, l, &s, &c128_arith, o);
    }
  });
  return cprod_finish(&s, o);
}

static SEXP cprod_finish(rsimd_cprod_state *s, const rsimd_opts *o) {
  Rcomplex r;
  int k;

  /* With missing values (and na.rm = FALSE), both parts are NA if one was
     NA, else NaN. */
  if (s->saw_nan && !o->na_rm) {
    double v = s->saw_na ? NA_REAL : R_NaN;
    return c128_scalar(v, v);
  }
  r = s->p;
  if (o->precision == RSIMD_PREC_COMPENSATED) {
    r.r += s->lo.r;
    r.i += s->lo.i;
  } else if (o->precision == RSIMD_PREC_PAIRWISE) {
    /* The unfinished leaf, then the leaves from the latest to the
       earliest. */
    for (k = 0; k < 64; k++) {
      if ((s->leaves >> k) & 1) {
        rsimd_active->ew2_c128(RSIMD_EW_MUL, s->pw + k, &r, 1, 0, &r, &c128_arith);
      }
    }
  }
  return c128_scalar(r.r, r.i);
}

static double part_value(const rsimd_reduce_result *r, int precision);

SEXP rsimd_c128_mean(const rsimd_in *in, const rsimd_opts *o) {
  rsimd_reduce_result r[2];
  int k;
  double v[2];

  rsimd_reduce_result_init(&r[0], RSIMD_RED_SUM);
  rsimd_reduce_result_init(&r[1], RSIMD_RED_SUM);
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, { rsimd_active->sum_c128(px, len, r, o); });
  for (k = 0; k < 2; k++) {
    v[k] = part_value(&r[k], o->precision);
    if (!isnan(v[k])) v[k] /= (double) r[k].count;
  }
  return c128_scalar(v[0], v[1]);
}

/* Base R's NA/NaN fix-up of a complex cumsum or cumprod (cum.c's
   chandleNaN): when the last value's real (imaginary) part is NaN, that
   part of every value from the first element with a NaN part on becomes
   NA if an NA part has been seen by then, else NaN. */
static void c128_scan_nan(const rsimd_in *in, Rcomplex *po) {
  const int r_nan = isnan(po[in->n - 1].r), i_nan = isnan(po[in->n - 1].i);
  int has_nan = 0, has_na = 0;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i++) {
      Rcomplex *v = po + off + i;
      has_nan = has_nan || isnan(px[i].r) || isnan(px[i].i);
      has_na = has_na || (has_nan && (R_IsNA(px[i].r) || R_IsNA(px[i].i)));
      if (has_na) {
        if (r_nan) v->r = NA_REAL;
        if (i_nan) v->i = NA_REAL;
      } else if (has_nan) {
        if (r_nan) v->r = R_NaN;
        if (i_nan) v->i = R_NaN;
      }
    }
  });
}

SEXP rsimd_c128_scan(const rsimd_in *in, int op) {
  SEXP out = PROTECT(rsimd_alloc_like(RSIMD_C128, in->n));
  Rcomplex *po = (Rcomplex *) rsimd_out_ptr(out), s;

  s.r = op == 0 ? 0.0 : 1.0;
  s.i = 0.0;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off,
                      { rsimd_active->scan_c128(op, px, len, po + off, &s, &c128_arith); });
  if (in->n > 0 && (isnan(s.r) || isnan(s.i))) c128_scan_nan(in, po);
  UNPROTECT(1);
  return out;
}

void rsimd_c128_cmp(int op, SEXP x, SEXP y, int *po) {
  rsimd_bin b;
  int flags;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
    Rf_error("internal error: complex comparison of %s and %s", rsimd_etype_names[b.x.type],
             rsimd_etype_names[b.y.type]);
  }
  flags = (b.x_scalar ? RSIMD_EW_SCALAR(0) : 0) | (b.y_scalar ? RSIMD_EW_SCALAR(1) : 0);
  RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off,
                       { rsimd_active->cmp_c128(op, px, py, len, flags, po + off); });
}

/* Mod(z) ("mod") or Arg(z) ("arg") in accuracy mode `accuracy` (0
   accurate, 1 fast). Complex z only; the R side handles the others. */
static SEXP simd_cmath_impl(SEXP z, SEXP op, SEXP accuracy) {
  static const char *const names[] = {"mod", "arg"};
  int code = rsimd_arg_choice(op, names, 2), a = rsimd_arg_accuracy(accuracy);
  SEXP out;
  double *po;
  rsimd_in in;

  rsimd_in_init(&in, z, "z");
  if (in.type != RSIMD_C128) Rf_error("internal error: Mod/Arg of %s", rsimd_etype_names[in.type]);
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
  po = (double *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off, {
    rsimd_active->math1_c128(code | (a ? RSIMD_MATH_FAST : 0), px, len, po + off);
  });
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cmath(SEXP z, SEXP op, SEXP accuracy) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_cmath_impl(z, op, accuracy), 0));
}

SEXP rsimd_c128_neg(const rsimd_in *in) {
  SEXP out = PROTECT(rsimd_alloc_like(RSIMD_C128, in->n));
  Rcomplex *po = (Rcomplex *) rsimd_out_ptr(out);
  rsimd_opts o;

  rsimd_opts_init_fixed(&o, 0, 1, 0);
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    rsimd_active->ew1_f64(RSIMD_EW_NEG, rsimd_c128_as_f64(px, len), 2 * len, 0,
                          (double *) (po + off), &o);
  });
  UNPROTECT(1);
  return out;
}

/* The value of one part of a finished sum: NA if that part saw an NA,
   else NaN if it saw a NaN, else the accumulated value (missing values
   are only tracked with na.rm = FALSE and na_check). */
static double part_value(const rsimd_reduce_result *r, int precision) {
  if (r->saw_na) return rsimd_na_real();
  if (r->saw_nan) return R_NaN;
  return rsimd_reduce_value(r, precision);
}

SEXP rsimd_c128_sum(const rsimd_in *in, const rsimd_opts *o) {
  rsimd_reduce_result r[2];
  SEXP out;
  Rcomplex *po;

  rsimd_reduce_result_init(&r[0], RSIMD_RED_SUM);
  rsimd_reduce_result_init(&r[1], RSIMD_RED_SUM);
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, { rsimd_active->sum_c128(px, len, r, o); });
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, 1));
  po = (Rcomplex *) rsimd_out_ptr(out);
  po->r = part_value(&r[0], o->precision);
  po->i = part_value(&r[1], o->precision);
  UNPROTECT(1);
  return out;
}

int rsimd_c128_pred(const rsimd_in *in, int code, int mode, int *po) {
  int res = mode == RSIMD_PRED_ALL;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    int r = rsimd_active->pred_c128(code, px, len, mode, mode == RSIMD_PRED_ELT ? po + off : NULL);
    if (mode != RSIMD_PRED_ELT && r != res) {
      res = r;
      break;
    }
  });
  return res;
}

/* Conj(z), Re(z) or Im(z) by name. Complex z gives a complex conjugate
   and double parts; double, integer and logical z give double, as in base
   R: Conj and Re are the values (a bare copy), Im is zeros. Raw is an
   error (rsimd_type_error()); integer64 is rejected on the R side. */
static SEXP simd_cplx_impl(SEXP z, SEXP op) {
  static const char *const names[] = {"conj", "re", "im"};
  int code = rsimd_arg_choice(op, names, 3);
  SEXP out;
  rsimd_in in;

  rsimd_in_init(&in, z, "z");
  switch (in.type) {
  case RSIMD_C128:
    if (code == 0) {
      Rcomplex *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_C128, in.n));
      po = (Rcomplex *) rsimd_out_ptr(out);
      RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off,
                          { rsimd_active->conj_c128(px, len, po + off); });
    } else {
      double *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
      po = (double *) rsimd_out_ptr(out);
      RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off,
                          { rsimd_active->part_c128(code == 2, px, len, po + off); });
    }
    break;
  case RSIMD_F64:
  case RSIMD_I32:
  case RSIMD_LGL: {
    double *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
    po = (double *) rsimd_out_ptr(out);
    if (code == 2) {
      if (in.n > 0) memset(po, 0, (size_t) in.n * sizeof(double));
    } else if (in.type == RSIMD_F64) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                          { memcpy(po + off, px, (size_t) len * sizeof(double)); });
    } else {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        rsimd_active->convert(RSIMD_CVT_I32_F64, RSIMD_CVT_CHECKED, px, len, po + off);
      });
    }
    break;
  }
  case RSIMD_U8: rsimd_type_error(RSIMD_U8, "z", "non-numeric argument to function");
  default: rsimd_type_error(in.type, "z", NULL);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cplx(SEXP z, SEXP op) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_cplx_impl(z, op), 1));
}
