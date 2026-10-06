/* Complex elementary functions: base R's own computation of each function
 * for one element (the none tier, and the vector tiers' fallback), and the
 * .Call entry points.
 *
 * The functions below repeat base R's src/main/complex.c: most call the
 * platform's C99 functions (csqrt, cexp, clog, ...), some are R's own code
 * (z_tan, z_asin, z_acos, z_atan, the inverse hyperbolic functions through
 * the trigonometric ones, and the power, log-with-base and atan2 code).
 * This file is compiled with R's compiler flags rather than a tier's (in
 * particular with the compiler's default floating-point contraction), as
 * base R's own file is, so the same expressions round the same way.
 * Products and quotients in whole-number powers use the multiply and divide
 * chosen at load time (rsimd_c128_arith), as the kernels do, which
 * reproduce base R's results whichever compiler built R.
 *
 * acosh is the exception: base R's acosh(z) is acos(z) * i, which is the
 * negative of the principal value wherever its real part is negative (the
 * lower half-plane and real z > 1). rsimd uses C99's cacosh instead; on
 * Windows, whose cacosh loses the sign and accuracy near the imaginary axis,
 * it is built from cacos as the vector kernels build it.
 *
 * On Windows, base R has no working ctanh and computes tanh as
 * -i tan(iz); so does rsimd there. */

#include <complex.h>
#include <math.h>
#include <string.h>
#include <Rmath.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"
#include "api_complex.h"

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

/* ---- Base R's functions --------------------------------------------------- */

/* z * w and z / w with base R's operators, as the kernels compute them:
   the load-time variant's formula (formula_c128, compiled without
   contraction), and mul1 or div1 where that gives a NaN part. mul1 and
   div1 alone are this compiler's operators, which round differently from
   base R's when R was built by another compiler. */
static double complex carith(int op, double complex z, double complex w) {
  const rsimd_c128_arith *ar = rsimd_c128_arith_get();
  const int v1 = op == RSIMD_EW_MUL ? ar->mul_re : ar->div;
  Rcomplex a, b, r;
  from_c99(z, &a);
  from_c99(w, &b);
  if (op == RSIMD_EW_MUL ? v1 != RSIMD_CMUL_SCALAR : v1 != RSIMD_CDIV_SCALAR) {
    rsimd_active->formula_c128(op, v1, ar->mul_im, &a, &b, 1, &r);
    if (!isnan(r.r) && !isnan(r.i)) return to_c99(&r);
  }
  (op == RSIMD_EW_MUL ? ar->mul1 : ar->div1)(&a, &b, &r);
  return to_c99(&r);
}

static double complex cmul(double complex z, double complex w) {
  return carith(RSIMD_EW_MUL, z, w);
}

static double complex R_cpow_n(double complex X, int k) {
  if (k == 0) return (double complex) 1.;
  else if (k == 1) return X;
  else if (k < 0) return carith(RSIMD_EW_DIV, 1., R_cpow_n(X, -k));
  else { /* k > 0 */
    double complex z = (double complex) 1.;
    while (k > 0) {
      if (k & 1) z = cmul(z, X);
      if (k == 1) break;
      k >>= 1; /* efficient division by 2; now have k >= 1 */
      X = cmul(X, X);
    }
    return z;
  }
}

/* Base R's X^Y. Windows builds of R have no cpow and use the polar form. */
static double complex mycpow(double complex X, double complex Y) {
  double complex Z;
  double yr = creal(Y), yi = cimag(Y);
  int k;
  if (X == 0.0) {
    if (yi == 0.0) Z = R_pow(0.0, yr);
    else Z = R_NaN + R_NaN * I;
  } else if (yi == 0.0 && fabs(yr) <= 65536 && yr == (k = (int) yr)) {
    /* range check before the cast: (int) of NaN or a huge yr is UB */
    Z = R_cpow_n(X, k);
  } else {
#ifndef _WIN32
    Z = cpow(X, Y);
#else
    double rho, r, i, theta;
    r = hypot(creal(X), cimag(X));
    i = atan2(cimag(X), creal(X));
    theta = i * yr;
    if (yi == 0.0) rho = pow(r, yr);
    else {
      /* rearrangement of cexp(X * clog(Y)) */
      r = log(r);
      theta += r * yi;
      rho = exp(r * yr - i * yi);
    }
    __real__ Z = rho * cos(theta);
    __imag__ Z = rho * sin(theta);
#endif
  }
  return Z;
}

static double complex z_tan(double complex z) {
  double y = cimag(z);
  double complex r = ctan(z);
  if (R_FINITE(y) && fabs(y) > 25.0) {
    /* at this point the real part is nearly zero, and the
       imaginary part is one: but some OSes get the imag as NaN */
    __imag__ r = y < 0 ? -1.0 : 1.0;
  }
  return r;
}

static double complex z_asin(double complex z) {
  if (cimag(z) == 0 && fabs(creal(z)) > 1) {
    double alpha, t1, t2, x = creal(z), ri;
    t1 = 0.5 * fabs(x + 1);
    t2 = 0.5 * fabs(x - 1);
    alpha = t1 + t2;
    ri = log(alpha + sqrt(alpha * alpha - 1));
    if (x > 1) ri *= -1;
    return asin(t1 - t2) + ri * I;
  }
  return casin(z);
}

static double complex z_acos(double complex z) {
  if (cimag(z) == 0 && fabs(creal(z)) > 1) return M_PI_2 - z_asin(z);
  return cacos(z);
}

static double complex z_atan(double complex z) {
  if (creal(z) == 0 && fabs(cimag(z)) > 1) {
    double y = cimag(z), rr, ri;
    rr = (y > 0) ? M_PI_2 : -M_PI_2;
    ri = 0.25 * log(((y + 1) * (y + 1)) / ((y - 1) * (y - 1)));
    return rr + ri * I;
  }
  return catan(z);
}

/* Base R's ctanh where the platform's is not working (R's R_ctanh). */
#ifdef _WIN32
static double complex z_tanh(double complex z) { return -I * z_tan(z * I); /* A&S 4.5.9 */ }
#else
#define z_tanh ctanh
#endif

/* C99's cacosh: Re >= 0 and the sign of Im(z) in the imaginary part. */
#ifdef _WIN32
static double complex z_cacosh(double complex z) {
  double complex a = cacos(z), r;
  __real__ r = fabs(cimag(a));
  __imag__ r = copysign(creal(a), cimag(z));
  return r;
}
#else
#define z_cacosh cacosh
#endif

static double complex z_asinh(double complex z) { return -I * z_asin(z * I); }

static double complex z_atanh(double complex z) { return -I * z_atan(z * I); }

static void z_logbase(const Rcomplex *z, const Rcomplex *base, Rcomplex *r) {
  from_c99(clog(to_c99(z)) / clog(to_c99(base)), r);
}

static void z_atan2(const Rcomplex *csn, const Rcomplex *ccs, Rcomplex *r) {
  double complex dr, dcsn = to_c99(csn), dccs = to_c99(ccs);
  if (dccs == 0) {
    if (dcsn == 0) {
      r->r = NA_REAL;
      r->i = NA_REAL; /* Why not R_NaN? */
      return;
    } else {
      double y = creal(dcsn);
      if (ISNAN(y)) dr = y;
      else dr = ((y >= 0) ? M_PI_2 : -M_PI_2);
    }
  } else {
    dr = catan(dcsn / dccs);
    if (creal(dccs) < 0) dr += M_PI;
    if (creal(dr) > M_PI) dr -= 2 * M_PI;
  }
  from_c99(dr, r);
}

#define CBASE1(name, expr)                                                                   \
  static void name(const Rcomplex *x, Rcomplex *out) {                                     \
    double complex z = to_c99(x);                                                          \
    from_c99(expr, out);                                                                   \
  }

CBASE1(b_sqrt, csqrt(z))
CBASE1(b_exp, cexp(z))
CBASE1(b_log, clog(z))
CBASE1(b_sin, csin(z))
CBASE1(b_cos, ccos(z))
CBASE1(b_tan, z_tan(z))
CBASE1(b_sinh, csinh(z))
CBASE1(b_cosh, ccosh(z))
CBASE1(b_tanh, z_tanh(z))
CBASE1(b_asin, z_asin(z))
CBASE1(b_acos, z_acos(z))
CBASE1(b_atan, z_atan(z))
CBASE1(b_asinh, z_asinh(z))
CBASE1(b_acosh, z_cacosh(z))
CBASE1(b_atanh, z_atanh(z))

static void b_pow(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  from_c99(mycpow(to_c99(x), to_c99(y)), out);
}

/* asin_cut is found by cmath_base_get(). */
static rsimd_cmath_base cmath_base = {
  {b_sqrt, b_exp, b_log, b_sin, b_cos, b_tan, b_sinh, b_cosh, b_tanh, b_asin, b_acos, b_atan,
   b_asinh, b_acosh, b_atanh},
  {b_pow, z_logbase, z_atan2},
  RSIMD_CUT_UNKNOWN};

/* How z_asin, as this file's compiler built it, rounds alpha * alpha - 1
   on its branch cut (rsimd_cmath_base in kernel_types.h): b_asin is
   compared with both roundings on real inputs just above 1, where they
   give different results. */
static int probe_asin_cut(void) {
  int k, fused = 1, unfused = 1;
  for (k = 1; k <= 8; k++) {
    const double x = 1 + ldexp((double) (2 * k + 1), -30);
    const double alpha = 0.5 * fabs(x + 1) + 0.5 * fabs(x - 1);
    volatile double sq = alpha * alpha; /* rounded on its own */
    Rcomplex z, r;
    z.r = x;
    z.i = 0;
    b_asin(&z, &r);
    if (r.i != -log(alpha + sqrt(sq - 1))) unfused = 0;
    if (r.i != -log(alpha + sqrt(fma(alpha, alpha, -1)))) fused = 0;
  }
  return fused == unfused ? RSIMD_CUT_UNKNOWN : fused ? RSIMD_CUT_FUSED : RSIMD_CUT_UNFUSED;
}

static const rsimd_cmath_base *cmath_base_get(void) {
  static int probed = 0;
  if (!probed) {
    cmath_base.asin_cut = probe_asin_cut();
    probed = 1;
  }
  return &cmath_base;
}

/* ---- Entry points --------------------------------------------------------- */

static int lookup_op(SEXP op, const char *const *names, int count) {
  const char *name = rsimd_arg_str(op, "op");
  int i;
  for (i = 0; i < count; i++) {
    if (strcmp(name, names[i]) == 0) return i;
  }
  Rf_error("internal error: unknown op '%s'", name);
  return -1; /* not reached */
}

static int fast_bit(SEXP accuracy) {
  int a = rsimd_arg_int1(accuracy, "accuracy");
  if (a != 0 && a != 1) Rf_error("internal error: invalid accuracy code %d", a);
  return a ? RSIMD_MATH_FAST : 0;
}

/* Base R's names of the functions, for the warning. */
static const char *const cm1_names[] = {"sqrt", "exp",  "log",  "sin",  "cos",
                                        "tan",  "sinh", "cosh", "tanh", "asin",
                                        "acos", "atan", "asinh", "acosh", "atanh"};

/* The function `op` (by name) of complex z, in accuracy mode `accuracy`. */
static SEXP simd_cmath1_impl(SEXP z, SEXP op, SEXP accuracy) {
  int code = lookup_op(op, cm1_names, RSIMD_CM_COUNT), fast = fast_bit(accuracy), nan = 0;
  SEXP out;
  Rcomplex *po;
  rsimd_in in;

  rsimd_in_init(&in, z, "z");
  if (in.type != RSIMD_C128) Rf_error("internal error: complex %s of %s", cm1_names[code],
                                      rsimd_etype_names[in.type]);
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, in.n));
  po = (Rcomplex *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off, {
    nan |= rsimd_active->cmath1_c128(code | fast, px, len, po + off, cmath_base_get());
  });
  if (nan) rsimd_warn("NaNs produced in function \"%s\"", cm1_names[code]);
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cmath1(SEXP z, SEXP op, SEXP accuracy) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_cmath1_impl(z, op, accuracy), 0));
}

/* x^y ("pow"), log(x, base = y) ("logb") or atan2(x, y) ("atan2") of two
   complex vectors with the length-1 broadcast rule; `name` is the
   function's name in base R's warning. */
static SEXP simd_cmath2_impl(SEXP x, SEXP y, SEXP op, SEXP accuracy, SEXP name) {
  static const char *const names[] = {"pow", "logb", "atan2"};
  int code = lookup_op(op, names, RSIMD_CM2_COUNT), fast = fast_bit(accuracy), nan = 0, flags;
  const char *fname = rsimd_arg_str(name, "name");
  SEXP out;
  Rcomplex *po;
  rsimd_bin b;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
    Rf_error("internal error: complex %s of %s and %s", names[code], rsimd_etype_names[b.x.type],
             rsimd_etype_names[b.y.type]);
  }
  flags = (b.x_scalar ? RSIMD_EW_SCALAR(0) : 0) | (b.y_scalar ? RSIMD_EW_SCALAR(1) : 0);
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, b.n));
  po = (Rcomplex *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off, {
    nan |= rsimd_active->cmath2_c128(code | fast, px, py, len, flags, po + off, cmath_base_get(),
                                     rsimd_c128_arith_get());
  });
  if (nan) rsimd_warn("NaNs produced in function \"%s\"", fname);
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cmath2(SEXP x, SEXP y, SEXP op, SEXP accuracy, SEXP name) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_cmath2_impl(x, y, op, accuracy, name), 0));
}
