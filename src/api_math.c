/* .Call entry points of the elementary functions: classify the operands
   (double, integer or logical, read lane-wise by the kernels), apply the
   length-1 broadcast rule to the two-operand functions, run the kernel
   chunk by chunk through the active implementation and warn "NaNs
   produced" as base R does. Results are bare double vectors. */

#include <math.h>
#include <stdint.h>
#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

static const char *const math_msg = "non-numeric argument to mathematical function";

static int lookup_op(SEXP op, const char *const *names, int count) {
  const char *name = rsimd_arg_str(op, "op");
  int i;
  for (i = 0; i < count; i++) {
    if (strcmp(name, names[i]) == 0) return i;
  }
  Rf_error("internal error: unknown op '%s'", name);
  return -1; /* not reached */
}

/* Errors with `msg` unless every operand is double, integer or logical
   (complex and integer64 operands are rejected on the R side, which names
   the function), and returns the RSIMD_EW_I32 flags of the int32 ones. */
static int check_operands(const rsimd_ew *e, const char *msg) {
  int i, f = 0;
  for (i = 0; i < e->k; i++) {
    rsimd_etype t = e->in[i].type;
    if (t == RSIMD_I32 || t == RSIMD_LGL) {
      f |= RSIMD_EW_I32(i);
    } else if (t != RSIMD_F64) {
      Rf_error("%s", msg);
    }
  }
  return f;
}

static void warn_nan(int st) {
  if (st & RSIMD_EW_NAN_PRODUCED) Rf_warning("NaNs produced");
}

/* What run_math1() does with the kernel's result: keep it and warn as the
   kernel says, keep it without warning, or ignore the kernel and give NA
   everywhere. */
enum { MATH1_WARN, MATH1_QUIET, MATH1_ALL_NA };

/* math1_f64 op `code` (with the LOGB divisor p) over x. */
static SEXP run_math1(SEXP x, int code, double p, int mode) {
  int flags, st = 0;
  SEXP out;
  double *po;
  rsimd_ew e;

  rsimd_ew_init(&e, 1, &x, (const char *const[]){"x"});
  flags = check_operands(&e, math_msg);
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  po = (double *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK_EW(&e, pp, len, off, {
    if (mode == MATH1_ALL_NA) {
      R_xlen_t i;
      for (i = 0; i < len; i++) po[off + i] = NA_REAL;
    } else {
      st |= rsimd_active->math1_f64(code, pp[0], len, flags, p, po + off);
    }
  });
  if (mode == MATH1_WARN) warn_nan(st);
  UNPROTECT(1);
  return out;
}

/* Unary functions by name (the RSIMD_MATH_* order of kernel_types.h,
   without LOGB). */
SEXP C_simd_math1(SEXP x, SEXP op) {
  static const char *const names[] = {
    "exp",  "exp2", "exp10", "expm1", "log",   "log2",  "log10", "log1p", "",
    "cbrt", "sin",  "cos",   "tan",   "asin",  "acos",  "atan",  "sinpi", "cospi",
    "tanpi", "sinh", "cosh", "tanh",  "asinh", "acosh", "atanh"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0]));
  if (code == RSIMD_MATH_LOGB) Rf_error("internal error: unknown op ''");
  return run_math1(x, code, 1.0, MATH1_WARN);
}

/* log(x, base) for a single double, integer or logical base, as base R's
   logbase(): log10 for base
   10, log2 for base 2, otherwise log(x) / log(base) with log(base) NaN for
   a negative base; with base NA every result is NA, with base NaN every
   result is NaN except NA where x is NA, and neither warns. */
SEXP C_simd_log(SEXP x, SEXP base) {
  double b;
  switch (TYPEOF(base)) {
  case REALSXP:
  case INTSXP:
  case LGLSXP: break;
  default: Rf_error("%s", math_msg);
  }
  if (Rf_xlength(base) != 1) Rf_error("'base' must be a single number");
  b = Rf_asReal(base);
  /* log(x) / NaN is NaN, and the input's NaN where x is NA or NaN. */
  if (ISNA(b)) return run_math1(x, RSIMD_MATH_LOGB, b, MATH1_ALL_NA);
  if (ISNAN(b)) return run_math1(x, RSIMD_MATH_LOGB, b, MATH1_QUIET);
  if (b == 10) return run_math1(x, RSIMD_MATH_LOG10, 1.0, MATH1_WARN);
  if (b == 2) return run_math1(x, RSIMD_MATH_LOG2, 1.0, MATH1_WARN);
  return run_math1(x, RSIMD_MATH_LOGB, b > 0 ? log(b) : b == 0 ? R_NegInf : R_NaN, MATH1_WARN);
}

/* Binary functions by name: pow(x, y), atan2(y, x), hypot(x, y). */
SEXP C_simd_math2(SEXP x, SEXP y, SEXP op) {
  static const char *const names[] = {"pow", "atan2", "hypot"};
  static const char *const xy[] = {"x", "y"}, *const yx[] = {"y", "x"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0])), flags, st = 0;
  SEXP sargs[2], out;
  double *po;
  rsimd_ew e;

  sargs[0] = x;
  sargs[1] = y;
  rsimd_ew_init(&e, 2, sargs, code == RSIMD_MATH_ATAN2 ? yx : xy);
  flags = e.flags | check_operands(&e, code == RSIMD_MATH_POW
                                         ? "non-numeric argument to binary operator"
                                         : math_msg);
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  po = (double *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
    st |= rsimd_active->math2_f64(code, p[0], p[1], len, flags, po + off);
  });
  warn_nan(st);
  UNPROTECT(1);
  return out;
}

/* list(sin = sin(x), cos = cos(x)), one warning for both. */
SEXP C_simd_sincos(SEXP x) {
  int flags, st = 0;
  SEXP s, c, out, names;
  double *ps, *pc;
  rsimd_ew e;

  rsimd_ew_init(&e, 1, &x, (const char *const[]){"x"});
  flags = check_operands(&e, math_msg);
  s = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  c = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  ps = (double *) rsimd_out_ptr(s);
  pc = (double *) rsimd_out_ptr(c);
  RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
    st |= rsimd_active->sincos_f64(p[0], len, flags, ps + off, pc + off);
  });
  out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, s);
  SET_VECTOR_ELT(out, 1, c);
  names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("sin"));
  SET_STRING_ELT(names, 1, Rf_mkChar("cos"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  warn_nan(st);
  UNPROTECT(4);
  return out;
}

/* The bits of x as an integer ordered like the doubles: -0 and +0 both
   map to 0, and adjacent doubles differ by 1. */
static int64_t ordered_bits(double x) {
  int64_t i;
  memcpy(&i, &x, sizeof i);
  return i < 0 ? INT64_MIN - i : i;
}

/* The distance in units in the last place between a[i] and b[i], as
   doubles: 0 for two NAs, two other NaNs or the same infinity, Inf when
   only one is missing or one is NA and the other NaN; +0 and -0 are 0
   apart. Used by the accuracy tests. */
SEXP C_simd_ulp_dist(SEXP a, SEXP b) {
  static const char *const names[] = {"a", "b"};
  SEXP sargs[2], out;
  double *po;
  rsimd_ew e;

  sargs[0] = a;
  sargs[1] = b;
  rsimd_ew_init(&e, 2, sargs, names);
  if (e.in[0].type != RSIMD_F64 || e.in[1].type != RSIMD_F64) {
    Rf_error("'a' and 'b' must be double vectors");
  }
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  po = (double *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i++) {
      double x = ((const double *) p[0])[(e.flags & RSIMD_EW_SCALAR(0)) ? 0 : i];
      double y = ((const double *) p[1])[(e.flags & RSIMD_EW_SCALAR(1)) ? 0 : i];
      if (ISNAN(x) || ISNAN(y)) {
        po[off + i] = (ISNAN(x) && ISNAN(y) && ISNA(x) == ISNA(y)) ? 0 : R_PosInf;
      } else if (x == y) {
        po[off + i] = 0;
      } else {
        /* Below 2^64 for any two non-NaN doubles, so exact in uint64_t. */
        int64_t ox = ordered_bits(x), oy = ordered_bits(y);
        po[off + i] =
          (double) (ox > oy ? (uint64_t) ox - (uint64_t) oy : (uint64_t) oy - (uint64_t) ox);
      }
    }
  });
  UNPROTECT(1);
  return out;
}
