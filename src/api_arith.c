/* .Call entry points of the elementwise operations: classify the operands,
   apply the length-1 broadcast rule, pick the double or the integer kernel
   from the operand types (as base R's arithmetic does: logical and integer
   give integer, anything with a double gives double), run it chunk by
   chunk through the active implementation, and turn its status bits into
   base R's warnings. Results are bare vectors. */

#include <string.h>
#include <Rmath.h>
#include "rsimd.h"
#include "dispatch.h"
#include "api_complex.h"

static int is_int_like(rsimd_etype t) { return t == RSIMD_I32 || t == RSIMD_LGL; }

static int lookup_op(SEXP op, const char *const *names, int count) {
  const char *name = rsimd_arg_str(op, "op");
  int i;
  for (i = 0; i < count; i++) {
    if (strcmp(name, names[i]) == 0) return i;
  }
  Rf_error("internal error: unknown op '%s'", name);
  return -1; /* not reached */
}

/* Errors for an operand type the elementwise ops do not take, with base
   R's message for raw (`raw_msg`); complex operands are allowed when
   `cplx` is set. Returns 1 if an operand is complex. The R side rejects
   complex for the other ops and integer64, naming the function. */
static int check_numeric(const rsimd_ew *e, const char *raw_msg, int cplx) {
  int i, any_c128 = 0;
  for (i = 0; i < e->k; i++) {
    if (e->in[i].type == RSIMD_U8) Rf_error("%s", raw_msg);
  }
  for (i = 0; i < e->k; i++) {
    rsimd_etype t = e->in[i].type;
    if (t == RSIMD_C128 && cplx) any_c128 = 1;
    else if (t != RSIMD_F64 && !is_int_like(t)) {
      Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[t]);
    }
  }
  return any_c128;
}

/* RSIMD_EW_I32(k) for each int32 operand of a double kernel. */
static int i32_flags(const rsimd_ew *e) {
  int i, f = 0;
  for (i = 0; i < e->k; i++) {
    if (is_int_like(e->in[i].type)) f |= RSIMD_EW_I32(i);
  }
  return f;
}

static int all_int_like(const rsimd_ew *e) {
  int i;
  for (i = 0; i < e->k; i++) {
    if (!is_int_like(e->in[i].type)) return 0;
  }
  return 1;
}

static void warn_status(int st) {
  if (st & RSIMD_EW_OVERFLOW) rsimd_warn_int_overflow();
  if (st & RSIMD_EW_NAN_PRODUCED) Rf_warning("NaNs produced");
}

static const char *const binop_msg = "non-numeric argument to binary operator";

/* Binary ops by name (kernel_types.h): x + y, x - y, x * y, x / y,
   x %/% y, x %% y, pmin, pmax, their na.rm = TRUE forms, copysign and the
   wrapping integer ops. Integer and logical operands give an integer
   result, except for / and copysign; the R side ensures that the _wrap ops
   only see integer and logical operands. + and - also take complex
   operands, which the R side has made both complex. */
SEXP C_simd_ew2(SEXP x, SEXP y, SEXP op, SEXP na_check) {
  static const char *const names[] = {"add",      "sub",      "mul",      "div",     "idiv",
                                      "mod",      "pmin",     "pmax",     "pmin_num", "pmax_num",
                                      "copysign", "add_wrap", "sub_wrap", "mul_wrap"};
  static const char *const args[] = {"x", "y"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0])), st = 0;
  SEXP sargs[2], out;
  rsimd_opts o;
  rsimd_ew e;

  sargs[0] = x;
  sargs[1] = y;
  rsimd_ew_init(&e, 2, sargs, args);
  if (check_numeric(&e, binop_msg, code == RSIMD_EW_ADD || code == RSIMD_EW_SUB)) {
    rsimd_opts_init(&o, R_NilValue, na_check, R_NilValue, 0);
    return rsimd_c128_add(code, x, y, &o);
  }
  rsimd_opts_init(&o, R_NilValue, na_check, R_NilValue, e.no_na_hint);
  if (all_int_like(&e) && code != RSIMD_EW_DIV && code != RSIMD_EW_COPYSIGN) {
    int *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_I32, e.n));
    po = (int *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew2_i32(code, (const int *) p[0], (const int *) p[1], len, e.flags,
                                  po + off, &o);
    });
  } else {
    double *po;
    int flags = e.flags | i32_flags(&e);
    if (code >= RSIMD_EW_ADD_WRAP) Rf_error("internal error: wrapping op on doubles");
    out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
    po = (double *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew2_f64(code, p[0], p[1], len, flags, po + off, &o);
    });
  }
  warn_status(st);
  UNPROTECT(1);
  return out;
}

/* Ternary ops by name: fma(x, y, z), mul_add, add_mul, lerp(x, y, t) and
   clamp(x, lo, hi). mul_add, add_mul and clamp of integer and logical
   operands give an integer result; everything else is double. clamp
   errors when lo > hi anywhere. */
SEXP C_simd_ew3(SEXP x, SEXP y, SEXP z, SEXP op, SEXP na_check) {
  static const char *const names[] = {"fma", "mul_add", "add_mul", "lerp", "clamp"};
  static const char *const xyz[] = {"x", "y", "z"}, *const xyt[] = {"x", "y", "t"},
                           *const xlohi[] = {"x", "lo", "hi"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0])), st = 0;
  SEXP sargs[3], out;
  rsimd_opts o;
  rsimd_ew e;

  sargs[0] = x;
  sargs[1] = y;
  sargs[2] = z;
  rsimd_ew_init(&e, 3, sargs,
                code == RSIMD_EW_CLAMP ? xlohi : code == RSIMD_EW_LERP ? xyt : xyz);
  check_numeric(&e, binop_msg, 0);
  rsimd_opts_init(&o, R_NilValue, na_check, R_NilValue, e.no_na_hint);
  if (all_int_like(&e) && code != RSIMD_EW_FMA && code != RSIMD_EW_LERP) {
    int *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_I32, e.n));
    po = (int *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew3_i32(code, (const int *) p[0], (const int *) p[1],
                                  (const int *) p[2], len, e.flags, po + off, &o);
    });
  } else {
    double *po;
    int flags = e.flags | i32_flags(&e);
    out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
    po = (double *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew3_f64(code, p[0], p[1], p[2], len, flags, po + off, &o);
    });
  }
  if (st & RSIMD_EW_LO_GT_HI) Rf_error("'lo' must not be greater than 'hi'");
  warn_status(st);
  UNPROTECT(1);
  return out;
}

/* Unary ops by name: neg and abs keep integer and logical input integer
   (they wrap, so NA stays NA); sign, recip, sqrt, floor, ceiling, trunc
   and round (half to even) are double for every input, as in base R.
   sqrt warns "NaNs produced" for a negative number. neg also takes
   complex input. */
SEXP C_simd_ew1(SEXP x, SEXP op) {
  static const char *const names[] = {"neg",   "abs",     "sign",  "recip", "sqrt",
                                      "floor", "ceiling", "trunc", "round"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0])), st = 0;
  SEXP out;
  rsimd_opts o;
  rsimd_ew e;

  rsimd_ew_init(&e, 1, &x, (const char *const[]){"x"});
  if (check_numeric(&e,
                    code == RSIMD_EW_NEG ? "invalid argument to unary operator"
                                         : "non-numeric argument to mathematical function",
                    code == RSIMD_EW_NEG)) {
    return rsimd_c128_neg(&e.in[0]);
  }
  rsimd_opts_init(&o, R_NilValue, R_NilValue, R_NilValue, e.no_na_hint);
  if (is_int_like(e.in[0].type) && (code == RSIMD_EW_NEG || code == RSIMD_EW_ABS)) {
    int *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_I32, e.n));
    po = (int *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew1_i32(code, (const int *) p[0], len, po + off, &o);
    });
  } else {
    double *po;
    int flags = i32_flags(&e);
    out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
    po = (double *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      st |= rsimd_active->ew1_f64(code, p[0], len, flags, po + off, &o);
    });
  }
  warn_status(st);
  UNPROTECT(1);
  return out;
}

/* round(x, digits) for digits other than 0, through R's own fround() one
   element at a time (not vectorised), with base R's missing-value rule
   for two-argument math functions: NA if x or digits is NA, else NaN if
   either is NaN. */
SEXP C_simd_round_digits(SEXP x, SEXP digits) {
  double d = rsimd_arg_num1(digits, "digits");
  SEXP out;
  double *po;
  rsimd_ew e;

  rsimd_ew_init(&e, 1, &x, (const char *const[]){"x"});
  check_numeric(&e, "non-numeric argument to mathematical function", 0);
  out = PROTECT(rsimd_alloc_like(RSIMD_F64, e.n));
  po = (double *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
    R_xlen_t i;
    for (i = 0; i < len; i++) {
      double v;
      if (e.in[0].type == RSIMD_F64) {
        v = ((const double *) p[0])[i];
      } else {
        int iv = ((const int *) p[0])[i];
        v = iv == NA_INTEGER ? NA_REAL : (double) iv;
      }
      if (ISNA(v) || ISNA(d)) po[off + i] = NA_REAL;
      else if (ISNAN(v) || ISNAN(d)) po[off + i] = R_NaN;
      else po[off + i] = fround(v, d);
    }
  });
  UNPROTECT(1);
  return out;
}
