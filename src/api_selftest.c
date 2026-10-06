/* Internal .Call routines that run the self-test kernel slots through the
   active implementation, so the tests can check the NA, precision and
   overflow rules of na.h on every tier. Not part of the package's
   interface. */

#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"
#include "selftest.h"

static int lookup_name(const char *name, const char *const *table, int count, const char *what) {
  int i;
  for (i = 0; i < count; i++) {
    if (table[i] != NULL && strcmp(name, table[i]) == 0) return i;
  }
  Rf_error("unknown %s '%s'", what, name);
  return -1; /* not reached */
}

static SEXP finished(SEXP value, const rsimd_reduce_result *r) {
  const char *names[] = {"value", "count", "saw_na", "saw_nan", "any_true", "any_false", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, value);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal((double) r->count));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(r->saw_na));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r->saw_nan));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(r->any_true));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(r->any_false));
  UNPROTECT(1);
  return out;
}

/* A sum-like reduction of x (and y for term "xy") folded chunk by chunk:
   list(value, count, saw_na, saw_nan, any_true, any_false). term is "x"
   (sum), "sq" (sum_sq), "abs" (sum_abs) or "xy" (dot) for doubles; integer
   and logical x take "x" only. precision is the integer code. */
static SEXP simd_debug_fold_impl(SEXP x, SEXP y, SEXP term, SEXP na_rm, SEXP na_check,
                                 SEXP precision) {
  static const char *const terms[] = {"x", "sq", "abs", "xy"};
  static const int ops[] = {RSIMD_RED_SUM, RSIMD_RED_SUM_SQ, RSIMD_RED_SUM_ABS, RSIMD_RED_DOT};
  int t = lookup_name(rsimd_arg_str(term, "term"), terms, 4, "term");
  rsimd_reduce_result r;
  rsimd_opts o;
  SEXP value;

  rsimd_reduce_result_init(&r, ops[t]);
  if (t == RSIMD_TERM_XY) {
    rsimd_bin b;
    rsimd_bin_init(&b, x, y);
    if (b.x.type != RSIMD_F64 || b.y.type != RSIMD_F64 || b.x_scalar || b.y_scalar) {
      Rf_error("'x' and 'y' must be double vectors of equal length");
    }
    rsimd_opts_init(&o, na_rm, na_check, precision, b.x.no_na_hint && b.y.no_na_hint);
    RSIMD_FOREACH_CHUNK2(&b, double, px, py, len, off, {
      rsimd_active->selftest_fold_f64(px, py, len, t, &r, &o);
    });
    value = PROTECT(rsimd_reduce_finish(ops[t], RSIMD_F64, b.n, &r, &o));
  } else {
    rsimd_in in;
    rsimd_in_init(&in, x, "x");
    rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
    if (in.type == RSIMD_F64) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
        rsimd_active->selftest_fold_f64(px, NULL, len, t, &r, &o);
      });
    } else if ((in.type == RSIMD_I32 || in.type == RSIMD_LGL) && t == RSIMD_TERM_X) {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        rsimd_active->selftest_sum_i32(px, len, &r, &o);
      });
    } else {
      Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
    }
    value = PROTECT(rsimd_reduce_finish(ops[t], in.type, in.n, &r, &o));
  }
  value = finished(value, &r);
  UNPROTECT(1);
  return value;
}

SEXP C_simd_debug_fold(SEXP x, SEXP y, SEXP term, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_debug_fold_impl(x, y, term, na_rm, na_check, precision));
}

/* any(x) (op "any", stopping at the first TRUE), all(x) (op "all",
   stopping at the first FALSE) or the flags of a full scan finished as any
   (op "scan"), for logical or integer x. */
static SEXP simd_debug_lgl_impl(SEXP x, SEXP op, SEXP na_rm, SEXP na_check) {
  static const char *const ops[] = {"any", "all", "scan"};
  int k = lookup_name(rsimd_arg_str(op, "op"), ops, 3, "op");
  int stop = k == 0 ? RSIMD_STOP_TRUE : k == 1 ? RSIMD_STOP_FALSE : RSIMD_STOP_NONE;
  int red = k == 1 ? RSIMD_RED_ALL : RSIMD_RED_ANY;
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;
  SEXP value;

  rsimd_in_init(&in, x, "x");
  if (in.type != RSIMD_LGL && in.type != RSIMD_I32) {
    Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  rsimd_opts_init(&o, na_rm, na_check, R_NilValue, in.no_na_hint);
  rsimd_reduce_result_init(&r, red);
  RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
    rsimd_active->selftest_lgl(px, len, stop, &r, &o);
    if ((stop == RSIMD_STOP_TRUE && r.any_true) || (stop == RSIMD_STOP_FALSE && r.any_false)) {
      break;
    }
  });
  value = PROTECT(rsimd_reduce_finish(red, RSIMD_LGL, in.n, &r, &o));
  value = finished(value, &r);
  UNPROTECT(1);
  return value;
}

SEXP C_simd_debug_lgl(SEXP x, SEXP op, SEXP na_rm, SEXP na_check) {
  rsimd_entry();
  return rsimd_exit(simd_debug_lgl_impl(x, op, na_rm, na_check));
}

/* Elementwise op on x and y (same type, equal lengths or one of length 1):
   integer ops "add", "sub", "mul", "neg", "abs", their "_wrap" forms,
   "idiv" and "mod"; double ops "add", "sub", "mul", "div", "pmin",
   "pmax". Unary ops ignore y. Warns once on integer overflow. */
static SEXP simd_debug_arith_impl(SEXP x, SEXP y, SEXP op, SEXP na_check) {
  static const char *const i32_ops[RSIMD_ST_I32_COUNT] = {
    "add", "sub", "mul", "neg", "abs", "add_wrap", "sub_wrap", "mul_wrap", "neg_wrap", "abs_wrap",
    "idiv", "mod"};
  static const char *const f64_ops[RSIMD_ST_F64_COUNT] = {"add", "sub", "mul",
                                                          "div", "pmin", "pmax"};
  const char *name = rsimd_arg_str(op, "op");
  rsimd_bin b;
  rsimd_opts o;
  SEXP out;
  int code, overflow = 0;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != b.y.type) Rf_error("'x' and 'y' must have the same type");
  rsimd_opts_init(&o, R_NilValue, na_check, R_NilValue, b.x.no_na_hint && b.y.no_na_hint);
  out = PROTECT(rsimd_alloc_like(b.x.type == RSIMD_LGL ? RSIMD_I32 : b.x.type, b.n));
  switch (b.x.type) {
  case RSIMD_I32:
  case RSIMD_LGL: {
    int *po = (int *) rsimd_out_ptr(out);
    code = lookup_name(name, i32_ops, RSIMD_ST_I32_COUNT, "integer op");
    RSIMD_FOREACH_CHUNK2(&b, int, px, py, len, off, {
      overflow |= rsimd_active->selftest_arith_i32(code, px, py, len, b.x_scalar, b.y_scalar,
                                                   po + off, &o);
    });
    break;
  }
  case RSIMD_F64: {
    double *po = (double *) rsimd_out_ptr(out);
    code = lookup_name(name, f64_ops, RSIMD_ST_F64_COUNT, "double op");
    RSIMD_FOREACH_CHUNK2(&b, double, px, py, len, off, {
      rsimd_active->selftest_arith_f64(code, px, py, len, b.x_scalar, b.y_scalar, po + off, &o);
    });
    break;
  }
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[b.x.type]);
  }
  if (overflow) rsimd_warn_int_overflow();
  UNPROTECT(1);
  return out;
}

SEXP C_simd_debug_arith(SEXP x, SEXP y, SEXP op, SEXP na_check) {
  rsimd_entry();
  return rsimd_exit(simd_debug_arith_impl(x, y, op, na_check));
}
