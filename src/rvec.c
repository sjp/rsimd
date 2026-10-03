/* R vector access layer; see rvec.h. This is the only file that uses R's
   raw vector accessors and R_CheckUserInterrupt (tools/lint_c.sh). */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "rvec.h"
#include <R_ext/Utils.h>

_Static_assert(sizeof(Rcomplex) == 2 * sizeof(double), "Rcomplex must be two doubles");

const char *const rsimd_etype_names[RSIMD_BAD + 1] = {
  "double", "integer", "logical", "raw", "integer64", "complex", "unsupported"
};

rsimd_etype rsimd_etype_of(SEXP x) {
  if (Rf_isS4(x)) return RSIMD_BAD;
  switch (TYPEOF(x)) {
  case REALSXP: return Rf_inherits(x, "integer64") ? RSIMD_I64 : RSIMD_F64;
  case INTSXP: return Rf_isFactor(x) ? RSIMD_BAD : RSIMD_I32;
  case LGLSXP: return RSIMD_LGL;
  case RAWSXP: return RSIMD_U8;
  case CPLXSXP: return RSIMD_C128;
  default: return RSIMD_BAD;
  }
}

rsimd_etype rsimd_check_atomic(SEXP x, const char *arg) {
  rsimd_etype t = rsimd_etype_of(x);
  if (t == RSIMD_BAD) {
    /* The first class, or the type for an object without one (R_data_class
       is not part of R's API). */
    SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
    const char *name = TYPEOF(cls) == STRSXP && XLENGTH(cls) > 0
                           ? Rf_translateChar(STRING_ELT(cls, 0))
                           : (TYPEOF(x) == VECSXP ? "list" : Rf_type2char(TYPEOF(x)));
    Rf_error("'%s' must be an atomic vector (double, integer, logical, raw, complex or "
             "integer64), not %s",
             arg, name);
  }
  return t;
}

size_t rsimd_etype_size(rsimd_etype t) {
  switch (t) {
  case RSIMD_F64:
  case RSIMD_I64: return sizeof(double);
  case RSIMD_I32:
  case RSIMD_LGL: return sizeof(int);
  case RSIMD_U8: return sizeof(Rbyte);
  case RSIMD_C128: return sizeof(Rcomplex);
  default: return 0;
  }
}

/* ---- Promotion ---------------------------------------------------------- */

rsimd_etype rsimd_promote(rsimd_etype a, rsimd_etype b) {
  int has_i64 = a == RSIMD_I64 || b == RSIMD_I64;
  if (a == RSIMD_BAD || b == RSIMD_BAD) return RSIMD_BAD;
  if (a == b) return a;
  if (a == RSIMD_U8 || b == RSIMD_U8) return RSIMD_BAD;
  if (a == RSIMD_C128 || b == RSIMD_C128) return has_i64 ? RSIMD_BAD : RSIMD_C128;
  if (a == RSIMD_F64 || b == RSIMD_F64) return RSIMD_F64;
  if (has_i64) return RSIMD_I64;
  return RSIMD_I32; /* logical with integer */
}

int rsimd_promote_warns(rsimd_etype a, rsimd_etype b) {
  return (a == RSIMD_I64 && b == RSIMD_F64) || (a == RSIMD_F64 && b == RSIMD_I64);
}

/* ---- Chunked read access ------------------------------------------------ */

int rsimd_in_init(rsimd_in *v, SEXP x, const char *arg) {
  v->sx = x;
  v->type = rsimd_check_atomic(x, arg);
  v->n = XLENGTH(x);
  /* Never DATAPTR(): it would materialise a compact sequence. */
  v->ptr = DATAPTR_OR_NULL(x);
  switch (v->type) {
  case RSIMD_F64: v->no_na_hint = REAL_NO_NA(x); break;
  case RSIMD_I32: v->no_na_hint = INTEGER_NO_NA(x); break;
  case RSIMD_LGL: v->no_na_hint = LOGICAL_NO_NA(x); break;
  case RSIMD_U8: v->no_na_hint = 1; break; /* raw has no NA */
  default: v->no_na_hint = 0; break;       /* REAL_NO_NA says nothing about int64 NA */
  }
  return 0;
}

const void *rsimd_in_region(const rsimd_in *v, R_xlen_t i, R_xlen_t *len, void *buf) {
  R_xlen_t want, got = 0;
  if (v->ptr != NULL) {
    *len = v->n - i;
    return (const char *) v->ptr + (size_t) i * rsimd_etype_size(v->type);
  }
  want = v->n - i < RSIMD_CHUNK ? v->n - i : RSIMD_CHUNK;
  switch (v->type) {
  case RSIMD_F64:
  case RSIMD_I64: got = REAL_GET_REGION(v->sx, i, want, (double *) buf); break;
  case RSIMD_I32: got = INTEGER_GET_REGION(v->sx, i, want, (int *) buf); break;
  case RSIMD_LGL: got = LOGICAL_GET_REGION(v->sx, i, want, (int *) buf); break;
  case RSIMD_U8: got = RAW_GET_REGION(v->sx, i, want, (Rbyte *) buf); break;
  case RSIMD_C128: got = COMPLEX_GET_REGION(v->sx, i, want, (Rcomplex *) buf); break;
  default: break;
  }
  if (got <= 0) Rf_error("internal error: no elements read at index %.0f", (double) i);
  *len = got;
  return buf;
}

/* ---- Chunk loop --------------------------------------------------------- */

R_xlen_t rsimd_stride = RSIMD_INTERRUPT_STRIDE;

void rsimd_rvec_init(void) {
  /* Invalid values are ignored here; .onLoad warns about them. */
  const char *s = getenv("RSIMD_DEBUG_STRIDE"), *p;
  long long v;
  if (s == NULL || *s == '\0') return;
  for (p = s; *p != '\0'; p++) {
    if (*p < '0' || *p > '9') return;
  }
  v = strtoll(s, NULL, 10); /* LLONG_MAX on overflow: no interrupt checks */
  if (v < 1) return;
  if (v > LLONG_MAX - RSIMD_PAIRWISE_LEAF) v = LLONG_MAX - RSIMD_PAIRWISE_LEAF;
  v = (v + RSIMD_PAIRWISE_LEAF - 1) / RSIMD_PAIRWISE_LEAF * RSIMD_PAIRWISE_LEAF;
  rsimd_stride = (R_xlen_t) v;
}

void rsimd_check_interrupt(void) {
  R_CheckUserInterrupt();
}

/* ---- Binary operands ---------------------------------------------------- */

int rsimd_bin_init(rsimd_bin *b, SEXP x, SEXP y) {
  rsimd_in_init(&b->x, x, "x");
  rsimd_in_init(&b->y, y, "y");
  b->x_scalar = b->y_scalar = 0;
  if (b->x.n == b->y.n) {
    b->n = b->x.n;
  } else if (b->y.n == 1) {
    b->y_scalar = 1;
    b->n = b->x.n;
  } else if (b->x.n == 1) {
    b->x_scalar = 1;
    b->n = b->y.n;
  } else {
    Rf_error("lengths of 'x' (%lld) and 'y' (%lld) must be equal or one of them must be 1",
             (long long) b->x.n, (long long) b->y.n);
  }
  return 0;
}

/* ---- Reductions --------------------------------------------------------- */

const char *const rsimd_reduce_op_names[RSIMD_RED_OP_COUNT] = {
  "sum", "prod", "mean", "min", "max", "which_min", "which_max", "any", "all", "any_na",
  "count_na", "sum_sq", "sum_abs", "dot", "norm", "dist", "cosine", "var", "sd"
};

void rsimd_reduce_result_init(rsimd_reduce_result *r, int op) {
  memset(r, 0, sizeof *r);
  r->idx = -1;
  switch (op) {
  case RSIMD_RED_PROD: r->f64 = 1.0; break;
  case RSIMD_RED_MIN:
  case RSIMD_RED_WHICH_MIN:
    r->f64 = R_PosInf;
    r->i64 = INT64_MAX;
    break;
  case RSIMD_RED_MAX:
  case RSIMD_RED_WHICH_MAX:
    r->f64 = R_NegInf;
    r->i64 = INT64_MIN;
    break;
  default: break;
  }
}

static SEXP scalar_i64(int64_t v) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  memcpy(REAL(out), &v, sizeof v);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("integer64"));
  UNPROTECT(1);
  return out;
}

/* An index or count: integer, or double for a long input. */
static SEXP scalar_index(double v, R_xlen_t n) {
  return n > INT_MAX ? Rf_ScalarReal(v) : Rf_ScalarInteger((int) v);
}

/* NA of the type a reduction returns. */
static SEXP na_result(int op, rsimd_etype type) {
  int is_int = type == RSIMD_I32 || type == RSIMD_LGL;
  if (op == RSIMD_RED_SUM || op == RSIMD_RED_MIN || op == RSIMD_RED_MAX) {
    if (is_int) return Rf_ScalarInteger(NA_INTEGER);
    if (type == RSIMD_I64) return scalar_i64(RSIMD_NA_I64);
  }
  return Rf_ScalarReal(NA_REAL);
}

SEXP rsimd_reduce_finish(int op, rsimd_etype type, R_xlen_t n, const rsimd_reduce_result *r,
                         const rsimd_opts *o) {
  int is_int = type == RSIMD_I32 || type == RSIMD_LGL;
  int is_real = type == RSIMD_F64 || is_int;
  /* Missing values decide the result of the NaN-propagating ops. */
  int missing = !o->na_rm && (r->saw_na || r->saw_nan);
  double value = rsimd_reduce_value(r, o->precision);

  switch (op) {
  case RSIMD_RED_SUM:
  case RSIMD_RED_PROD:
  case RSIMD_RED_MEAN:
  case RSIMD_RED_MIN:
  case RSIMD_RED_MAX:
  case RSIMD_RED_SUM_SQ:
  case RSIMD_RED_SUM_ABS:
  case RSIMD_RED_DOT:
  case RSIMD_RED_NORM:
  case RSIMD_RED_DIST:
  case RSIMD_RED_COSINE:
  case RSIMD_RED_VAR:
  case RSIMD_RED_SD:
    if (!(is_real || type == RSIMD_I64)) break;
    if (type == RSIMD_I64 && op != RSIMD_RED_SUM && op != RSIMD_RED_MIN && op != RSIMD_RED_MAX) {
      break;
    }
    if (missing) {
      /* NA wins over NaN, whatever the order; integer inputs have no NaN. */
      if (r->saw_na || is_int || type == RSIMD_I64) return na_result(op, type);
      return Rf_ScalarReal(R_NaN);
    }
    break;
  default: break;
  }

  switch (op) {
  case RSIMD_RED_SUM:
    if (type == RSIMD_F64) return Rf_ScalarReal(value);
    if (is_int) {
      /* Base R returns a double when an integer sum overflows. */
      if (r->overflow) return Rf_ScalarReal(r->f64 + (double) r->i64);
      if (r->i64 > INT_MAX || r->i64 < -INT_MAX) return Rf_ScalarReal((double) r->i64);
      return Rf_ScalarInteger((int) r->i64);
    }
    if (type == RSIMD_I64) return scalar_i64(r->i64);
    break;
  case RSIMD_RED_MEAN:
    if (is_real) return Rf_ScalarReal(r->count == 0 ? R_NaN : value);
    break;
  case RSIMD_RED_VAR:
  case RSIMD_RED_SD:
    if (is_real) return Rf_ScalarReal(r->count < 2 ? NA_REAL : value);
    break;
  case RSIMD_RED_PROD:
  case RSIMD_RED_SUM_SQ:
  case RSIMD_RED_SUM_ABS:
  case RSIMD_RED_DOT:
  case RSIMD_RED_NORM:
  case RSIMD_RED_DIST:
  case RSIMD_RED_COSINE:
    if (is_real) return Rf_ScalarReal(value);
    break;
  case RSIMD_RED_MIN:
  case RSIMD_RED_MAX:
    if (is_real && r->count == 0) {
      int is_min = op == RSIMD_RED_MIN;
      Rf_warning("no non-missing arguments to %s; returning %s", is_min ? "min" : "max",
                 is_min ? "Inf" : "-Inf");
      return Rf_ScalarReal(is_min ? R_PosInf : R_NegInf);
    }
    if (type == RSIMD_F64) return Rf_ScalarReal(value);
    if (is_int) return Rf_ScalarInteger((int) r->i64);
    if (type == RSIMD_I64) return scalar_i64(r->i64);
    break;
  case RSIMD_RED_WHICH_MIN:
  case RSIMD_RED_WHICH_MAX:
    if (is_real || type == RSIMD_I64) {
      if (r->idx < 0) return Rf_allocVector(n > INT_MAX ? REALSXP : INTSXP, 0);
      return scalar_index((double) r->idx + 1, n);
    }
    break;
  case RSIMD_RED_ANY:
    if (r->any_true) return Rf_ScalarLogical(TRUE);
    return Rf_ScalarLogical(!o->na_rm && r->saw_na ? NA_LOGICAL : FALSE);
  case RSIMD_RED_ALL:
    if (r->any_false) return Rf_ScalarLogical(FALSE);
    return Rf_ScalarLogical(!o->na_rm && r->saw_na ? NA_LOGICAL : TRUE);
  case RSIMD_RED_ANY_NA: return Rf_ScalarLogical(r->saw_na || r->saw_nan);
  case RSIMD_RED_COUNT_NA: return scalar_index((double) r->i64, n);
  default: Rf_error("internal error: unknown reduction %d", op);
  }
  Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[type]);
  return R_NilValue; /* not reached */
}

void rsimd_warn_int_overflow(void) {
  Rf_warning("NAs produced by integer overflow");
}

/* ---- Results ------------------------------------------------------------ */

SEXP rsimd_alloc_like(rsimd_etype t, R_xlen_t n) {
  static const SEXPTYPE sexptype[RSIMD_BAD] = {REALSXP, INTSXP, LGLSXP, RAWSXP, REALSXP, CPLXSXP};
  if (t < 0 || t >= RSIMD_BAD) Rf_error("internal error: cannot allocate type %d", (int) t);
  return Rf_allocVector(sexptype[t], n);
}

void *rsimd_out_ptr(SEXP out) {
  /* A fresh result is never ALTREP, so this does not materialise anything.
     (DATAPTR itself is not part of R's API.) */
  switch (TYPEOF(out)) {
  case REALSXP: return REAL(out);
  case INTSXP: return INTEGER(out);
  case LGLSXP: return LOGICAL(out);
  case RAWSXP: return RAW(out);
  case CPLXSXP: return COMPLEX(out);
  default: Rf_error("internal error: no data pointer for type %s", Rf_type2char(TYPEOF(out)));
  }
  return NULL; /* not reached */
}

void rsimd_copy_class(SEXP x, SEXP out, int no_na) {
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  const char *name;
  if (TYPEOF(cls) != STRSXP || XLENGTH(cls) != 1) return;
  name = CHAR(STRING_ELT(cls, 0));
  if (strcmp(name, "integer64") == 0) {
    if (TYPEOF(out) == REALSXP) Rf_setAttrib(out, R_ClassSymbol, cls);
  } else if (strcmp(name, "simd_vec") == 0) {
    SEXP impl = Rf_install("rsimd_impl");
    Rf_setAttrib(out, R_ClassSymbol, cls);
    Rf_setAttrib(out, impl, Rf_getAttrib(x, impl));
    if (no_na) {
      SEXP sym = Rf_install("rsimd_no_na");
      Rf_setAttrib(out, sym, Rf_getAttrib(x, sym));
    }
  }
}

/* ---- Options and arguments ---------------------------------------------- */

void rsimd_opts_init(rsimd_opts *o, SEXP na_rm, SEXP na_check, SEXP precision, int no_na_hint) {
  o->na_rm = Rf_isNull(na_rm) ? 0 : rsimd_arg_lgl1(na_rm, "na.rm");
  o->na_check = (Rf_isNull(na_check) ? 1 : rsimd_arg_lgl1(na_check, "na_check")) && !no_na_hint;
  o->precision = Rf_isNull(precision) ? RSIMD_PREC_FAST : rsimd_arg_int1(precision, "precision");
  if (o->precision < RSIMD_PREC_FAST || o->precision > RSIMD_PREC_COMPENSATED) {
    Rf_error("internal error: invalid precision code %d", o->precision);
  }
}

int rsimd_arg_lgl1(SEXP x, const char *name) {
  if (TYPEOF(x) != LGLSXP || Rf_xlength(x) != 1 || LOGICAL_ELT(x, 0) == NA_LOGICAL) {
    Rf_error("'%s' must be TRUE or FALSE", name);
  }
  return LOGICAL_ELT(x, 0);
}

int rsimd_arg_int1(SEXP x, const char *name) {
  if (Rf_xlength(x) == 1) {
    if (TYPEOF(x) == INTSXP && INTEGER_ELT(x, 0) != NA_INTEGER) return INTEGER_ELT(x, 0);
    if (TYPEOF(x) == REALSXP) {
      double v = REAL_ELT(x, 0);
      if (v > INT_MIN && v <= INT_MAX && v == (double) (int) v) return (int) v;
    }
  }
  Rf_error("'%s' must be a single integer", name);
  return 0; /* not reached */
}

double rsimd_arg_dbl1(SEXP x, const char *name) {
  if (Rf_xlength(x) == 1) {
    if (TYPEOF(x) == REALSXP && !ISNAN(REAL_ELT(x, 0))) return REAL_ELT(x, 0);
    if (TYPEOF(x) == INTSXP && INTEGER_ELT(x, 0) != NA_INTEGER) return INTEGER_ELT(x, 0);
  }
  Rf_error("'%s' must be a single number", name);
  return 0; /* not reached */
}

const char *rsimd_arg_str(SEXP x, const char *name) {
  if (TYPEOF(x) != STRSXP || Rf_xlength(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("'%s' must be a single string", name);
  }
  return CHAR(STRING_ELT(x, 0));
}
