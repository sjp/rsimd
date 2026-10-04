/* R vector access layer; see rvec.h. This is the only file that uses R's
   raw vector accessors and R_CheckUserInterrupt (tools/lint_c.sh). */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "rvec.h"
#include "dispatch.h"
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
                           : (TYPEOF(x) == VECSXP ? "list" : Rf_type2char((SEXPTYPE) TYPEOF(x)));
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

/* ---- simd_vec operands ---------------------------------------------- */

/* What the current .Call has seen of its simd_vec operands. */
static struct {
  int seen;        /* an operand is a simd_vec */
  rsimd_tier pin;  /* their common pinned tier, or RSIMD_TIER_COUNT */
  int na_free;     /* every operand is known NA-free */
} sv_call = {0, RSIMD_TIER_COUNT, 1};

void rsimd_entry(void) {
  rsimd_active = rsimd_selected;
  sv_call.seen = 0;
  sv_call.pin = RSIMD_TIER_COUNT;
  sv_call.na_free = 1;
}

/* 1 if the single element of v (contiguous, n == 1) is missing. */
static int scalar_is_na(const rsimd_in *v) {
  switch (v->type) {
  case RSIMD_F64: return ISNAN(*(const double *) v->ptr);
  case RSIMD_I32:
  case RSIMD_LGL: return *(const int *) v->ptr == NA_INTEGER;
  case RSIMD_I64: {
    int64_t b;
    memcpy(&b, v->ptr, sizeof b);
    return b == INT64_MIN;
  }
  case RSIMD_C128: {
    const Rcomplex *z = (const Rcomplex *) v->ptr;
    return ISNAN(z->r) || ISNAN(z->i);
  }
  default: return 0;
  }
}

/* Records operand v (named arg) in sv_call: a simd_vec's rsimd_na_free
   flag sets v->no_na_hint, and its rsimd_impl pin switches rsimd_active to
   the pinned tier's table, erroring if the tier is not available or two
   operands are pinned to different tiers. */
static void sv_note(rsimd_in *v, const char *arg) {
  SEXP x = v->sx, flag, impl;
  if (Rf_inherits(x, "simd_vec")) {
    sv_call.seen = 1;
    flag = Rf_getAttrib(x, Rf_install("rsimd_na_free"));
    if (TYPEOF(flag) == LGLSXP && XLENGTH(flag) == 1 && LOGICAL_ELT(flag, 0) == TRUE) {
      v->no_na_hint = 1;
    }
    impl = Rf_getAttrib(x, Rf_install("rsimd_impl"));
    if (impl != R_NilValue) {
      rsimd_tier t = RSIMD_TIER_COUNT;
      if (TYPEOF(impl) == STRSXP && XLENGTH(impl) == 1 && STRING_ELT(impl, 0) != NA_STRING) {
        t = rsimd_tier_from_name(CHAR(STRING_ELT(impl, 0)));
      }
      if (t == RSIMD_TIER_COUNT) {
        Rf_errorcall(R_NilValue,
                     "'%s' has an invalid 'rsimd_impl' attribute; reset it with "
                     "simd_impl(%s) <- NULL",
                     arg, arg);
      }
      if (sv_call.pin != RSIMD_TIER_COUNT && sv_call.pin != t) {
        Rf_errorcall(R_NilValue,
                     "operands pinned to different implementations ('%s' vs '%s'); "
                     "unpin one with simd_impl(x) <- NULL",
                     rsimd_tier_names[sv_call.pin], rsimd_tier_names[t]);
      }
      if (rsimd_tier_resolved(t) == NULL) {
        Rf_errorcall(R_NilValue,
                     "'%s' is pinned to implementation '%s', which is not available on "
                     "this machine; unpin it with simd_impl(%s) <- NULL",
                     arg, rsimd_tier_names[t], arg);
      }
      sv_call.pin = t;
      rsimd_active = rsimd_tier_resolved(t);
    }
  }
  sv_call.na_free = sv_call.na_free && v->no_na_hint;
}

SEXP rsimd_sv_result(SEXP out, int keeps_na_free) {
  SEXP cls, impl, flag, impl_sym, flag_sym;
  int na_free;
  if (!sv_call.seen) return out;
  if (MAYBE_REFERENCED(out)) out = Rf_shallow_duplicate(out);
  PROTECT(out);
  if (Rf_inherits(out, "integer64")) {
    cls = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(cls, 0, Rf_mkChar("simd_vec"));
    SET_STRING_ELT(cls, 1, Rf_mkChar("integer64"));
  } else {
    cls = PROTECT(Rf_mkString("simd_vec"));
  }
  Rf_setAttrib(out, R_ClassSymbol, cls);
  impl_sym = Rf_install("rsimd_impl");
  flag_sym = Rf_install("rsimd_na_free");
  impl = PROTECT(sv_call.pin == RSIMD_TIER_COUNT ? R_NilValue
                                                 : Rf_mkString(rsimd_tier_names[sv_call.pin]));
  Rf_setAttrib(out, impl_sym, impl);
  na_free = TYPEOF(out) == RAWSXP || (keeps_na_free && sv_call.na_free);
  flag = PROTECT(na_free ? Rf_ScalarLogical(TRUE) : R_NilValue);
  Rf_setAttrib(out, flag_sym, flag);
  UNPROTECT(4);
  return out;
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
  if (!v->no_na_hint && v->n == 1 && v->ptr != NULL) v->no_na_hint = !scalar_is_na(v);
  sv_note(v, arg);
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

/* ---- Elementwise operands ----------------------------------------------- */

int rsimd_ew_init(rsimd_ew *e, int k, const SEXP *args, const char *const *names) {
  R_xlen_t n = 1;
  int i, have_n = 0, ok = 1;
  if (k < 1 || k > RSIMD_EW_MAX_ARGS) Rf_error("internal error: %d operands", k);
  e->k = k;
  e->flags = 0;
  e->no_na_hint = 1;
  for (i = 0; i < k; i++) {
    rsimd_in_init(&e->in[i], args[i], names[i]);
    e->no_na_hint = e->no_na_hint && e->in[i].no_na_hint;
    if (e->in[i].n == 1) continue;
    if (!have_n) {
      n = e->in[i].n;
      have_n = 1;
    } else if (e->in[i].n != n) {
      ok = 0;
    }
  }
  if (!ok) {
    if (k == 2) {
      Rf_error("lengths of '%s' (%lld) and '%s' (%lld) must be equal or one of them must be 1",
               names[0], (long long) e->in[0].n, names[1], (long long) e->in[1].n);
    }
    Rf_error("lengths of '%s' (%lld), '%s' (%lld) and '%s' (%lld) must be equal or 1", names[0],
             (long long) e->in[0].n, names[1], (long long) e->in[1].n, names[2],
             (long long) e->in[2].n);
  }
  e->n = n;
  for (i = 0; i < k; i++) {
    if (e->in[i].n == 1 && n != 1) e->flags |= RSIMD_EW_SCALAR(i);
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
  /* min/max kernels compute both extrema: the maximum starts here. */
  r->f64_hi = R_NegInf;
  r->i64_hi = INT64_MIN;
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

void rsimd_set_i64_class(SEXP out) {
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("integer64"));
}

SEXP rsimd_scalar_i64(int64_t v) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  memcpy(REAL(out), &v, sizeof v);
  rsimd_set_i64_class(out);
  UNPROTECT(1);
  return out;
}

void rsimd_warn_i64_overflow(void) {
  Rf_warning("NAs produced by integer64 overflow");
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
    if (type == RSIMD_I64) return rsimd_scalar_i64(RSIMD_NA_I64);
  }
  return Rf_ScalarReal(NA_REAL);
}

SEXP rsimd_range_pair(SEXP lo, SEXP hi) {
  SEXP out;
  if (TYPEOF(lo) == INTSXP && TYPEOF(hi) == INTSXP) {
    out = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(out)[0] = INTEGER(lo)[0];
    INTEGER(out)[1] = INTEGER(hi)[0];
  } else {
    out = PROTECT(Rf_allocVector(REALSXP, 2));
    REAL(out)[0] = Rf_asReal(lo);
    REAL(out)[1] = Rf_asReal(hi);
  }
  UNPROTECT(1);
  return out;
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
      if (r->saw_na || is_int || type == RSIMD_I64 || op == RSIMD_RED_VAR || op == RSIMD_RED_SD) {
        return na_result(op, type);
      }
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
    if (type == RSIMD_I64) {
      if (r->carry != 0 || r->i64 == RSIMD_NA_I64) {
        rsimd_warn_i64_overflow();
        return rsimd_scalar_i64(RSIMD_NA_I64);
      }
      return rsimd_scalar_i64(r->i64);
    }
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
    if (type == RSIMD_I64 && r->count == 0) {
      int is_min = op == RSIMD_RED_MIN;
      Rf_warning("no non-NA value, returning the %s possible integer64 value %s9223372036854775807",
                 is_min ? "highest" : "lowest", is_min ? "+" : "-");
      return rsimd_scalar_i64(is_min ? INT64_MAX : -INT64_MAX);
    }
    if (type == RSIMD_F64) return Rf_ScalarReal(value);
    if (is_int) return Rf_ScalarInteger((int) r->i64);
    if (type == RSIMD_I64) return rsimd_scalar_i64(r->i64);
    break;
  case RSIMD_RED_WHICH_MIN:
  case RSIMD_RED_WHICH_MAX:
    if (is_real || type == RSIMD_I64 || type == RSIMD_U8) {
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
  case RSIMD_RED_COUNT_NA: return Rf_ScalarReal((double) r->i64);
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
  default:
    Rf_error("internal error: no data pointer for type %s", Rf_type2char((SEXPTYPE) TYPEOF(out)));
  }
  return NULL; /* not reached */
}

void rsimd_copy_class(SEXP x, SEXP out) {
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (TYPEOF(cls) != STRSXP || XLENGTH(cls) != 1) return;
  if (strcmp(CHAR(STRING_ELT(cls, 0)), "integer64") == 0 && TYPEOF(out) == REALSXP) {
    Rf_setAttrib(out, R_ClassSymbol, cls);
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

double rsimd_arg_num1(SEXP x, const char *name) {
  if (Rf_xlength(x) == 1) {
    if (TYPEOF(x) == REALSXP) return REAL_ELT(x, 0);
    if (TYPEOF(x) == INTSXP) return INTEGER_ELT(x, 0) == NA_INTEGER ? NA_REAL : INTEGER_ELT(x, 0);
    if (TYPEOF(x) == LGLSXP && LOGICAL_ELT(x, 0) == NA_LOGICAL) return NA_REAL;
  }
  Rf_error("'%s' must be a single number", name);
  return 0; /* not reached */
}

int rsimd_str_in(SEXP x, const char *const *names) {
  const char *s;
  if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || STRING_ELT(x, 0) == NA_STRING) return 0;
  s = CHAR(STRING_ELT(x, 0));
  for (; *names != NULL; names++) {
    if (strcmp(s, *names) == 0) return 1;
  }
  return 0;
}

const char *rsimd_arg_str(SEXP x, const char *name) {
  if (TYPEOF(x) != STRSXP || Rf_xlength(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("'%s' must be a single string", name);
  }
  return CHAR(STRING_ELT(x, 0));
}
