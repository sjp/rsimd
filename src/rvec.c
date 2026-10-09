/* R vector access layer; see rvec.h. This is the only file that uses R's
   raw vector accessors and R_CheckUserInterrupt (tools/lint_c.sh). */

#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rvec.h"
#include "dispatch.h"
#include <R_ext/Utils.h>
#include <Rversion.h>

_Static_assert(sizeof(Rcomplex) == 2 * sizeof(double), "Rcomplex must be two doubles");

const char *const rsimd_etype_names[RSIMD_BAD + 1] = {
  "double", "integer", "logical", "raw", "integer64", "complex", "unsupported"
};

rsimd_etype rsimd_etype_of(SEXP x) {
  if (Rf_isS4(x)) return RSIMD_BAD;
  /* Classed data (Date, POSIXct, difftime, factor ...) would lose its class
     and meaning; only integer64 and simd_vec are taken. Rf_isObject() only
     tests a bit, so plain input pays nothing for this. */
  if (Rf_isObject(x) && !Rf_inherits(x, "integer64") && !Rf_inherits(x, "simd_vec")) return RSIMD_BAD;
  switch (TYPEOF(x)) {
  case REALSXP: return Rf_inherits(x, "integer64") ? RSIMD_I64 : RSIMD_F64;
  case INTSXP: return RSIMD_I32;
  case LGLSXP: return RSIMD_LGL;
  case RAWSXP: return RSIMD_U8;
  case CPLXSXP: return RSIMD_C128;
  default: return RSIMD_BAD;
  }
}

rsimd_etype rsimd_check_atomic(SEXP x, const char *arg) {
  rsimd_etype t = rsimd_etype_of(x);
  if (t == RSIMD_BAD) {
    /* NULL is rejected, not taken as a zero-length vector. */
    if (Rf_isNull(x)) Rf_error("'%s' is NULL", arg);
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

/* ---- simd_vec operands ---------------------------------------------- */

/* What the current .Call has seen of its simd_vec operands. */
static struct {
  int seen;        /* an operand is a simd_vec */
  rsimd_tier pin;  /* their common pinned tier, or RSIMD_TIER_COUNT */
  int na_free;     /* every operand is known NA-free */
} sv_call = {0, RSIMD_TIER_COUNT, 1};

/* The current .Call's deferred warnings. */
#define RSIMD_WARN_MAX 8
#define RSIMD_WARN_LEN 256
static struct {
  int n;
  char msg[RSIMD_WARN_MAX][RSIMD_WARN_LEN];
} warn_pending;

void rsimd_warn(const char *fmt, ...) {
  va_list ap;
  if (warn_pending.n == RSIMD_WARN_MAX) return;
  va_start(ap, fmt);
  vsnprintf(warn_pending.msg[warn_pending.n++], RSIMD_WARN_LEN, fmt, ap);
  va_end(ap);
}

void rsimd_warn_flush(void) {
  char msg[RSIMD_WARN_MAX][RSIMD_WARN_LEN];
  int i, n = warn_pending.n;
  /* Taken off the list first: a handler may call rsimd again. */
  memcpy(msg, warn_pending.msg, (size_t) n * sizeof msg[0]);
  warn_pending.n = 0;
  for (i = 0; i < n; i++) Rf_warning("%s", msg[i]);
}

SEXP rsimd_exit(SEXP out) {
  if (warn_pending.n == 0) return out;
  PROTECT(out);
  rsimd_warn_flush();
  UNPROTECT(1);
  return out;
}

/* ---- Options ------------------------------------------------------------ */

/* Option and attribute symbols, installed by rsimd_rvec_init(). */
static SEXP sym_impl, sym_precision, sym_math_accuracy, sym_na_check;
static SEXP sym_sv_impl, sym_sv_na_token;

/* The class attributes of a simd_vec result, and the values of its
   rsimd_impl attribute by tier: made once (preserved) and shared by every
   result. rsimd_sv_stamp() gives an object it flags a class vector of its
   own. */
static SEXP sv_cls, sv_cls_i64, sv_impl_value[RSIMD_TIER_COUNT];

/* The value of option rsimd.impl that rsimd_entry() last synced with: its
   CHARSXP (preserved, so that its address is not reused), or R_NilValue
   when the option was unset. impl_synced is 0 before the first sync and
   after every selection (rsimd_impl_selected()). */
static SEXP impl_seen = NULL;
static int impl_synced = 0, impl_syncing = 0;

/* The CHARSXP of option rsimd.impl, R_NilValue when it is unset, or
   R_UnboundValue for any other value (which never counts as synced). */
static SEXP impl_option(void) {
  SEXP v = Rf_GetOption1(sym_impl);
  if (TYPEOF(v) == STRSXP && XLENGTH(v) == 1) return STRING_ELT(v, 0);
  return v == R_NilValue ? R_NilValue : R_UnboundValue;
}

/* Evaluates call (protected by the caller) in the namespace. */
static SEXP ns_eval(SEXP call) {
  SEXP ns = PROTECT(R_FindNamespace(PROTECT(Rf_mkString("rsimd")))), out;
  out = Rf_eval(call, ns);
  UNPROTECT(2);
  return out;
}

static SEXP impl_sync_eval(void *data) {
  SEXP call = PROTECT(Rf_lang1(Rf_install(".sync_option"))), problem;
  (void) data;
  problem = ns_eval(call);
  UNPROTECT(1);
  return problem;
}

static void impl_sync_done(void *data, Rboolean jump) {
  (void) data;
  (void) jump;
  impl_syncing = 0;
}

/* Honours a change to option rsimd.impl made without simd_use(): when the
   option differs from the value last synced, evaluates .sync_option() in
   the namespace, which selects it, or returns why it cannot be selected,
   and then raises that error, naming the option, through .sync_error().
   One option lookup when nothing changed. R code run by the sync cannot
   sync again; the error comes after it, so its calling handlers can. An
   error leaves nothing set, so the next call tries again. */
static void impl_sync(void) {
  SEXP c = impl_option(), cont, problem;
  if (impl_syncing || (impl_synced && c == impl_seen)) return;
  impl_syncing = 1;
  cont = PROTECT(R_MakeUnwindCont());
  problem = PROTECT(R_UnwindProtect(impl_sync_eval, NULL, impl_sync_done, NULL, cont));
  if (problem != R_NilValue) {
    ns_eval(PROTECT(Rf_lang2(Rf_install(".sync_error"), problem)));
    UNPROTECT(1); /* not reached: .sync_error() always errors */
  }
  UNPROTECT(2);
  c = impl_option();
  if (c == R_UnboundValue) return;
  if (impl_seen != NULL && impl_seen != R_NilValue) R_ReleaseObject(impl_seen);
  if (c != R_NilValue) R_PreserveObject(c);
  impl_seen = c;
  impl_synced = 1;
}

void rsimd_impl_selected(void) {
  impl_synced = 0;
}

/* The index of the single string v in names (count entries), or -1. */
static int option_choice(SEXP v, const char *const *names, int count) {
  int i;
  if (TYPEOF(v) != STRSXP || XLENGTH(v) != 1 || STRING_ELT(v, 0) == NA_STRING) return -1;
  for (i = 0; i < count; i++) {
    if (strcmp(CHAR(STRING_ELT(v, 0)), names[i]) == 0) return i;
  }
  return -1;
}

int rsimd_arg_precision(SEXP precision) {
  static const char *const modes[] = {"fast", "pairwise", "compensated"};
  int p;
  if (precision == R_NilValue) {
    SEXP v = Rf_GetOption1(sym_precision);
    if (v == R_NilValue) return RSIMD_PREC_FAST;
    p = option_choice(v, modes, 3);
    if (p < 0) {
      Rf_errorcall(R_NilValue, "invalid option rsimd.precision: precision mode must be one "
                               "of \"fast\", \"pairwise\", \"compensated\"; reset it "
                               "with simd_precision()");
    }
    return p;
  }
  p = rsimd_arg_int1(precision, "precision");
  if (p < RSIMD_PREC_FAST || p > RSIMD_PREC_COMPENSATED) {
    Rf_error("internal error: invalid precision code %d", p);
  }
  return p;
}

int rsimd_arg_accuracy(SEXP accuracy) {
  static const char *const modes[] = {"accurate", "fast"};
  int a;
  if (accuracy == R_NilValue) {
    SEXP v = Rf_GetOption1(sym_math_accuracy);
    if (v == R_NilValue) return 0;
    a = option_choice(v, modes, 2);
    if (a < 0) {
      Rf_errorcall(R_NilValue, "invalid option rsimd.math_accuracy: math accuracy mode must "
                               "be one of \"accurate\", \"fast\"; reset it with "
                               "simd_math_accuracy()");
    }
    return a;
  }
  a = rsimd_arg_int1(accuracy, "accuracy");
  if (a != 0 && a != 1) Rf_error("internal error: invalid accuracy code %d", a);
  return a;
}

/* Option rsimd.na_check: TRUE when unset. */
static int na_check_option(void) {
  SEXP v = Rf_GetOption1(sym_na_check);
  if (v == R_NilValue) return 1;
  if (TYPEOF(v) != LGLSXP || XLENGTH(v) != 1 || LOGICAL_ELT(v, 0) == NA_LOGICAL) {
    Rf_errorcall(R_NilValue, "invalid option rsimd.na_check: must be TRUE or FALSE; reset it "
                             "with options(rsimd.na_check = TRUE)");
  }
  return LOGICAL_ELT(v, 0);
}

void rsimd_entry(void) {
  impl_sync();
  rsimd_active = rsimd_selected;
  warn_pending.n = 0;
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

/* Records operand v (named arg) in sv_call: a simd_vec's valid NA-free
   flag sets v->no_na_hint, and its rsimd_impl pin switches rsimd_active to
   the pinned tier's table, erroring if the tier is not available or two
   operands are pinned to different tiers. */
static void sv_note(rsimd_in *v, const char *arg) {
  SEXP x = v->sx, impl;
  if (Rf_isObject(x) && Rf_inherits(x, "simd_vec")) {
    sv_call.seen = 1;
    if (rsimd_sv_flag(x) == 1) v->no_na_hint = 1;
    impl = Rf_getAttrib(x, sym_sv_impl);
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

/* Attribute `name` of x without Rf_getAttrib(), which marks the value it
   returns as shared (ENSURE_NAMEDMAX) and so would void every token it
   read. R_mapAttrib() is new in R 4.6.0; ATTRIB() was the way before.
   R CMD check reports ATTRIB() as non-API only from R 4.6.0 on, which
   never compiles that branch, so no supported R version flags it. */
#if R_VERSION >= R_Version(4, 6, 0)
struct sv_attr_find {
  SEXP name, value;
};

static SEXP sv_attr_visit(SEXP tag, SEXP value, void *data) {
  struct sv_attr_find *f = data;
  if (tag == f->name) f->value = value;
  return NULL;
}

static SEXP sv_attr(SEXP x, SEXP name) {
  struct sv_attr_find f = {name, R_NilValue};
  R_mapAttrib(x, sv_attr_visit, &f);
  return f.value;
}
#else
static SEXP sv_attr(SEXP x, SEXP name) {
  SEXP a;
  for (a = ATTRIB(x); a != R_NilValue; a = CDR(a)) {
    if (TAG(a) == name) return CAR(a);
  }
  return R_NilValue;
}
#endif

int rsimd_sv_flag(SEXP x) {
  SEXP token, flag;
  if (TYPEOF(x) == RAWSXP) return 1;
  token = sv_attr(x, sym_sv_na_token);
  if (TYPEOF(token) != EXTPTRSXP || R_ExternalPtrAddr(token) != (void *) x ||
      MAYBE_SHARED(token) || R_ExternalPtrTag(token) != sv_attr(x, R_ClassSymbol)) {
    return -1;
  }
  flag = R_ExternalPtrProtected(token);
  if (TYPEOF(flag) != LGLSXP || XLENGTH(flag) != 1 || LOGICAL_ELT(flag, 0) == NA_LOGICAL) {
    return -1;
  }
  return LOGICAL_ELT(flag, 0);
}

/* Removes the token of x. The value is replaced first: that drops the
   token's reference count, which removing the attribute alone would leave
   counted by the unlinked cell, so that the object a copy x shares the
   token with keeps its flag. */
static void sv_drop_token(SEXP x, SEXP token_sym) {
  if (sv_attr(x, token_sym) == R_NilValue) return;
  Rf_setAttrib(x, token_sym, Rf_ScalarLogical(FALSE));
  Rf_setAttrib(x, token_sym, R_NilValue);
}

SEXP rsimd_sv_stamp(SEXP x, int flag) {
  SEXP token_sym = sym_sv_na_token;
  PROTECT(x);
  if (flag < 0) {
    sv_drop_token(x, token_sym);
  } else {
    SEXP value, cls = sv_attr(x, R_ClassSymbol);
    /* The tag must be a class vector no other object has (as a shared
       result class or a constant of R code has): class(x) <- class(y),
       after an edit made while x had no class, would put it back. */
    if (MAYBE_SHARED(cls)) {
      cls = PROTECT(Rf_duplicate(cls));
      Rf_setAttrib(x, R_ClassSymbol, cls);
      UNPROTECT(1);
    }
    value = PROTECT(Rf_ScalarLogical(flag != 0));
    Rf_setAttrib(x, token_sym, PROTECT(R_MakeExternalPtr((void *) x, cls, value)));
    UNPROTECT(2);
  }
  UNPROTECT(1);
  return x;
}

SEXP rsimd_sv_release(SEXP x) {
  SEXP token_sym = sym_sv_na_token, token = sv_attr(x, token_sym);
  if (token == R_NilValue || (TYPEOF(token) == EXTPTRSXP && R_ExternalPtrAddr(token) == x)) {
    return x;
  }
  PROTECT(x);
  sv_drop_token(x, token_sym);
  UNPROTECT(1);
  return x;
}

SEXP rsimd_sv_bare(SEXP x) {
  SEXP out;
  const void *src;
  size_t size;
  R_xlen_t n = XLENGTH(x);
  switch (TYPEOF(x)) {
  case REALSXP: size = sizeof(double); break;
  case INTSXP:
  case LGLSXP: size = sizeof(int); break;
  case RAWSXP: size = 1; break;
  case CPLXSXP: size = sizeof(Rcomplex); break;
  default: return R_NilValue;
  }
  src = DATAPTR_OR_NULL(x);
  if (src == NULL && n > 0) return R_NilValue;
  out = Rf_allocVector(TYPEOF(x), n);
  if (n > 0) memcpy(rsimd_out_ptr(out), src, (size_t) n * size);
  return out;
}

SEXP rsimd_sv_result(SEXP out, int keeps_na_free) {
  SEXP impl_sym = sym_sv_impl;
  int na_free;
  if (!sv_call.seen) return out;
  if (MAYBE_REFERENCED(out)) out = Rf_shallow_duplicate(out);
  PROTECT(out);
  Rf_setAttrib(out, R_ClassSymbol,
               Rf_isObject(out) && Rf_inherits(out, "integer64") ? sv_cls_i64 : sv_cls);
  if (sv_call.pin != RSIMD_TIER_COUNT) {
    SEXP *impl = &sv_impl_value[sv_call.pin];
    if (*impl == NULL) {
      *impl = Rf_mkString(rsimd_tier_names[sv_call.pin]);
      R_PreserveObject(*impl);
    }
    Rf_setAttrib(out, impl_sym, *impl);
  } else if (sv_attr(out, impl_sym) != R_NilValue) {
    Rf_setAttrib(out, impl_sym, R_NilValue);
  }
  na_free = TYPEOF(out) == RAWSXP || (keeps_na_free && sv_call.na_free);
  rsimd_sv_stamp(out, na_free ? 1 : -1);
  UNPROTECT(1);
  return out;
}

/* ---- Chunked read access ------------------------------------------------ */

int rsimd_in_init(rsimd_in *v, SEXP x, const char *arg) {
  v->sx = x;
  v->type = rsimd_check_atomic(x, arg);
  v->n = XLENGTH(x);
  /* The hint is read before any data pointer: from R 4.6.1-patched on, a
     wrapper clears its sortedness and no-NA metadata when a pointer to its
     data is taken. */
  switch (v->type) {
  case RSIMD_F64: v->no_na_hint = REAL_NO_NA(x); break;
  case RSIMD_I32: v->no_na_hint = INTEGER_NO_NA(x); break;
  case RSIMD_LGL: v->no_na_hint = LOGICAL_NO_NA(x); break;
  case RSIMD_U8: v->no_na_hint = 1; break; /* raw has no NA */
  default: v->no_na_hint = 0; break;       /* REAL_NO_NA says nothing about int64 NA */
  }
  /* Never DATAPTR(): it would materialise a compact sequence. */
  v->ptr = DATAPTR_OR_NULL(x);
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
  sym_impl = Rf_install("rsimd.impl");
  sym_precision = Rf_install("rsimd.precision");
  sym_math_accuracy = Rf_install("rsimd.math_accuracy");
  sym_na_check = Rf_install("rsimd.na_check");
  sym_sv_impl = Rf_install("rsimd_impl");
  sym_sv_na_token = Rf_install("rsimd_na_token");
  sv_cls = Rf_mkString("simd_vec");
  R_PreserveObject(sv_cls);
  sv_cls_i64 = Rf_allocVector(STRSXP, 2);
  R_PreserveObject(sv_cls_i64);
  SET_STRING_ELT(sv_cls_i64, 0, Rf_mkChar("simd_vec"));
  SET_STRING_ELT(sv_cls_i64, 1, Rf_mkChar("integer64"));
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

void rsimd_ew_no_c128(const rsimd_ew *e) {
  int i;
  for (i = 0; i < e->k; i++) {
    if (e->in[i].type == RSIMD_C128) Rf_error("internal error: complex operand in a chunk loop");
  }
}

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

/* 1 while bit64's namespace is known to be loaded: set when
   rsimd_set_i64_class() loads it, cleared by a hook on bit64's unloading
   (.onLoad), so that a loaded bit64 costs nothing per call. */
static int bit64_loaded = 0;

SEXP C_simd_bit64_unloaded(void) {
  bit64_loaded = 0;
  return R_NilValue;
}

void rsimd_set_i64_class(SEXP out) {
  /* An integer64 object without bit64's methods registered sorts, prints
     and compares as the doubles its bits make, so bit64 is loaded first.
     Evaluating R code here is safe where allocating is: out is protected
     by the caller. */
  if (!bit64_loaded) {
    SEXP pkg = PROTECT(Rf_mkString("bit64"));
    SEXP call = PROTECT(Rf_lang3(Rf_install("requireNamespace"), pkg,
                                 Rf_ScalarLogical(1)));
    SET_TAG(CDDR(call), Rf_install("quietly"));
    bit64_loaded = Rf_asLogical(Rf_eval(call, R_BaseEnv)) == 1;
    UNPROTECT(2);
  }
  if (!bit64_loaded) {
    Rf_error("integer64 results need package 'bit64'; install it with "
             "install.packages(\"bit64\")");
  }
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
  rsimd_warn("NAs produced by integer64 overflow");
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
      rsimd_warn("no non-missing arguments to %s; returning %s", is_min ? "min" : "max",
                 is_min ? "Inf" : "-Inf");
      return Rf_ScalarReal(is_min ? R_PosInf : R_NegInf);
    }
    if (type == RSIMD_I64 && r->count == 0) {
      int is_min = op == RSIMD_RED_MIN;
      rsimd_warn("no non-NA value, returning the %s possible integer64 value %s9223372036854775807",
                 is_min ? "highest" : "lowest", is_min ? "+" : "-");
      return rsimd_scalar_i64(is_min ? INT64_MAX : -INT64_MAX);
    }
    if (type == RSIMD_F64) return Rf_ScalarReal(value);
    if (is_int) return Rf_ScalarInteger((int) r->i64);
    /* No INT64_MIN check as for the sum: INT64_MIN is NA on input, so
       an extremum is never it. */
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
  rsimd_warn("NAs produced by integer overflow");
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

/* ---- Options and arguments ---------------------------------------------- */

void rsimd_opts_init(rsimd_opts *o, SEXP na_rm, SEXP na_check, int no_na_hint) {
  rsimd_opts_init_fixed(
      o, rsimd_arg_na_rm(na_rm),
      Rf_isNull(na_check) ? na_check_option() : rsimd_arg_lgl1(na_check, "na_check"), no_na_hint);
}

void rsimd_opts_init_fixed(rsimd_opts *o, int na_rm, int na_check, int no_na_hint) {
  o->na_rm = na_rm;
  o->na_check = na_check && !no_na_hint;
  o->precision = RSIMD_PREC_FAST;
  o->extrema = RSIMD_EXT_BOTH;
}

int rsimd_arg_na_rm(SEXP na_rm) {
  return Rf_isNull(na_rm) ? 0 : rsimd_arg_lgl1(na_rm, "na.rm");
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

/* The index of name in table, or -1. */
static int find_name(const char *name, const char *const *table, int count) {
  int i;
  for (i = 0; i < count; i++) {
    if (table[i] != NULL && strcmp(name, table[i]) == 0) return i;
  }
  return -1;
}

int rsimd_lookup(const char *name, const char *const *table, int count, const char *what) {
  int i = find_name(name, table, count);
  if (i < 0) Rf_error("unknown %s '%s'", what, name);
  return i;
}

int rsimd_arg_choice(SEXP op, const char *const *names, int count) {
  const char *name = rsimd_arg_str(op, "op");
  int i = find_name(name, names, count);
  if (i < 0) Rf_error("internal error: unknown op '%s'", name);
  return i;
}
