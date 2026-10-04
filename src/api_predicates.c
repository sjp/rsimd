/* .Call entry points of the predicates (is.na, is.nan, is.finite,
   is.infinite, negative, zero; elementwise or as any/all) and of the
   elementwise comparisons. Results are bare logical vectors, or a single
   logical value for any/all. */

#include <string.h>
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

/* The predicate's value for every element of a type it is constant on
   (FALSE): a logical vector of n FALSE, or the any/all answer. */
static SEXP constant_false(R_xlen_t n, int mode) {
  SEXP out;
  if (mode == RSIMD_PRED_ANY) return Rf_ScalarLogical(0);
  if (mode == RSIMD_PRED_ALL) return Rf_ScalarLogical(n == 0);
  out = PROTECT(rsimd_alloc_like(RSIMD_LGL, n));
  if (n > 0) memset(rsimd_out_ptr(out), 0, (size_t) n * sizeof(int));
  UNPROTECT(1);
  return out;
}

/* Predicate `op` by name: "na" (NA or NaN), "nan" (NaN, not NA), "finite",
   "infinite", "negative" (sign bit set, not NaN; for integers < 0) and
   "zero", of a double, integer, logical or raw vector, for mode 0
   (elementwise), 1 (any) or 2 (all). Integers and logicals are never NaN
   or infinite; raw bytes are never NA, NaN, finite or infinite (as base
   R's is.finite says), and "negative" and "zero" do not take them.
   Complex elements take the first four, with base R's either-part rules
   (finite: both parts). integer64 elements are like integers (NA is
   INT64_MIN). Empty input gives FALSE for any and TRUE for all. */
SEXP C_simd_pred(SEXP x, SEXP op, SEXP mode) {
  static const char *const names[] = {"na", "nan", "finite", "infinite", "negative", "zero"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0]));
  int m = rsimd_arg_int1(mode, "mode"), res = m == RSIMD_PRED_ALL;
  SEXP out = R_NilValue;
  int *po = NULL;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  if (m < RSIMD_PRED_ELT || m > RSIMD_PRED_ALL) Rf_error("internal error: mode %d", m);
  switch (in.type) {
  case RSIMD_F64: break;
  case RSIMD_I32:
  case RSIMD_LGL:
  case RSIMD_I64:
    if (code == RSIMD_PRED_NAN || code == RSIMD_PRED_INFINITE) return constant_false(in.n, m);
    break;
  case RSIMD_U8:
    if (code == RSIMD_PRED_NEGATIVE || code == RSIMD_PRED_ZERO) {
      Rf_error("invalid 'type' (raw) of argument");
    }
    return constant_false(in.n, m);
  case RSIMD_C128:
    if (code == RSIMD_PRED_NEGATIVE || code == RSIMD_PRED_ZERO) {
      Rf_error("invalid 'type' (complex) of argument");
    }
    break;
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  if (m == RSIMD_PRED_ELT) {
    out = PROTECT(rsimd_alloc_like(RSIMD_LGL, in.n));
    po = (int *) rsimd_out_ptr(out);
  }
  /* any stops at the first chunk with a TRUE, all at the first with a
     FALSE; the kernels also stop inside the chunk. */
  if (in.type == RSIMD_C128) {
    res = rsimd_c128_pred(&in, code, m, po);
  } else if (in.type == RSIMD_I64) {
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      int r = rsimd_active->pred_i64(code, (const int64_t *) px, len, m,
                                     m == RSIMD_PRED_ELT ? po + off : NULL);
      if (m != RSIMD_PRED_ELT && r != res) {
        res = r;
        break;
      }
    });
  } else if (in.type == RSIMD_F64) {
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      int r = rsimd_active->pred_f64(code, px, len, m, m == RSIMD_PRED_ELT ? po + off : NULL);
      if (m != RSIMD_PRED_ELT && r != res) {
        res = r;
        break;
      }
    });
  } else {
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      int r = rsimd_active->pred_i32(code, px, len, m, m == RSIMD_PRED_ELT ? po + off : NULL);
      if (m != RSIMD_PRED_ELT && r != res) {
        res = r;
        break;
      }
    });
  }
  if (m != RSIMD_PRED_ELT) return Rf_ScalarLogical(res);
  UNPROTECT(1);
  return out;
}

/* Comparison `op` by name ("eq", "ne", "lt", "le", "gt", "ge") of two
   operands under the length-1 broadcast rule, as base R's ==, != ...:
   integer and logical operands compare as integers, anything with a
   double as doubles (integers converted exactly), raw with raw as unsigned
   bytes, integer64 with integer64, integer or logical as 64-bit integers.
   A missing operand gives NA. The R side converts a raw operand compared
   with a non-raw one to integer first, and an integer64 operand compared
   with a double to double. */
SEXP C_simd_cmp(SEXP x, SEXP y, SEXP op) {
  static const char *const names[] = {"eq", "ne", "lt", "le", "gt", "ge"};
  static const char *const args[] = {"x", "y"};
  int code = lookup_op(op, names, (int) (sizeof names / sizeof names[0]));
  SEXP sargs[2], out;
  rsimd_etype tx, ty;
  int *po;
  rsimd_ew e;

  sargs[0] = x;
  sargs[1] = y;
  rsimd_ew_init(&e, 2, sargs, args);
  tx = e.in[0].type;
  ty = e.in[1].type;
  out = PROTECT(rsimd_alloc_like(RSIMD_LGL, e.n));
  po = (int *) rsimd_out_ptr(out);
  if (tx == RSIMD_U8 && ty == RSIMD_U8) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_u8(code, (const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags,
                           po + off);
    });
  } else if ((tx == RSIMD_I64 || ty == RSIMD_I64) && (tx == RSIMD_I64 || is_int_like(tx)) &&
             (ty == RSIMD_I64 || is_int_like(ty))) {
    int flags = e.flags | (is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                (is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_i64(code, p[0], p[1], len, flags, po + off);
    });
  } else if (is_int_like(tx) && is_int_like(ty)) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_i32(code, (const int *) p[0], (const int *) p[1], len, e.flags,
                            po + off);
    });
  } else if ((tx == RSIMD_F64 || is_int_like(tx)) && (ty == RSIMD_F64 || is_int_like(ty))) {
    int flags = e.flags | (is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                (is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_f64(code, p[0], p[1], len, flags, po + off);
    });
  } else {
    Rf_error("comparison of these types is not implemented");
  }
  UNPROTECT(1);
  return out;
}
