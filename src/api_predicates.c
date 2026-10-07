/* .Call entry points of the predicates (is.na, is.nan, is.finite,
   is.infinite, negative, zero, and the number classes normal, subnormal,
   whole, even, odd and pow2; elementwise or as any/all), of the
   elementwise comparisons and of the Hamming distances. Results are bare
   logical vectors, a single logical value for any/all, or a double count. */

#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "api_complex.h"

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
   "infinite", "negative" (sign bit set, not NaN; for integers < 0),
   "zero", "normal", "subnormal", "whole", "even", "odd" and "pow2" (see
   RSIMD_PRED_* in kernel_types.h), of a double, integer, logical or raw
   vector, for mode 0 (elementwise), 1 (any) or 2 (all). Integers and
   logicals are never NaN, infinite or subnormal; raw bytes are never NA,
   NaN, finite or infinite (as base R's is.finite says), and "negative" and
   the ops after it do not take them. Complex elements take the first
   four, with base R's either-part rules (finite: both parts). integer64
   elements are like integers (NA is INT64_MIN). Empty input gives FALSE
   for any and TRUE for all. */
static SEXP simd_pred_impl(SEXP x, SEXP op, SEXP mode) {
  static const char *const names[] = {"na",   "nan",  "finite", "infinite",  "negative", "zero",
                                      "normal", "subnormal", "whole", "even", "odd", "pow2"};
  int code = rsimd_arg_choice(op, names, (int) (sizeof names / sizeof names[0]));
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
    if (code == RSIMD_PRED_NAN || code == RSIMD_PRED_INFINITE || code == RSIMD_PRED_SUBNORMAL) {
      return constant_false(in.n, m);
    }
    break;
  case RSIMD_U8:
    if (code >= RSIMD_PRED_NEGATIVE) Rf_error("invalid 'type' (raw) of argument");
    return constant_false(in.n, m);
  case RSIMD_C128:
    if (code >= RSIMD_PRED_NEGATIVE) {
      Rf_error("invalid 'type' (complex) of argument");
    }
    break;
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  /* An input known to have no missing value (as for na_check) needs no
     scan for one. */
  if (code == RSIMD_PRED_NA && m != RSIMD_PRED_ELT && in.no_na_hint && in.type != RSIMD_C128) {
    return Rf_ScalarLogical(m == RSIMD_PRED_ALL && in.n == 0);
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

SEXP C_simd_pred(SEXP x, SEXP op, SEXP mode) {
  rsimd_entry();
  return rsimd_exit(simd_pred_impl(x, op, mode));
}

/* Comparison `op` by name ("eq", "ne", "lt", "le", "gt", "ge") of two
   operands under the length-1 broadcast rule, as base R's ==, != ...:
   integer and logical operands compare as integers, anything with a
   double as doubles (integers converted exactly), raw with raw as unsigned
   bytes, integer64 with integer64, integer or logical as 64-bit integers.
   A missing operand gives NA. Complex operands (both complex: the R side
   converts the other) compare both parts for eq and ne, NA when any part
   is NA or NaN; the other comparisons give base R's error. The R side
   converts a raw operand compared with a non-raw one to integer first,
   and an integer64 operand compared with a double to double. */
static SEXP simd_cmp_impl(SEXP x, SEXP y, SEXP op) {
  static const char *const names[] = {"eq", "ne", "lt", "le", "gt", "ge"};
  static const char *const args[] = {"x", "y"};
  int code = rsimd_arg_choice(op, names, (int) (sizeof names / sizeof names[0]));
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
  if (tx == RSIMD_C128 || ty == RSIMD_C128) {
    if (code != RSIMD_CMP_EQ && code != RSIMD_CMP_NE) {
      Rf_error("invalid comparison with complex values");
    }
    rsimd_c128_cmp(code, x, y, po);
  } else if (tx == RSIMD_U8 && ty == RSIMD_U8) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_u8(code, (const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags,
                           po + off);
    });
  } else if ((tx == RSIMD_I64 || ty == RSIMD_I64) && (tx == RSIMD_I64 || rsimd_is_int_like(tx)) &&
             (ty == RSIMD_I64 || rsimd_is_int_like(ty))) {
    int flags = e.flags | (rsimd_is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                (rsimd_is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_i64(code, p[0], p[1], len, flags, po + off);
    });
  } else if (rsimd_is_int_like(tx) && rsimd_is_int_like(ty)) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_i32(code, (const int *) p[0], (const int *) p[1], len, e.flags,
                            po + off);
    });
  } else if ((tx == RSIMD_F64 || rsimd_is_int_like(tx)) && (ty == RSIMD_F64 || rsimd_is_int_like(ty))) {
    int flags = e.flags | (rsimd_is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                (rsimd_is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->cmp_f64(code, p[0], p[1], len, flags, po + off);
    });
  } else {
    Rf_error("comparison of these types is not implemented");
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cmp(SEXP x, SEXP y, SEXP op) {
  rsimd_entry();
  return rsimd_exit(simd_cmp_impl(x, y, op));
}

/* The number of pairs with x != y (simd_hamming), as a double: NA when a
   pair has a missing operand, unless na.rm = TRUE skips those pairs.
   Operands follow the comparison rules of simd_cmp_impl (the R side has
   already converted raw with non-raw, integer64 with double, and anything
   with complex); complex pairs differ when either part does. */
static SEXP simd_hamming_impl(SEXP x, SEXP y, SEXP na_rm) {
  static const char *const args[] = {"x", "y"};
  rsimd_reduce_result r;
  rsimd_opts o;
  SEXP sargs[2];
  rsimd_etype tx, ty;

  rsimd_opts_init(&o, na_rm, Rf_ScalarLogical(TRUE), 0);
  memset(&r, 0, sizeof r);
  if (TYPEOF(x) == CPLXSXP || TYPEOF(y) == CPLXSXP) {
    rsimd_bin b;
    int flags;
    rsimd_bin_init(&b, x, y);
    if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
      Rf_error("internal error: complex Hamming distance of %s and %s",
               rsimd_etype_names[b.x.type], rsimd_etype_names[b.y.type]);
    }
    flags = (b.x_scalar ? RSIMD_EW_SCALAR(0) : 0) | (b.y_scalar ? RSIMD_EW_SCALAR(1) : 0);
    RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off, {
      rsimd_active->hamming_c128(px, py, len, flags, &r, &o);
      if (r.saw_na && !o.na_rm) break;
    });
  } else {
    rsimd_ew e;
    sargs[0] = x;
    sargs[1] = y;
    rsimd_ew_init(&e, 2, sargs, args);
    tx = e.in[0].type;
    ty = e.in[1].type;
    if (tx == RSIMD_U8 && ty == RSIMD_U8) {
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->hamming_u8((const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags, &r);
      });
    } else if ((tx == RSIMD_I64 || ty == RSIMD_I64) && (tx == RSIMD_I64 || rsimd_is_int_like(tx)) &&
               (ty == RSIMD_I64 || rsimd_is_int_like(ty))) {
      int flags = e.flags | (rsimd_is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                  (rsimd_is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->hamming_i64(p[0], p[1], len, flags, &r, &o);
        if (r.saw_na && !o.na_rm) break;
      });
    } else if (rsimd_is_int_like(tx) && rsimd_is_int_like(ty)) {
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->hamming_i32((const int *) p[0], (const int *) p[1], len, e.flags, &r, &o);
        if (r.saw_na && !o.na_rm) break;
      });
    } else if ((tx == RSIMD_F64 || rsimd_is_int_like(tx)) && (ty == RSIMD_F64 || rsimd_is_int_like(ty))) {
      int flags = e.flags | (rsimd_is_int_like(tx) ? RSIMD_EW_I32(0) : 0) |
                  (rsimd_is_int_like(ty) ? RSIMD_EW_I32(1) : 0);
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->hamming_f64(p[0], p[1], len, flags, &r, &o);
        if (r.saw_na && !o.na_rm) break;
      });
    } else {
      Rf_error("Hamming distance of these types is not implemented");
    }
  }
  return Rf_ScalarReal(r.saw_na && !o.na_rm ? NA_REAL : (double) r.i64);
}

SEXP C_simd_hamming(SEXP x, SEXP y, SEXP na_rm) {
  rsimd_entry();
  return rsimd_exit(simd_hamming_impl(x, y, na_rm));
}

/* The number of differing bits of x and y (simd_hamming_bits), as a
   double: the population count of x ^ y over the pairs, NA elements by
   their bit pattern. Both operands are integer or logical, both raw or
   both integer64 (the R side checks). */
static SEXP simd_hamming_bits_impl(SEXP x, SEXP y) {
  static const char *const args[] = {"x", "y"};
  rsimd_reduce_result r;
  SEXP sargs[2];
  rsimd_etype tx, ty;
  rsimd_ew e;

  memset(&r, 0, sizeof r);
  sargs[0] = x;
  sargs[1] = y;
  rsimd_ew_init(&e, 2, sargs, args);
  tx = e.in[0].type;
  ty = e.in[1].type;
  if (tx == RSIMD_U8 && ty == RSIMD_U8) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->hamming_bits_u8((const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags, &r);
    });
  } else if (rsimd_is_int_like(tx) && rsimd_is_int_like(ty)) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->hamming_bits_i32((const int *) p[0], (const int *) p[1], len, e.flags, &r);
    });
  } else if (tx == RSIMD_I64 && ty == RSIMD_I64) {
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->hamming_bits_i64((const int64_t *) p[0], (const int64_t *) p[1], len, e.flags,
                                     &r);
    });
  } else {
    Rf_error("internal error: bit Hamming distance of %s and %s", rsimd_etype_names[tx],
             rsimd_etype_names[ty]);
  }
  return Rf_ScalarReal((double) r.i64);
}

SEXP C_simd_hamming_bits(SEXP x, SEXP y) {
  rsimd_entry();
  return rsimd_exit(simd_hamming_bits_impl(x, y));
}
