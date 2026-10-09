/* .Call entry points of three-valued logic (&, |, xor, !), the bitwise
   ops, shifts, rotates and bit counts, and the population count total.
   Raw operands give raw results (bytewise), except for the counts, which
   are integer; raw cannot be mixed with other types (base R's error). */

#include <math.h>
#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

static const char *const raw_mix_msg =
  "operations are possible only for numeric, logical or complex types";

/* 1 if every operand is raw, 0 if none is; errors with base R's message
   for a mix, and for types other than `allowed` operands. */
static int all_raw(const rsimd_ew *e, int allow_double, int allow_i64) {
  int i, raw = 0;
  for (i = 0; i < e->k; i++) {
    rsimd_etype t = e->in[i].type;
    if (t == RSIMD_U8) raw++;
    else if (!rsimd_is_int_like(t) && !(allow_double && t == RSIMD_F64) &&
             !(allow_i64 && t == RSIMD_I64)) {
      Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[t]);
    }
  }
  if (raw != 0 && raw != e->k) Rf_error("%s", raw_mix_msg);
  return raw != 0;
}

static int init_operands(rsimd_ew *e, SEXP x, SEXP y) {
  static const char *const args[] = {"x", "y"};
  SEXP sargs[2];
  sargs[0] = x;
  sargs[1] = y;
  return rsimd_ew_init(e, Rf_isNull(y) ? 1 : 2, sargs, args);
}

/* Three-valued logic by name: "and", "or", "xor" or "not" (y NULL).
   Double, integer and logical operands are read as logical values and
   give a logical result; two raw operands give raw (bytewise). */
static SEXP simd_logic_impl(SEXP x, SEXP y, SEXP op) {
  static const char *const names[] = {"and", "or", "xor", "not"};
  int code = rsimd_arg_choice(op, names, (int) (sizeof names / sizeof names[0]));
  SEXP out;
  rsimd_ew e;

  if ((code == RSIMD_LOGIC_NOT) != Rf_isNull(y)) {
    Rf_error("internal error: operands of '%s'", names[code]);
  }
  init_operands(&e, x, y);
  if (all_raw(&e, 1, 0)) {
    /* The logic and bit op codes of and, or, xor and not agree. */
    Rbyte *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_U8, e.n));
    po = (Rbyte *) rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      rsimd_active->bit_u8(code, (const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags, 0,
                           po + off);
    });
  } else {
    int *po, i, dbl = 0, flags = e.flags;
    for (i = 0; i < e.k; i++) {
      if (e.in[i].type == RSIMD_F64) dbl = 1;
      else flags |= RSIMD_EW_I32(i);
    }
    out = PROTECT(rsimd_alloc_like(RSIMD_LGL, e.n));
    po = (int *) rsimd_out_ptr(out);
    if (dbl) {
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->logic_f64(code, p[0], p[1], len, flags, po + off);
      });
    } else {
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->logic_i32(code, (const int *) p[0], (const int *) p[1], len, e.flags,
                                po + off);
      });
    }
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_logic(SEXP x, SEXP y, SEXP op) {
  rsimd_entry();
  return rsimd_exit(simd_logic_impl(x, y, op));
}

/* Bitwise op by name: "and", "or", "xor" (x and y), "not", "shl", "shr",
   "sar", "rotl", "rotr" (x and the count k), "popcount", "lzcnt", "tzcnt"
   (x). Integer and logical operands give integer results; raw operands
   give raw, or integer counts; an integer64 operand (the other integer64,
   integer or logical) gives integer64, or integer counts. The R side has
   checked k: a shift by NA_integer_ (a count outside 0..31, or 0..63 for
   integer64) gives NA everywhere; raw counts are in range, and rotate
   counts are reduced modulo the width. */
static SEXP simd_bit_impl(SEXP x, SEXP y, SEXP op, SEXP k, SEXP na_check) {
  static const char *const names[] = {"and", "or",   "xor",      "not",   "shl",  "shr",
                                      "sar", "rotl", "rotr", "popcount", "lzcnt", "tzcnt"};
  int code = rsimd_arg_choice(op, names, (int) (sizeof names / sizeof names[0]));
  /* NA (a shift count out of range) is allowed. */
  double kd = Rf_isNull(k) ? 0.0 : rsimd_arg_num1(k, "n");
  int count = isnan(kd) ? NA_INTEGER : (int) kd;
  int counts = code >= RSIMD_BIT_POPCNT;
  SEXP out;
  rsimd_opts o;
  rsimd_ew e;

  if ((code <= RSIMD_BIT_XOR) == Rf_isNull(y)) {
    Rf_error("internal error: operands of '%s'", names[code]);
  }
  init_operands(&e, x, y);
  rsimd_opts_init(&o, R_NilValue, na_check, e.no_na_hint);
  if (e.in[0].type == RSIMD_I64 || (e.k == 2 && e.in[1].type == RSIMD_I64)) {
    int i, flags = e.flags;
    all_raw(&e, 0, 1);
    for (i = 0; i < e.k; i++) {
      if (rsimd_is_int_like(e.in[i].type)) flags |= RSIMD_EW_I32(i);
    }
    if (counts) {
      int *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_I32, e.n));
      po = (int *) rsimd_out_ptr(out);
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->bit_i64(code, p[0], p[1], len, flags, 0, po + off, &o);
      });
    } else {
      int64_t *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_I64, e.n));
      po = (int64_t *) rsimd_out_ptr(out);
      if (count == NA_INTEGER) {
        R_xlen_t j;
        for (j = 0; j < e.n; j++) po[j] = RSIMD_NA_I64;
      } else {
        RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
          rsimd_active->bit_i64(code, p[0], p[1], len, flags, count, po + off, &o);
        });
      }
      rsimd_set_i64_class(out);
    }
    UNPROTECT(1);
    return out;
  }
  if (all_raw(&e, 0, 0)) {
    void *po;
    if (code == RSIMD_BIT_SAR) Rf_error("invalid 'type' (raw) of argument");
    out = PROTECT(rsimd_alloc_like(counts ? RSIMD_I32 : RSIMD_U8, e.n));
    po = rsimd_out_ptr(out);
    RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
      void *dst = counts ? (void *) ((int *) po + off) : (void *) ((Rbyte *) po + off);
      rsimd_active->bit_u8(code, (const Rbyte *) p[0], (const Rbyte *) p[1], len, e.flags,
                           count, dst);
    });
  } else {
    int *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_I32, e.n));
    po = (int *) rsimd_out_ptr(out);
    if (count == NA_INTEGER) {
      R_xlen_t i;
      for (i = 0; i < e.n; i++) po[i] = NA_INTEGER;
    } else {
      RSIMD_FOREACH_CHUNK_EW(&e, p, len, off, {
        rsimd_active->bit_i32(code, (const int *) p[0], (const int *) p[1], len, e.flags,
                              count, po + off, &o);
      });
    }
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_bit(SEXP x, SEXP y, SEXP op, SEXP k, SEXP na_check) {
  static const char *const keeps[] = {"popcount", "lzcnt", "tzcnt", NULL};
  rsimd_entry();
  return rsimd_exit(rsimd_sv_result(simd_bit_impl(x, y, op, k, na_check), rsimd_str_in(op, keeps)));
}

/* The total number of set bits of an integer, logical, integer64 or raw
   vector, as a double: NA if an element is NA (unless na.rm = TRUE, which
   skips it). */
static SEXP simd_popcount_total_impl(SEXP x, SEXP na_rm, SEXP na_check) {
  rsimd_reduce_result r;
  rsimd_opts o;
  rsimd_in in;

  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  memset(&r, 0, sizeof r);
  if (in.type == RSIMD_U8) {
    RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off, { rsimd_active->popcnt_sum_u8(px, len, &r); });
  } else if (rsimd_is_int_like(in.type)) {
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      rsimd_active->popcnt_sum_i32(px, len, &r, &o);
      if (r.saw_na && !o.na_rm) break;
    });
  } else if (in.type == RSIMD_I64) {
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      rsimd_active->popcnt_sum_i64((const int64_t *) px, len, &r, &o);
      if (r.saw_na && !o.na_rm) break;
    });
  } else {
    Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  return Rf_ScalarReal(r.saw_na && !o.na_rm ? NA_REAL : (double) r.i64);
}

SEXP C_simd_popcount_total(SEXP x, SEXP na_rm, SEXP na_check) {
  rsimd_entry();
  return rsimd_exit(simd_popcount_total_impl(x, na_rm, na_check));
}
