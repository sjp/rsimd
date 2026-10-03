/* .Call entry points of the access layer: the operand promotion used by the
   R side, and internal routines that let the tests drive the layer
   directly (they are not part of the package's interface). */

#include <string.h>
#include "rsimd.h"
#include "rvec.h"

/* c(<common type>, <type of x>, <type of y>) as rsimd_etype_names. */
SEXP C_simd_promote(SEXP x, SEXP y) {
  rsimd_etype tx = rsimd_check_atomic(x, "x"), ty = rsimd_check_atomic(y, "y");
  rsimd_etype to = rsimd_promote(tx, ty);
  SEXP out;
  if (to == RSIMD_BAD) {
    Rf_error("cannot combine 'x' of type %s with 'y' of type %s", rsimd_etype_names[tx],
             rsimd_etype_names[ty]);
  }
  out = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(out, 0, Rf_mkChar(rsimd_etype_names[to]));
  SET_STRING_ELT(out, 1, Rf_mkChar(rsimd_etype_names[tx]));
  SET_STRING_ELT(out, 2, Rf_mkChar(rsimd_etype_names[ty]));
  UNPROTECT(1);
  return out;
}

/* Reads x chunk by chunk and reports list(path, regions, sum, na,
   no_na_hint, stride, type): "contiguous" or "regions", the number of
   chunks, a sum of the elements in index order (NA skipped and counted for
   integer, logical and integer64; real plus imaginary parts for complex),
   so every path gives the identical value. Allocates nothing proportional
   to the input. */
SEXP C_simd_debug_regions(SEXP x) {
  const char *names[] = {"path", "regions", "sum", "na", "no_na_hint", "stride", "type", ""};
  rsimd_in in;
  double regions = 0, sum = 0, na = 0;
  SEXP out;

  rsimd_in_init(&in, x, "x");
  switch (in.type) {
  case RSIMD_F64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) sum += px[i];
    });
    break;
  case RSIMD_I64:
    RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) {
        int64_t v;
        memcpy(&v, px + i, sizeof v);
        if (v == INT64_MIN) na++;
        else sum += (double) v;
      }
    });
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) {
        if (px[i] == NA_INTEGER) na++;
        else sum += px[i];
      }
    });
    break;
  case RSIMD_U8:
    RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) sum += px[i];
    });
    break;
  case RSIMD_C128:
    RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off, {
      const double *d = rsimd_c128_as_f64(px, len);
      regions++;
      for (R_xlen_t i = 0; i < 2 * len; i++) sum += d[i];
    });
    break;
  default: break;
  }

  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, Rf_mkString(in.ptr != NULL ? "contiguous" : "regions"));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(regions));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(sum));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(na));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(in.no_na_hint));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal((double) rsimd_stride));
  SET_VECTOR_ELT(out, 6, Rf_mkString(rsimd_etype_names[in.type]));
  UNPROTECT(1);
  return out;
}

#define COPY_CHUNKS(T)                                                                    \
  RSIMD_FOREACH_CHUNK(&in, T, px, len, off, {                                             \
    memcpy((T *) dst + off, px, (size_t) len * sizeof(T));                                \
  })

/* A copy of x read through the chunk loop, allocated and given attributes
   as a result would be (no_na: the op guarantees an NA-free result). */
SEXP C_simd_debug_copy(SEXP x, SEXP no_na) {
  rsimd_in in;
  SEXP out;
  void *dst;
  rsimd_in_init(&in, x, "x");
  out = PROTECT(rsimd_alloc_like(in.type, in.n));
  dst = rsimd_out_ptr(out);
  switch (in.type) {
  case RSIMD_F64:
  case RSIMD_I64: COPY_CHUNKS(double); break;
  case RSIMD_I32:
  case RSIMD_LGL: COPY_CHUNKS(int); break;
  case RSIMD_U8: COPY_CHUNKS(Rbyte); break;
  case RSIMD_C128: COPY_CHUNKS(Rcomplex); break;
  default: break;
  }
  rsimd_copy_class(x, out, rsimd_arg_lgl1(no_na, "no_na"));
  UNPROTECT(1);
  return out;
}

/* x + y as double through the binary chunk loop, for double or
   integer/logical operands of one type (integer NA gives NA). Returns
   list(value, n, x_scalar, y_scalar, regions). */
SEXP C_simd_debug_bin(SEXP x, SEXP y) {
  const char *names[] = {"value", "n", "x_scalar", "y_scalar", "regions", ""};
  rsimd_bin b;
  double regions = 0, *dst;
  SEXP out, value;

  rsimd_bin_init(&b, x, y);
  value = PROTECT(rsimd_alloc_like(RSIMD_F64, b.n));
  dst = rsimd_out_ptr(value);
  if (b.x.type == RSIMD_F64 && b.y.type == RSIMD_F64) {
    RSIMD_FOREACH_CHUNK2(&b, double, px, py, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) {
        dst[off + i] = px[b.x_scalar ? 0 : i] + py[b.y_scalar ? 0 : i];
      }
    });
  } else if ((b.x.type == RSIMD_I32 || b.x.type == RSIMD_LGL) &&
             (b.y.type == RSIMD_I32 || b.y.type == RSIMD_LGL)) {
    RSIMD_FOREACH_CHUNK2(&b, int, px, py, len, off, {
      regions++;
      for (R_xlen_t i = 0; i < len; i++) {
        int a = px[b.x_scalar ? 0 : i], c = py[b.y_scalar ? 0 : i];
        dst[off + i] = a == NA_INTEGER || c == NA_INTEGER ? NA_REAL : (double) a + c;
      }
    });
  } else {
    Rf_error("debug_bin takes two double or two integer/logical operands");
  }

  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, value);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal((double) b.n));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(b.x_scalar));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(b.y_scalar));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(regions));
  UNPROTECT(2);
  return out;
}

static int lookup(const char *name, const char *const *table, int count, const char *what) {
  int i;
  for (i = 0; i < count; i++) {
    if (strcmp(name, table[i]) == 0) return i;
  }
  Rf_error("unknown %s '%s'", what, name);
  return -1; /* not reached */
}

static SEXP result_fields(const rsimd_reduce_result *r) {
  const char *names[] = {"f64",    "comp",    "i64",      "idx",      "count",
                         "saw_na", "saw_nan", "overflow", "any_true", "any_false", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(r->f64));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(r->comp));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal((double) r->i64));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal((double) r->idx));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal((double) r->count));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(r->saw_na));
  SET_VECTOR_ELT(out, 6, Rf_ScalarLogical(r->saw_nan));
  SET_VECTOR_ELT(out, 7, Rf_ScalarLogical(r->overflow));
  SET_VECTOR_ELT(out, 8, Rf_ScalarLogical(r->any_true));
  SET_VECTOR_ELT(out, 9, Rf_ScalarLogical(r->any_false));
  UNPROTECT(1);
  return out;
}

/* Starts reduction `op` from its identity, overrides the fields named in
   the list `fields` (f64, comp, i64, idx, count, saw_na, saw_nan, overflow,
   any_true, any_false) and finishes it as a reduction over n elements of
   `type` with precision code `precision`. Returns list(init, value): the
   identity's fields and the finished R value. */
SEXP C_simd_debug_finish(SEXP op, SEXP type, SEXP n, SEXP fields, SEXP precision) {
  const char *names[] = {"init", "value", ""};
  int o_idx = lookup(rsimd_arg_str(op, "op"), rsimd_reduce_op_names, RSIMD_RED_OP_COUNT,
                     "reduction");
  rsimd_etype t = (rsimd_etype) lookup(rsimd_arg_str(type, "type"), rsimd_etype_names,
                                       RSIMD_BAD, "type");
  R_xlen_t len = (R_xlen_t) rsimd_arg_dbl1(n, "n");
  SEXP fnames = Rf_getAttrib(fields, R_NamesSymbol), out;
  rsimd_reduce_result r;
  rsimd_opts o;
  R_xlen_t i;

  rsimd_opts_init(&o, R_NilValue, R_NilValue, precision, 0);
  rsimd_reduce_result_init(&r, o_idx);
  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, result_fields(&r));

  if (TYPEOF(fields) != VECSXP || (Rf_xlength(fields) > 0 && Rf_isNull(fnames))) {
    Rf_error("'fields' must be a named list");
  }
  for (i = 0; i < Rf_xlength(fields); i++) {
    const char *f = CHAR(STRING_ELT(fnames, i));
    double v = rsimd_arg_dbl1(VECTOR_ELT(fields, i), f);
    if (strcmp(f, "f64") == 0) r.f64 = v;
    else if (strcmp(f, "comp") == 0) r.comp = v;
    else if (strcmp(f, "i64") == 0) r.i64 = (int64_t) v;
    else if (strcmp(f, "idx") == 0) r.idx = (R_xlen_t) v;
    else if (strcmp(f, "count") == 0) r.count = (R_xlen_t) v;
    else if (strcmp(f, "saw_na") == 0) r.saw_na = v != 0;
    else if (strcmp(f, "saw_nan") == 0) r.saw_nan = v != 0;
    else if (strcmp(f, "overflow") == 0) r.overflow = v != 0;
    else if (strcmp(f, "any_true") == 0) r.any_true = v != 0;
    else if (strcmp(f, "any_false") == 0) r.any_false = v != 0;
    else Rf_error("unknown field '%s'", f);
  }
  SET_VECTOR_ELT(out, 1, rsimd_reduce_finish(o_idx, t, len, &r, &o));
  UNPROTECT(1);
  return out;
}

/* The rsimd_opts an entry point taking x and these arguments would use
   (NULL arguments take their defaults): list(na_rm, na_check, precision). */
SEXP C_simd_debug_opts(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  const char *names[] = {"na_rm", "na_check", "precision", ""};
  rsimd_in in;
  rsimd_opts o;
  SEXP out;
  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, precision, in.no_na_hint);
  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(o.na_rm));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(o.na_check));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(o.precision));
  UNPROTECT(1);
  return out;
}
