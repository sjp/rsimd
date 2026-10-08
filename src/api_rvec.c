/* .Call entry points of the access layer: the operand promotion used by the
   R side, and internal routines that let the tests drive the layer
   directly (they are not part of the package's interface). */

#include <string.h>
#include "rsimd.h"
#include "rvec.h"
#include "dispatch.h"

/* c(<common type>, <type of x>, <type of y>) as rsimd_etype_names. */
static SEXP simd_promote_impl(SEXP x, SEXP y) {
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

SEXP C_simd_promote(SEXP x, SEXP y) {
  rsimd_entry();
  return rsimd_exit(simd_promote_impl(x, y));
}

/* Reads x chunk by chunk and reports list(path, regions, sum, na,
   no_na_hint, stride, type): "contiguous" or "regions", the number of
   chunks, a sum of the elements in index order (NA skipped and counted for
   integer, logical and integer64; real plus imaginary parts for complex),
   so every path gives the identical value. Allocates nothing proportional
   to the input. */
static SEXP simd_debug_regions_impl(SEXP x) {
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

SEXP C_simd_debug_regions(SEXP x) {
  rsimd_entry();
  return rsimd_exit(simd_debug_regions_impl(x));
}

#define COPY_CHUNKS(T)                                                                    \
  RSIMD_FOREACH_CHUNK(&in, T, px, len, off, {                                             \
    memcpy((T *) dst + off, px, (size_t) len * sizeof(T));                                \
  })

/* A copy of x read through the chunk loop, allocated and given attributes
   as a result would be (no_na: the op guarantees an NA-free result): bare,
   except that an integer64 x (a plain or a simd_vec one) gives class
   "integer64", and a simd_vec x a simd_vec (rsimd_sv_result()). */
static SEXP simd_debug_copy_impl(SEXP x, SEXP no_na) {
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
  if (in.type == RSIMD_I64) rsimd_set_i64_class(out);
  out = rsimd_sv_result(out, rsimd_arg_lgl1(no_na, "no_na"));
  UNPROTECT(1);
  return out;
}

SEXP C_simd_debug_copy(SEXP x, SEXP no_na) {
  rsimd_entry();
  return rsimd_exit(simd_debug_copy_impl(x, no_na));
}

/* x + y as double through the binary chunk loop, for double or
   integer/logical operands of one type (integer NA gives NA). Returns
   list(value, n, x_scalar, y_scalar, regions). */
static SEXP simd_debug_bin_impl(SEXP x, SEXP y) {
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

SEXP C_simd_debug_bin(SEXP x, SEXP y) {
  rsimd_entry();
  return rsimd_exit(simd_debug_bin_impl(x, y));
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
   `type` with precision code `precision` and na.rm `na_rm`. Returns
   list(init, value): the identity's fields and the finished R value. */
static SEXP simd_debug_finish_impl(SEXP op, SEXP type, SEXP n, SEXP fields, SEXP precision,
                                   SEXP na_rm) {
  const char *names[] = {"init", "value", ""};
  int o_idx = rsimd_lookup(rsimd_arg_str(op, "op"), rsimd_reduce_op_names, RSIMD_RED_OP_COUNT,
                     "reduction");
  rsimd_etype t = (rsimd_etype) rsimd_lookup(rsimd_arg_str(type, "type"), rsimd_etype_names,
                                       RSIMD_BAD, "type");
  R_xlen_t len = (R_xlen_t) rsimd_arg_dbl1(n, "n");
  SEXP fnames = Rf_getAttrib(fields, R_NamesSymbol), out;
  rsimd_reduce_result r;
  rsimd_opts o;
  R_xlen_t i;

  rsimd_opts_init_fixed(&o, rsimd_arg_na_rm(na_rm), 1, 0);
  o.precision = rsimd_arg_precision(precision);
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

SEXP C_simd_debug_finish(SEXP op, SEXP type, SEXP n, SEXP fields, SEXP precision, SEXP na_rm) {
  rsimd_entry();
  return rsimd_exit(simd_debug_finish_impl(op, type, n, fields, precision, na_rm));
}

/* The rsimd_opts an entry point taking x and these arguments would use
   (NULL arguments take their defaults): list(na_rm, na_check, precision). */
static SEXP simd_debug_opts_impl(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  const char *names[] = {"na_rm", "na_check", "precision", ""};
  rsimd_in in;
  rsimd_opts o;
  SEXP out;
  rsimd_in_init(&in, x, "x");
  rsimd_opts_init(&o, na_rm, na_check, in.no_na_hint);
  o.precision = rsimd_arg_precision(precision);
  out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(o.na_rm));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(o.na_check));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(o.precision));
  UNPROTECT(1);
  return out;
}

SEXP C_simd_debug_opts(SEXP x, SEXP na_rm, SEXP na_check, SEXP precision) {
  rsimd_entry();
  return rsimd_exit(simd_debug_opts_impl(x, na_rm, na_check, precision));
}

/* The tier whose resolved table a call with operands x and y (y may be
   NULL) runs on: the selected tier, or the operands' pinned one. */
SEXP C_simd_debug_active(SEXP x, SEXP y) {
  rsimd_in in;
  int t;
  rsimd_entry();
  rsimd_in_init(&in, x, "x");
  if (!Rf_isNull(y)) rsimd_in_init(&in, y, "y");
  for (t = 0; t < RSIMD_TIER_COUNT; t++) {
    if (rsimd_tier_resolved((rsimd_tier) t) == rsimd_active) {
      return rsimd_exit(Rf_mkString(rsimd_tier_names[t]));
    }
  }
  Rf_error("internal error: the active table belongs to no tier");
  return R_NilValue; /* not reached */
}

/* The NA-free flag of simd_vec x (TRUE, FALSE or NULL), valid or not as
   rsimd_sv_flag() decides. */
SEXP C_simd_sv_flag(SEXP x) {
  int flag;
  rsimd_entry();
  flag = rsimd_sv_flag(x);
  return rsimd_exit(flag < 0 ? R_NilValue : Rf_ScalarLogical(flag));
}

/* x with its NA-free flag set to flag (TRUE, FALSE or NULL) and a fresh
   token; x is changed in place unless it is shared, as attr<- does. */
SEXP C_simd_sv_stamp(SEXP x, SEXP flag) {
  int f = Rf_isNull(flag) ? NA_LOGICAL : Rf_asLogical(flag);
  rsimd_entry();
  if (MAYBE_SHARED(x)) x = Rf_shallow_duplicate(x);
  return rsimd_exit(rsimd_sv_stamp(x, f == NA_LOGICAL ? -1 : f));
}

/* x without a token that belongs to another object (one R copied x from),
   changed in place. Replacement methods call it first: for f(y) <- v on a
   shared y, R hands them a copy of y, whose shared token would otherwise
   void the flag of y's other references for good. */
SEXP C_simd_sv_release(SEXP x) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_release(x));
}

/* TRUE when x is an ALTREP object, such as the wrapper R makes when it
   sets attributes on a shared long vector: unlist() reads its elements one
   at a time. */
SEXP C_simd_sv_altrep(SEXP x) {
  rsimd_entry();
  return rsimd_exit(Rf_ScalarLogical(ALTREP(x) != 0));
}

/* The data of x without attributes, or NULL; see rsimd_sv_bare(). */
SEXP C_simd_sv_bare(SEXP x) {
  rsimd_entry();
  return rsimd_exit(rsimd_sv_bare(x));
}
