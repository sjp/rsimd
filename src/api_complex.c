/* Complex numbers: the .Call entry point of Conj, Re and Im, and the
   complex paths of the add/sub, neg, sum and predicate entry points.

   Addition, subtraction and negation reuse the double kernels on the 2n
   doubles of the interleaved Rcomplex layout, so each part follows the
   double rules (NA payloads, na_check). A complex scalar operand is two
   doubles, which the double kernels' single-scalar broadcast cannot
   express: it is repeated into a small block that the kernel reads as an
   ordinary operand. Results are bare vectors. */

#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "api_complex.h"

/* Complex numbers per broadcast block: 4 KiB on the stack. */
#define RSIMD_C128_BCAST 256

static int lookup_op(SEXP op, const char *const *names, int count) {
  const char *name = rsimd_arg_str(op, "op");
  int i;
  for (i = 0; i < count; i++) {
    if (strcmp(name, names[i]) == 0) return i;
  }
  Rf_error("internal error: unknown op '%s'", name);
  return -1; /* not reached */
}

/* x op y (RSIMD_EW_ADD or RSIMD_EW_SUB) for one chunk of len elements;
   xs or ys is set when that operand is a scalar. */
static void c128_add_chunk(int op, const Rcomplex *px, const Rcomplex *py, R_xlen_t len, int xs,
                           int ys, Rcomplex *out, const rsimd_opts *o) {
  Rcomplex blk[RSIMD_C128_BCAST];
  R_xlen_t i;
  if (!xs && !ys) {
    rsimd_active->ew2_f64(op, rsimd_c128_as_f64(px, len), rsimd_c128_as_f64(py, len), 2 * len, 0,
                          (double *) out, o);
    return;
  }
  for (i = 0; i < RSIMD_C128_BCAST; i++) blk[i] = xs ? px[0] : py[0];
  for (i = 0; i < len; i += RSIMD_C128_BCAST) {
    R_xlen_t l = len - i < RSIMD_C128_BCAST ? len - i : RSIMD_C128_BCAST;
    rsimd_active->ew2_f64(op, (const double *) (xs ? blk : px + i),
                          (const double *) (ys ? blk : py + i), 2 * l, 0, (double *) (out + i), o);
  }
}

SEXP rsimd_c128_add(int op, SEXP x, SEXP y, const rsimd_opts *o) {
  rsimd_bin b;
  SEXP out;
  Rcomplex *po;

  rsimd_bin_init(&b, x, y);
  if (b.x.type != RSIMD_C128 || b.y.type != RSIMD_C128) {
    Rf_error("internal error: complex arithmetic on %s and %s", rsimd_etype_names[b.x.type],
             rsimd_etype_names[b.y.type]);
  }
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, b.n));
  po = (Rcomplex *) rsimd_out_ptr(out);
  RSIMD_FOREACH_CHUNK2(&b, Rcomplex, px, py, len, off, {
    c128_add_chunk(op, px, py, len, b.x_scalar, b.y_scalar, po + off, o);
  });
  UNPROTECT(1);
  return out;
}

SEXP rsimd_c128_neg(const rsimd_in *in) {
  SEXP out = PROTECT(rsimd_alloc_like(RSIMD_C128, in->n));
  Rcomplex *po = (Rcomplex *) rsimd_out_ptr(out);
  rsimd_opts o;

  rsimd_opts_init(&o, R_NilValue, R_NilValue, R_NilValue, 0);
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    rsimd_active->ew1_f64(RSIMD_EW_NEG, rsimd_c128_as_f64(px, len), 2 * len, 0,
                          (double *) (po + off), &o);
  });
  UNPROTECT(1);
  return out;
}

/* The value of one part of a finished sum: NA if that part saw an NA,
   else NaN if it saw a NaN, else the accumulated value (missing values
   are only tracked with na.rm = FALSE and na_check). */
static double part_value(const rsimd_reduce_result *r, int precision) {
  if (r->saw_na) return rsimd_na_real();
  if (r->saw_nan) return R_NaN;
  return rsimd_reduce_value(r, precision);
}

SEXP rsimd_c128_sum(const rsimd_in *in, const rsimd_opts *o) {
  rsimd_reduce_result r[2];
  SEXP out;
  Rcomplex *po;

  rsimd_reduce_result_init(&r[0], RSIMD_RED_SUM);
  rsimd_reduce_result_init(&r[1], RSIMD_RED_SUM);
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, { rsimd_active->sum_c128(px, len, r, o); });
  out = PROTECT(rsimd_alloc_like(RSIMD_C128, 1));
  po = (Rcomplex *) rsimd_out_ptr(out);
  po->r = part_value(&r[0], o->precision);
  po->i = part_value(&r[1], o->precision);
  UNPROTECT(1);
  return out;
}

int rsimd_c128_pred(const rsimd_in *in, int code, int mode, int *po) {
  int res = mode == RSIMD_PRED_ALL;
  RSIMD_FOREACH_CHUNK(in, Rcomplex, px, len, off, {
    int r = rsimd_active->pred_c128(code, px, len, mode, mode == RSIMD_PRED_ELT ? po + off : NULL);
    if (mode != RSIMD_PRED_ELT && r != res) {
      res = r;
      break;
    }
  });
  return res;
}

/* Conj(z), Re(z) or Im(z) by name. Complex z gives a complex conjugate
   and double parts; double, integer and logical z give double, as in base
   R: Conj and Re are the values (a bare copy), Im is zeros. Raw is an
   error with base R's message; integer64 is rejected on the R side. */
static SEXP simd_cplx_impl(SEXP z, SEXP op) {
  static const char *const names[] = {"conj", "re", "im"};
  int code = lookup_op(op, names, 3);
  SEXP out;
  rsimd_in in;

  rsimd_in_init(&in, z, "z");
  switch (in.type) {
  case RSIMD_C128:
    if (code == 0) {
      Rcomplex *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_C128, in.n));
      po = (Rcomplex *) rsimd_out_ptr(out);
      RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off,
                          { rsimd_active->conj_c128(px, len, po + off); });
    } else {
      double *po;
      out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
      po = (double *) rsimd_out_ptr(out);
      RSIMD_FOREACH_CHUNK(&in, Rcomplex, px, len, off,
                          { rsimd_active->part_c128(code == 2, px, len, po + off); });
    }
    break;
  case RSIMD_F64:
  case RSIMD_I32:
  case RSIMD_LGL: {
    double *po;
    out = PROTECT(rsimd_alloc_like(RSIMD_F64, in.n));
    po = (double *) rsimd_out_ptr(out);
    if (code == 2) {
      if (in.n > 0) memset(po, 0, (size_t) in.n * sizeof(double));
    } else if (in.type == RSIMD_F64) {
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                          { memcpy(po + off, px, (size_t) len * sizeof(double)); });
    } else {
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        rsimd_active->convert(RSIMD_CVT_I32_F64, RSIMD_CVT_CHECKED, px, len, po + off);
      });
    }
    break;
  }
  case RSIMD_U8: Rf_error("non-numeric argument to function");
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_cplx(SEXP z, SEXP op) {
  rsimd_entry();
  return rsimd_sv_result(simd_cplx_impl(z, op), 1);
}
