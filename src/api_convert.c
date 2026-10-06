/* .Call entry point of the type conversions simd_as_integer(),
   simd_as_double(), simd_as_logical(), simd_as_raw() and
   simd_as_integer64(): picks the conversion kernel from the input type and
   the target, runs it chunk by chunk and raises base R's (or, for
   integer64, bit64's) coercion warnings once per call. Results are bare
   vectors, with class "integer64" for that target. */

#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

/* x converted to `to` ("integer", "double", "logical", "raw" or
   "integer64") in `mode` (RSIMD_CVT_CHECKED, _SATURATING or _TRUNCATING;
   it applies to doubles and integer64 converted to integer, doubles
   converted to integer64 and everything converted to raw). Input of the
   target type is copied. `quiet` (TRUE when the R side converts an
   integer64 operand mixed with a double, having warned itself) drops the
   integer64 precision warning. */
/* 1 if converting x to `to` in `mode` cannot make a missing value from a
   non-missing one: conversions to double, logical and raw, saturating ones,
   and widening ones (logical, raw and integer to integer or integer64;
   integer64 to integer64). */
static int convert_keeps_na_free(SEXP x, SEXP to, SEXP mode) {
  const char *name = rsimd_arg_str(to, "to");
  rsimd_etype from = rsimd_etype_of(x);
  if (strcmp(name, "double") == 0 || strcmp(name, "logical") == 0 || strcmp(name, "raw") == 0) {
    return 1;
  }
  if (rsimd_arg_int1(mode, "mode") == RSIMD_CVT_SATURATING) return 1;
  return from == RSIMD_LGL || from == RSIMD_U8 || from == RSIMD_I32 ||
         (from == RSIMD_I64 && strcmp(name, "integer64") == 0);
}

static SEXP simd_convert_impl(SEXP x, SEXP to, SEXP mode, SEXP quiet) {
  static const char *const targets[] = {"integer", "double", "logical", "raw", "integer64"};
  static const rsimd_etype target_type[] = {RSIMD_I32, RSIMD_F64, RSIMD_LGL, RSIMD_U8, RSIMD_I64};
  const char *name = rsimd_arg_str(to, "to");
  int m = rsimd_arg_int1(mode, "mode"), st = 0, target = -1, op = -1, i;
  int q = rsimd_arg_lgl1(quiet, "quiet");
  rsimd_etype from;
  SEXP out;
  char *po;
  size_t out_size;
  rsimd_in in;

  for (i = 0; i < 5; i++) {
    if (strcmp(name, targets[i]) == 0) target = i;
  }
  if (target < 0) Rf_error("internal error: unknown target '%s'", name);
  if (m < RSIMD_CVT_CHECKED || m > RSIMD_CVT_TRUNCATING) Rf_error("internal error: mode %d", m);
  rsimd_in_init(&in, x, "x");
  from = in.type;
  /* Logical storage is int32, so logical to integer is a copy. */
  if (from == RSIMD_LGL && target_type[target] == RSIMD_I32) from = RSIMD_I32;
  switch (from) {
  case RSIMD_F64:
    op = target == 0 ? RSIMD_CVT_F64_I32 : target == 2 ? RSIMD_CVT_F64_LGL
       : target == 3 ? RSIMD_CVT_F64_U8 : target == 4 ? RSIMD_CVT_F64_I64 : -1;
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    op = target == 1 ? RSIMD_CVT_I32_F64 : target == 2 ? RSIMD_CVT_I32_LGL
       : target == 3 ? RSIMD_CVT_I32_U8 : target == 4 ? RSIMD_CVT_I32_I64 : -1;
    if (in.type == RSIMD_LGL && target == 2) op = -1;
    break;
  case RSIMD_U8:
    op = target == 0 ? RSIMD_CVT_U8_I32 : target == 1 ? RSIMD_CVT_U8_F64
       : target == 2 ? RSIMD_CVT_U8_LGL : target == 4 ? RSIMD_CVT_U8_I64 : -1;
    break;
  case RSIMD_I64:
    op = target == 0 ? RSIMD_CVT_I64_I32 : target == 1 ? RSIMD_CVT_I64_F64
       : target == 2 ? RSIMD_CVT_I64_LGL : target == 3 ? RSIMD_CVT_I64_U8 : -1;
    break;
  default: Rf_error("invalid 'type' (%s) of argument", rsimd_etype_names[in.type]);
  }
  out = PROTECT(rsimd_alloc_like(target_type[target], in.n));
  po = (char *) rsimd_out_ptr(out);
  out_size = rsimd_etype_size(target_type[target]);
  if (op < 0) {
    /* Same element type: copy. */
    const size_t size = rsimd_etype_size(in.type);
    switch (in.type) {
    case RSIMD_F64:
    case RSIMD_I64:
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off,
                          { memcpy(po + (size_t) off * size, px, (size_t) len * size); });
      break;
    case RSIMD_U8:
      RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off,
                          { memcpy(po + (size_t) off * size, px, (size_t) len * size); });
      break;
    default:
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off,
                          { memcpy(po + (size_t) off * size, px, (size_t) len * size); });
      break;
    }
  } else {
    switch (in.type) {
    case RSIMD_F64:
    case RSIMD_I64:
      RSIMD_FOREACH_CHUNK(&in, double, px, len, off, {
        st |= rsimd_active->convert(op, m, px, len, po + (size_t) off * out_size);
      });
      break;
    case RSIMD_U8:
      RSIMD_FOREACH_CHUNK(&in, Rbyte, px, len, off, {
        st |= rsimd_active->convert(op, m, px, len, po + (size_t) off * out_size);
      });
      break;
    default:
      RSIMD_FOREACH_CHUNK(&in, int, px, len, off, {
        st |= rsimd_active->convert(op, m, px, len, po + (size_t) off * out_size);
      });
      break;
    }
  }
  if (target_type[target] == RSIMD_I64) rsimd_set_i64_class(out);
  /* In base R's order (CoercionWarning), then bit64's. */
  if (st & RSIMD_CVT_WARN_INT) rsimd_warn("NAs introduced by coercion to integer range");
  if (st & RSIMD_CVT_WARN_RAW) rsimd_warn("out-of-range values treated as 0 in coercion to raw");
  if (st & RSIMD_CVT_WARN_I64) rsimd_warn_i64_overflow();
  if (st & RSIMD_CVT_WARN_I32_OVF) rsimd_warn_int_overflow();
  if ((st & RSIMD_CVT_WARN_PRECISION) && !q) {
    rsimd_warn("integer precision lost while converting to double");
  }
  UNPROTECT(1);
  return out;
}

SEXP C_simd_convert(SEXP x, SEXP to, SEXP mode, SEXP quiet) {
  rsimd_entry();
  return rsimd_exit(
    rsimd_sv_result(simd_convert_impl(x, to, mode, quiet), convert_keeps_na_free(x, to, mode)));
}
