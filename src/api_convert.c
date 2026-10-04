/* .Call entry point of the type conversions simd_as_integer(),
   simd_as_double(), simd_as_logical() and simd_as_raw(): picks the
   conversion kernel from the input type and the target, runs it chunk by
   chunk and raises base R's coercion warnings once per call. Results are
   bare vectors. */

#include <string.h>
#include "rsimd.h"
#include "dispatch.h"
#include "rvec.h"

/* x converted to `to` ("integer", "double", "logical" or "raw") in `mode`
   (RSIMD_CVT_CHECKED, _SATURATING or _TRUNCATING; it applies to doubles
   converted to integer and to everything converted to raw). Input of the
   target type is copied. */
SEXP C_simd_convert(SEXP x, SEXP to, SEXP mode) {
  static const char *const targets[] = {"integer", "double", "logical", "raw"};
  static const rsimd_etype target_type[] = {RSIMD_I32, RSIMD_F64, RSIMD_LGL, RSIMD_U8};
  const char *name = rsimd_arg_str(to, "to");
  int m = rsimd_arg_int1(mode, "mode"), st = 0, target = -1, op = -1, i;
  rsimd_etype from;
  SEXP out;
  char *po;
  size_t out_size;
  rsimd_in in;

  for (i = 0; i < 4; i++) {
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
       : target == 3 ? RSIMD_CVT_F64_U8 : -1;
    break;
  case RSIMD_I32:
  case RSIMD_LGL:
    op = target == 1 ? RSIMD_CVT_I32_F64 : target == 2 ? RSIMD_CVT_I32_LGL
       : target == 3 ? RSIMD_CVT_I32_U8 : -1;
    if (in.type == RSIMD_LGL && target == 2) op = -1;
    break;
  case RSIMD_U8:
    op = target == 0 ? RSIMD_CVT_U8_I32 : target == 1 ? RSIMD_CVT_U8_F64
       : target == 2 ? RSIMD_CVT_U8_LGL : -1;
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
  /* In base R's order (CoercionWarning). */
  if (st & RSIMD_CVT_WARN_INT) Rf_warning("NAs introduced by coercion to integer range");
  if (st & RSIMD_CVT_WARN_RAW) Rf_warning("out-of-range values treated as 0 in coercion to raw");
  UNPROTECT(1);
  return out;
}
