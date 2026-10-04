#ifndef RSIMD_API_COMPLEX_H
#define RSIMD_API_COMPLEX_H

/* The complex paths of the add/sub, neg, sum and predicate entry points
   (api_complex.c), called once those have classified their operands. */

#include "rvec.h"

/* x + y or x - y (op RSIMD_EW_ADD or RSIMD_EW_SUB) of two complex
   vectors, with the length-1 broadcast rule. */
SEXP rsimd_c128_add(int op, SEXP x, SEXP y, const rsimd_opts *o);
/* -x of a complex input. */
SEXP rsimd_c128_neg(const rsimd_in *in);
/* sum(x) of a complex input: each part is NA if it saw an NA, else NaN if
   it saw a NaN (na.rm = FALSE), else its sum in o->precision. */
SEXP rsimd_c128_sum(const rsimd_in *in, const rsimd_opts *o);
/* The predicate `code` over a complex input: writes po for mode
   RSIMD_PRED_ELT, else returns the any/all answer. */
int rsimd_c128_pred(const rsimd_in *in, int code, int mode, int *po);

#endif /* RSIMD_API_COMPLEX_H */
