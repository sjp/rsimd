#ifndef RSIMD_API_COMPLEX_H
#define RSIMD_API_COMPLEX_H

/* The complex paths of the add/sub, neg, sum and predicate entry points
   (api_complex.c), called once those have classified their operands. */

#include "rvec.h"

/* x + y or x - y (op RSIMD_EW_ADD or RSIMD_EW_SUB) of two complex
   vectors, with the length-1 broadcast rule. */
SEXP rsimd_c128_add(int op, SEXP x, SEXP y, const rsimd_opts *o);
/* x * y or x / y (op RSIMD_EW_MUL or RSIMD_EW_DIV) of two complex
   vectors, bit-identical to base R's operators. */
SEXP rsimd_c128_muldiv(int op, SEXP x, SEXP y);
/* prod(x) and mean(x) of a complex input. prod: in o->precision; with a
   missing element (na.rm = FALSE) both parts are NA if one was NA, else
   NaN. mean: sum's parts (as rsimd_c128_sum) divided by the count. */
SEXP rsimd_c128_prod(const rsimd_in *in, const rsimd_opts *o);
SEXP rsimd_c128_mean(const rsimd_in *in, const rsimd_opts *o);
/* cumsum (op 0) or cumprod (op 1) of a complex input, as base R. */
SEXP rsimd_c128_scan(const rsimd_in *in, int op);
/* x == y or x != y (RSIMD_CMP_EQ or RSIMD_CMP_NE) of two complex vectors
   into po (logical, NA where any part is NA or NaN). */
void rsimd_c128_cmp(int op, SEXP x, SEXP y, int *po);
/* The names of the arithmetic variants chosen at load time: a named
   character vector (mul, div, cumprod). */
SEXP c128_variants(void);
/* The arithmetic variants chosen at load time, for the complex functions
   that multiply or divide as base R does. */
const rsimd_c128_arith *rsimd_c128_arith_get(void);
/* -x of a complex input. */
SEXP rsimd_c128_neg(const rsimd_in *in);
/* sum(x) of a complex input: each part is NA if it saw an NA, else NaN if
   it saw a NaN (na.rm = FALSE), else its sum in o->precision. */
SEXP rsimd_c128_sum(const rsimd_in *in, const rsimd_opts *o);
/* The predicate `code` over a complex input: writes po for mode
   RSIMD_PRED_ELT, else returns the any/all answer. */
int rsimd_c128_pred(const rsimd_in *in, int code, int mode, int *po);

#endif /* RSIMD_API_COMPLEX_H */
