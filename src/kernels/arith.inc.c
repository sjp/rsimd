/* Elementwise kernels: out[i] = op(x[i], y[i], z[i]) with the op codes,
 * operand flags and status bits of kernel_types.h and the NA, overflow and
 * division rules of na.h.
 *
 * The none tier is plain C over the scalar helpers of na.h and is the
 * reference the vector tiers are tested against. The vector tiers process
 * full vectors and then one predicated vector. In that last vector the
 * inactive lanes of an operand repeat its first active element, so they
 * compute what an active lane computes and cannot raise a status bit
 * (overflow, lo > hi) that no element raised.
 *
 * Double kernels read each operand as doubles or as int32 elements
 * (RSIMD_EW_I32(k)); an int32 NA becomes NA_real_ unless na_check is off,
 * when the caller promises there is none.
 */

#define RSIMD_EW_IS_SCALAR(f, k) (((f) & RSIMD_EW_SCALAR(k)) != 0)
#define RSIMD_EW_IS_I32(f, k) (((f) & RSIMD_EW_I32(k)) != 0)

/* Element i of operand k as a double. */
static inline double rsimd_ew_get(const void *p, int flags, int k, R_xlen_t i, int check) {
  R_xlen_t j = RSIMD_EW_IS_SCALAR(flags, k) ? 0 : i;
  if (RSIMD_EW_IS_I32(flags, k)) {
    int v = ((const int *) p)[j];
    return check && v == RSIMD_NA_I32 ? rsimd_na_real() : (double) v;
  }
  return ((const double *) p)[j];
}

#if RSIMD_TIER_IS(none)

/* Runs `expr`, which sets r from a (and b, c), for every element, with
   the operands read by get(p, k, i). */
#define RSIMD_EW_NONE_ELT(nargs, get, expr)                                      \
  for (i = 0; i < n; i++) {                                                     \
    double a = get(x, 0, i), b = 0.0, c = 0.0, r;                                \
    if ((nargs) > 1) b = get(y, 1, i);                                           \
    if ((nargs) > 2) c = get(z, 2, i);                                           \
    (void) b;                                                                    \
    (void) c;                                                                    \
    expr;                                                                        \
    out[i] = r;                                                                  \
  }
#define RSIMD_EW_GET_F64(p, k, i) (((const double *) (p))[i])
#define RSIMD_EW_GET_ANY(p, k, i) rsimd_ew_get(p, flags, k, i, check)
/* Double vectors (flags 0), the common case, have their own loop without
   the per-element operand tests, so that the scalar reference runs about
   as fast as base R's own loops. */
#define RSIMD_EW_NONE_LOOP(nargs, expr)                                          \
  if (flags == 0) {                                                             \
    RSIMD_EW_NONE_ELT(nargs, RSIMD_EW_GET_F64, expr)                             \
  } else {                                                                      \
    RSIMD_EW_NONE_ELT(nargs, RSIMD_EW_GET_ANY, expr)                             \
  }

int RSIMD_KERNEL(ew1_f64)(int op, const void *x, R_xlen_t n, int flags, double *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew1_f64)(int op, const void *x, R_xlen_t n, int flags, double *out,
                          const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int check = o->na_check;
  int st = 0;
  R_xlen_t i;
  switch (op) {
  case RSIMD_EW_NEG: RSIMD_EW_NONE_LOOP(1, r = -a) break;
  case RSIMD_EW_ABS: RSIMD_EW_NONE_LOOP(1, r = fabs(a)) break;
  case RSIMD_EW_SIGN: RSIMD_EW_NONE_LOOP(1, r = rsimd_sign_f64(a)) break;
  case RSIMD_EW_RECIP: RSIMD_EW_NONE_LOOP(1, r = 1.0 / a) break;
  case RSIMD_EW_SQRT:
    RSIMD_EW_NONE_LOOP(1, {
      if (a < 0) st |= RSIMD_EW_NAN_PRODUCED;
      r = sqrt(a);
    })
    break;
  case RSIMD_EW_FLOOR: RSIMD_EW_NONE_LOOP(1, r = floor(a)) break;
  case RSIMD_EW_CEIL: RSIMD_EW_NONE_LOOP(1, r = ceil(a)) break;
  case RSIMD_EW_TRUNC: RSIMD_EW_NONE_LOOP(1, r = trunc(a)) break;
  case RSIMD_EW_ROUND: RSIMD_EW_NONE_LOOP(1, r = nearbyint(a)) break;
  default: break;
  }
  return st;
}

void RSIMD_KERNEL(round_digits_f64)(const void *x, R_xlen_t n, int flags, double p10, double big,
                                    double *out, const rsimd_opts *o);
void RSIMD_KERNEL(round_digits_f64)(const void *x, R_xlen_t n, int flags, double p10, double big,
                                    double *out, const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int check = o->na_check;
  R_xlen_t i;
  RSIMD_EW_NONE_LOOP(1, r = rsimd_round_digits_f64(a, p10, big))
}
int RSIMD_KERNEL(ew2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          double *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          double *out, const rsimd_opts *o) {
  const void *z = NULL;
  const int check = o->na_check;
  R_xlen_t i;
/* Every merged op gives NaN for a NaN operand, so only a NaN result can
   need the NA merge. */
#define RSIMD_EW_MERGED(e)                                                       \
  r = (e);                                                                       \
  if (isnan(r) && check) r = rsimd_na_merge_f64(r, a, b)
  switch (op) {
  case RSIMD_EW_ADD: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(a + b)) break;
  case RSIMD_EW_SUB: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(a - b)) break;
  case RSIMD_EW_MUL: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(a * b)) break;
  case RSIMD_EW_DIV: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(a / b)) break;
  case RSIMD_EW_IDIV: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(rsimd_idiv_f64(a, b))) break;
  case RSIMD_EW_MOD: RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(rsimd_mod_f64(a, b))) break;
  case RSIMD_EW_COPYSIGN:
    RSIMD_EW_NONE_LOOP(2, RSIMD_EW_MERGED(rsimd_copysign_f64(a, b)))
    break;
  case RSIMD_EW_PMIN: RSIMD_EW_NONE_LOOP(2, r = rsimd_pmin_f64(a, b)) break;
  case RSIMD_EW_PMAX: RSIMD_EW_NONE_LOOP(2, r = rsimd_pmax_f64(a, b)) break;
  case RSIMD_EW_PMIN_NUM: RSIMD_EW_NONE_LOOP(2, r = rsimd_pmin_num_f64(a, b)) break;
  case RSIMD_EW_PMAX_NUM: RSIMD_EW_NONE_LOOP(2, r = rsimd_pmax_num_f64(a, b)) break;
  default: break;
  }
#undef RSIMD_EW_MERGED
  return 0;
}

int RSIMD_KERNEL(ew3_f64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, double *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_f64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, double *out, const rsimd_opts *o) {
  const int check = o->na_check;
  int st = 0;
  R_xlen_t i;
#define RSIMD_EW_MERGED(e)                                                       \
  r = (e);                                                                       \
  if (isnan(r) && check) r = rsimd_na_merge3_f64(r, a, b, c)
  switch (op) {
  case RSIMD_EW_FMA: RSIMD_EW_NONE_LOOP(3, RSIMD_EW_MERGED(rsimd_fma(a, b, c))) break;
  case RSIMD_EW_MUL_ADD:
  case RSIMD_EW_MUL_ADD_APPROX: RSIMD_EW_NONE_LOOP(3, RSIMD_EW_MERGED(a * b + c)) break;
  case RSIMD_EW_ADD_MUL: RSIMD_EW_NONE_LOOP(3, RSIMD_EW_MERGED((a + b) * c)) break;
  case RSIMD_EW_LERP: RSIMD_EW_NONE_LOOP(3, RSIMD_EW_MERGED(rsimd_fma(c, b, (1.0 - c) * a))) break;
  case RSIMD_EW_CLAMP:
    RSIMD_EW_NONE_LOOP(3, {
      if (b > c) st |= RSIMD_EW_LO_GT_HI;
      r = rsimd_pmin_f64(rsimd_pmax_f64(a, b), c);
    })
    break;
  default: break;
  }
#undef RSIMD_EW_MERGED
  return st;
}

int RSIMD_KERNEL(ew1_i32)(int op, const int *x, R_xlen_t n, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew1_i32)(int op, const int *x, R_xlen_t n, int *out, const rsimd_opts *o) {
  R_xlen_t i;
  (void) o;
  for (i = 0; i < n; i++) out[i] = op == RSIMD_EW_ABS ? rsimd_abs_i32(x[i]) : rsimd_neg_i32(x[i]);
  return 0;
}

int RSIMD_KERNEL(ew2_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew2_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int *out,
                          const rsimd_opts *o) {
  const int check = o->na_check;
  const R_xlen_t sx = RSIMD_EW_IS_SCALAR(flags, 0) ? 0 : 1, sy = RSIMD_EW_IS_SCALAR(flags, 1) ? 0 : 1;
  int ovf = 0;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int32_t a = x[i * sx], b = y[i * sy], r;
    switch (op) {
    case RSIMD_EW_ADD: r = rsimd_add_i32(a, b, check, &ovf); break;
    case RSIMD_EW_SUB: r = rsimd_sub_i32(a, b, check, &ovf); break;
    case RSIMD_EW_MUL: r = rsimd_mul_i32(a, b, check, &ovf); break;
    case RSIMD_EW_ADD_WRAP: r = rsimd_add_wrap_i32(a, b, check); break;
    case RSIMD_EW_SUB_WRAP: r = rsimd_sub_wrap_i32(a, b, check); break;
    case RSIMD_EW_MUL_WRAP: r = rsimd_mul_wrap_i32(a, b, check); break;
    case RSIMD_EW_IDIV: r = rsimd_intdiv_i32(a, b, 0, check); break;
    case RSIMD_EW_MOD: r = rsimd_intdiv_i32(a, b, 1, check); break;
    case RSIMD_EW_PMIN: r = rsimd_pmin_i32(a, b); break;
    case RSIMD_EW_PMAX: r = rsimd_pmax_i32(a, b); break;
    case RSIMD_EW_PMIN_NUM: r = rsimd_pmin_num_i32(a, b); break;
    case RSIMD_EW_PMAX_NUM: r = rsimd_pmax_num_i32(a, b); break;
    default: r = RSIMD_NA_I32; break;
    }
    out[i] = r;
  }
  return ovf ? RSIMD_EW_OVERFLOW : 0;
}

int RSIMD_KERNEL(ew3_i32)(int op, const int *x, const int *y, const int *z, R_xlen_t n,
                          int flags, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_i32)(int op, const int *x, const int *y, const int *z, R_xlen_t n,
                          int flags, int *out, const rsimd_opts *o) {
  const int check = o->na_check;
  const R_xlen_t sx = RSIMD_EW_IS_SCALAR(flags, 0) ? 0 : 1, sy = RSIMD_EW_IS_SCALAR(flags, 1) ? 0 : 1,
                 sz = RSIMD_EW_IS_SCALAR(flags, 2) ? 0 : 1;
  int ovf = 0, st = 0;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int32_t a = x[i * sx], b = y[i * sy], c = z[i * sz], r;
    switch (op) {
    /* An NA intermediate (an overflow, or an NA operand) stays NA without
       counting as a second overflow, as in base R's x * y + z. */
    case RSIMD_EW_MUL_ADD:
      r = rsimd_mul_i32(a, b, check, &ovf);
      if (r != RSIMD_NA_I32) r = rsimd_add_i32(r, c, check, &ovf);
      break;
    case RSIMD_EW_ADD_MUL:
      r = rsimd_add_i32(a, b, check, &ovf);
      if (r != RSIMD_NA_I32) r = rsimd_mul_i32(r, c, check, &ovf);
      break;
    case RSIMD_EW_CLAMP:
      if (b != RSIMD_NA_I32 && c != RSIMD_NA_I32 && b > c) st |= RSIMD_EW_LO_GT_HI;
      r = rsimd_pmin_i32(rsimd_pmax_i32(a, b), c);
      break;
    default: r = RSIMD_NA_I32; break;
    }
    out[i] = r;
  }
  return st | (ovf ? RSIMD_EW_OVERFLOW : 0);
}

#undef RSIMD_EW_NONE_LOOP
#undef RSIMD_EW_NONE_ELT
#undef RSIMD_EW_GET_F64
#undef RSIMD_EW_GET_ANY

#else /* vector tiers */

/* ---- int32 kernels ------------------------------------------------------- */

/* Operand k of an int32 kernel at i: full vector, or the predicated last
   vector whose inactive lanes repeat element i. */
#define RSIMD_EW_LDI(p, k, i)                                                    \
  (RSIMD_EW_IS_SCALAR(flags, k) ? bc##k : rsimd_vi32_loadu((const int32_t *) (p) + (i)))
#define RSIMD_EW_LDI_P(p, k, i, pg)                                              \
  (RSIMD_EW_IS_SCALAR(flags, k)                                                  \
     ? bc##k                                                                     \
     : rsimd_vi32_loadu_p(pg, (const int32_t *) (p) + (i), ((const int32_t *) (p))[i]))

/* Runs `expr`, which sets the rsimd_vi32 r from a, b (and c), over the
   chunk. */
#define RSIMD_EW_I32_LOOP(nargs, expr)                                           \
  do {                                                                           \
    ptrdiff_t i = 0;                                                             \
    for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {                       \
      rsimd_vi32 a = RSIMD_EW_LDI(x, 0, i), b = bc1, c = bc2, r;                 \
      if ((nargs) > 1) b = RSIMD_EW_LDI(y, 1, i);                                \
      if ((nargs) > 2) c = RSIMD_EW_LDI(z, 2, i);                                \
      (void) b;                                                                  \
      (void) c;                                                                  \
      expr;                                                                      \
      rsimd_vi32_storeu((int32_t *) out + i, r);                                 \
    }                                                                            \
    if (i < n) {                                                                 \
      rsimd_p32 pg = rsimd_p32_while(i, n);                                      \
      rsimd_vi32 a = RSIMD_EW_LDI_P(x, 0, i, pg), b = bc1, c = bc2, r;           \
      if ((nargs) > 1) b = RSIMD_EW_LDI_P(y, 1, i, pg);                          \
      if ((nargs) > 2) c = RSIMD_EW_LDI_P(z, 2, i, pg);                          \
      (void) b;                                                                  \
      (void) c;                                                                  \
      expr;                                                                      \
      rsimd_vi32_storeu_p(pg, (int32_t *) out + i, r);                           \
    }                                                                            \
  } while (0)

/* Checked add, sub or mul of a and b into r, NA in the lanes of the mask
   na0, where an operand is NA (with check) or where the result is out of
   range; overflow outside the NA lanes sets ovf. */
#define RSIMD_EW_I32_CHECKED(r, a, b, opfn, ovffn, na0)                          \
  do {                                                                           \
    rsimd_mi32 na_ = check ? rsimd_mi32_or(na0, rsimd_vi32_na2(a, b)) : (na0);   \
    rsimd_mi32 o_;                                                               \
    r = opfn(a, b);                                                              \
    o_ = ovffn(a, b, r);                                                         \
    if (rsimd_mi32_any(rsimd_mi32_andnot(na_, o_))) ovf = 1;                     \
    r = rsimd_vi32_set_na(r, rsimd_mi32_or(o_, na_));                            \
  } while (0)

/* Wrapping op, NA where an operand is NA (with check). */
#define RSIMD_EW_I32_WRAP(r, a, b, opfn)                                         \
  do {                                                                           \
    r = opfn(a, b);                                                              \
    if (check) r = rsimd_vi32_set_na(r, rsimd_vi32_na2(a, b));                   \
  } while (0)

int RSIMD_KERNEL(ew1_i32)(int op, const int *x, R_xlen_t n, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew1_i32)(int op, const int *x, R_xlen_t n, int *out, const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int flags = 0;
  const rsimd_vi32 bc0 = rsimd_vi32_zero(), bc1 = bc0, bc2 = bc0;
  (void) o;
  (void) y;
  (void) z;
  (void) bc0;
  if (op == RSIMD_EW_ABS) {
    RSIMD_EW_I32_LOOP(1, r = rsimd_vi32_abs_wrap(a));
  } else {
    RSIMD_EW_I32_LOOP(1, r = rsimd_vi32_neg_wrap(a));
  }
  return 0;
}

/* %/% or %% of an int32 vector x by a scalar d other than 0 and
   INT32_MIN, by multiplication (rsimd_vi32_intdiv_const, na.h). */
static inline void rsimd_ew_intdiv_const_i32(const int *x, int32_t d, R_xlen_t n, int mod,
                                             int check, int *out) {
  const rsimd_divmagic_i32 g = rsimd_divmagic_i32_make(d);
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {
    rsimd_vi32_storeu(out + i, rsimd_vi32_intdiv_const(rsimd_vi32_loadu(x + i), &g, mod, check));
  }
  if (i < n) {
    rsimd_p32 pg = rsimd_p32_while(i, n);
    rsimd_vi32_storeu_p(pg, out + i,
                        rsimd_vi32_intdiv_const(rsimd_vi32_loadu_p(pg, x + i, 0), &g, mod, check));
  }
}

/* %/% or %% of int32 operands through double lanes (exact, na.h), two
   vectors per step. Inlined with constant arguments for the common case,
   so that the tests of xs, ys, mod and check leave the loop. */
RSIMD_ALWAYS_INLINE void rsimd_ew_intdiv_i32_(const int *x, const int *y, R_xlen_t n, const int xs,
                                              const int ys, const int mod, const int check,
                                              int *out) {
#ifdef RSIMD_NO_F64_SIMD
  const R_xlen_t sx = xs ? 0 : 1, sy = ys ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) out[i] = rsimd_intdiv_i32(x[i * sx], y[i * sy], mod, check);
#else
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 bx = rsimd_vf64_set1((double) x[0]), by = rsimd_vf64_set1((double) y[0]);
  ptrdiff_t i = 0;
  for (; i + 2 * W <= n; i += 2 * W) {
    rsimd_vf64 a0 = xs ? bx : rsimd_vf64_loadu_i32(x + i);
    rsimd_vf64 a1 = xs ? bx : rsimd_vf64_loadu_i32(x + i + W);
    rsimd_vf64 b0 = ys ? by : rsimd_vf64_loadu_i32(y + i);
    rsimd_vf64 b1 = ys ? by : rsimd_vf64_loadu_i32(y + i + W);
    rsimd_vf64_storeu_i32(out + i, rsimd_vf64_intdiv(a0, b0, mod, check));
    rsimd_vf64_storeu_i32(out + i + W, rsimd_vf64_intdiv(a1, b1, mod, check));
  }
  for (; i < n; i += W) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64 a = xs ? bx : rsimd_vf64_loadu_i32_p(pg, x + i, x[i]);
    rsimd_vf64 b = ys ? by : rsimd_vf64_loadu_i32_p(pg, y + i, y[i]);
    rsimd_vf64_storeu_i32_p(pg, out + i, rsimd_vf64_intdiv(a, b, mod, check));
  }
#endif
}
static void rsimd_ew_intdiv_i32(const int *x, const int *y, R_xlen_t n, int flags, int mod,
                                int check, int *out) {
  const int xs = RSIMD_EW_IS_SCALAR(flags, 0), ys = RSIMD_EW_IS_SCALAR(flags, 1);
  if (xs || ys) rsimd_ew_intdiv_i32_(x, y, n, xs, ys, mod, check, out);
  else if (mod) {
    if (check) rsimd_ew_intdiv_i32_(x, y, n, 0, 0, 1, 1, out);
    else rsimd_ew_intdiv_i32_(x, y, n, 0, 0, 1, 0, out);
  } else {
    if (check) rsimd_ew_intdiv_i32_(x, y, n, 0, 0, 0, 1, out);
    else rsimd_ew_intdiv_i32_(x, y, n, 0, 0, 0, 0, out);
  }
}

int RSIMD_KERNEL(ew2_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew2_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int *out,
                          const rsimd_opts *o) {
  const void *z = NULL;
  const int check = o->na_check;
  const rsimd_vi32 bc0 = rsimd_vi32_set1(x[0]), bc1 = rsimd_vi32_set1(y[0]),
                   bc2 = rsimd_vi32_zero();
  const rsimd_mi32 nomask = rsimd_mi32_none();
  int ovf = 0;
  (void) z;
  switch (op) {
  case RSIMD_EW_ADD:
    RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_CHECKED(r, a, b, rsimd_vi32_add, rsimd_vi32_add_ovf, nomask));
    break;
  case RSIMD_EW_SUB:
    RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_CHECKED(r, a, b, rsimd_vi32_sub, rsimd_vi32_sub_ovf, nomask));
    break;
  case RSIMD_EW_MUL:
    RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_CHECKED(r, a, b, rsimd_vi32_mul, rsimd_vi32_mul_ovf, nomask));
    break;
  case RSIMD_EW_ADD_WRAP: RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_WRAP(r, a, b, rsimd_vi32_add)); break;
  case RSIMD_EW_SUB_WRAP: RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_WRAP(r, a, b, rsimd_vi32_sub)); break;
  case RSIMD_EW_MUL_WRAP: RSIMD_EW_I32_LOOP(2, RSIMD_EW_I32_WRAP(r, a, b, rsimd_vi32_mul)); break;
  case RSIMD_EW_IDIV:
  case RSIMD_EW_MOD:
    if (!RSIMD_EW_IS_SCALAR(flags, 0) && RSIMD_EW_IS_SCALAR(flags, 1) && y[0] != 0 &&
        y[0] != RSIMD_NA_I32) {
      rsimd_ew_intdiv_const_i32(x, y[0], n, op == RSIMD_EW_MOD, check, out);
    } else {
      rsimd_ew_intdiv_i32(x, y, n, flags, op == RSIMD_EW_MOD, check, out);
    }
    break;
  /* NA is INT32_MIN, so min propagates it and max ignores it by itself. */
  case RSIMD_EW_PMIN: RSIMD_EW_I32_LOOP(2, r = rsimd_vi32_min(a, b)); break;
  case RSIMD_EW_PMAX:
    RSIMD_EW_I32_LOOP(2, r = rsimd_vi32_set_na(rsimd_vi32_max(a, b), rsimd_vi32_na2(a, b)));
    break;
  case RSIMD_EW_PMIN_NUM:
    RSIMD_EW_I32_LOOP(2, {
      r = rsimd_vi32_blend(rsimd_vi32_min(a, b), b, rsimd_vi32_is_na(a));
      r = rsimd_vi32_blend(r, a, rsimd_vi32_is_na(b));
    });
    break;
  case RSIMD_EW_PMAX_NUM: RSIMD_EW_I32_LOOP(2, r = rsimd_vi32_max(a, b)); break;
  default: break;
  }
  return ovf ? RSIMD_EW_OVERFLOW : 0;
}

int RSIMD_KERNEL(ew3_i32)(int op, const int *x, const int *y, const int *z, R_xlen_t n,
                          int flags, int *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_i32)(int op, const int *x, const int *y, const int *z, R_xlen_t n,
                          int flags, int *out, const rsimd_opts *o) {
  const int check = o->na_check;
  const rsimd_vi32 bc0 = rsimd_vi32_set1(x[0]), bc1 = rsimd_vi32_set1(y[0]),
                   bc2 = rsimd_vi32_set1(z[0]);
  const rsimd_mi32 nomask = rsimd_mi32_none();
  int ovf = 0, lohi = 0;
  switch (op) {
  case RSIMD_EW_MUL_ADD:
    RSIMD_EW_I32_LOOP(3, {
      rsimd_vi32 t;
      RSIMD_EW_I32_CHECKED(t, a, b, rsimd_vi32_mul, rsimd_vi32_mul_ovf, nomask);
      RSIMD_EW_I32_CHECKED(r, t, c, rsimd_vi32_add, rsimd_vi32_add_ovf, rsimd_vi32_is_na(t));
    });
    break;
  case RSIMD_EW_ADD_MUL:
    RSIMD_EW_I32_LOOP(3, {
      rsimd_vi32 t;
      RSIMD_EW_I32_CHECKED(t, a, b, rsimd_vi32_add, rsimd_vi32_add_ovf, nomask);
      RSIMD_EW_I32_CHECKED(r, t, c, rsimd_vi32_mul, rsimd_vi32_mul_ovf, rsimd_vi32_is_na(t));
    });
    break;
  case RSIMD_EW_CLAMP:
    RSIMD_EW_I32_LOOP(3, {
      rsimd_mi32 bad = rsimd_mi32_andnot(rsimd_vi32_na2(b, c), rsimd_vi32_cmp_gt(b, c));
      if (rsimd_mi32_any(bad)) lohi = 1;
      r = rsimd_vi32_set_na(rsimd_vi32_max(a, b), rsimd_vi32_na2(a, b));
      r = rsimd_vi32_min(r, c);
    });
    break;
  default: break;
  }
  return (ovf ? RSIMD_EW_OVERFLOW : 0) | (lohi ? RSIMD_EW_LO_GT_HI : 0);
}

/* ---- double kernels ------------------------------------------------------ */

#ifndef RSIMD_NO_F64_SIMD

/* int32 elements converted to doubles, NA to NA_real_ (with check). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_ew_from_i32(rsimd_vf64 v, int check) {
  if (!check) return v;
  return rsimd_vf64_blend(v, rsimd_vf64_set1(rsimd_na_real()),
                          rsimd_vf64_cmp_eq(v, rsimd_vf64_set1(RSIMD_NA_I32_AS_F64)));
}
/* Operand k as a broadcast vector of its element 0. */
static inline rsimd_vf64 rsimd_ew_bcast(const void *p, int flags, int k, int check) {
  double v;
  if (p == NULL) return rsimd_vf64_zero();
  if (RSIMD_EW_IS_I32(flags, k)) {
    int e = ((const int *) p)[0];
    v = check && e == RSIMD_NA_I32 ? rsimd_na_real() : (double) e;
  } else {
    v = ((const double *) p)[0];
  }
  return rsimd_vf64_set1(v);
}
/* Operand k at i: full vector, or the predicated last vector whose
   inactive lanes repeat element i. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_ew_ld(const void *p, int flags, int k, rsimd_vf64 bc,
                                           ptrdiff_t i, int check) {
  if (RSIMD_EW_IS_SCALAR(flags, k)) return bc;
  if (RSIMD_EW_IS_I32(flags, k)) {
    return rsimd_ew_from_i32(rsimd_vf64_loadu_i32((const int32_t *) p + i), check);
  }
  return rsimd_vf64_loadu((const double *) p + i);
}
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_ew_ld_p(const void *p, int flags, int k, rsimd_vf64 bc,
                                             ptrdiff_t i, rsimd_p64 pg, int check) {
  if (RSIMD_EW_IS_SCALAR(flags, k)) return bc;
  if (RSIMD_EW_IS_I32(flags, k)) {
    const int32_t *q = (const int32_t *) p + i;
    return rsimd_ew_from_i32(rsimd_vf64_loadu_i32_p(pg, q, q[0]), check);
  }
  return rsimd_vf64_loadu_p(pg, (const double *) p + i, ((const double *) p)[i]);
}

/* 1 if every int32 operand is a scalar (whose broadcast value
   rsimd_ew_bcast() has converted), so that all operands can be read as
   doubles through rsimd_ew_dptr(). */
#define RSIMD_EW_NO_I32_VECTOR(f) ((((f) >> 3) & ~(f) & 7) == 0)
/* Operand k as doubles read at a pointer that the loop advances by
   *step after each vector: a double vector's data (step RSIMD_LANES_64),
   or buf holding the broadcast value bc (step 0) for a scalar or an absent
   operand, so that the loop reads every operand the same way, without a
   per-vector test of its flags. */
static inline const double *rsimd_ew_dptr(const void *p, int flags, int k, rsimd_vf64 bc,
                                          double *buf, ptrdiff_t *step) {
  if (p != NULL && !RSIMD_EW_IS_SCALAR(flags, k)) {
    *step = (ptrdiff_t) RSIMD_LANES_64;
    return (const double *) p;
  }
  rsimd_vf64_storeu(buf, bc);
  *step = 0;
  return buf;
}
#define RSIMD_EW_DLD(k) rsimd_vf64_loadu(p##k##_)
#define RSIMD_EW_DLD_P(k) rsimd_vf64_loadu_p(pg, p##k##_, p##k##_[0])
/* Declares the pointers and steps of rsimd_ew_dptr() for operands x, y
   and z (NULL when absent), with broadcast values bc0, bc1 and bc2, and
   RSIMD_EW_DNEXT() advances them. */
#define RSIMD_EW_DPTRS(x, y, z)                                                  \
  double b0_[RSIMD_MAX_LANES_64], b1_[RSIMD_MAX_LANES_64], b2_[RSIMD_MAX_LANES_64]; \
  ptrdiff_t s0_, s1_, s2_;                                                       \
  const double *p0_ = rsimd_ew_dptr(x, flags, 0, bc0, b0_, &s0_);               \
  const double *p1_ = rsimd_ew_dptr(y, flags, 1, bc1, b1_, &s1_);               \
  const double *p2_ = rsimd_ew_dptr(z, flags, 2, bc2, b2_, &s2_)
#define RSIMD_EW_DNEXT()                                                         \
  do {                                                                           \
    p0_ += s0_;                                                                  \
    p1_ += s1_;                                                                  \
    p2_ += s2_;                                                                  \
  } while (0)

/* Runs `expr`, which sets the rsimd_vf64 r from a, b (and c), over the
   chunk with operand flags fl; `lanes` is the number of active lanes. */
#define RSIMD_EW_F64_LOOP_(fl, nargs, expr)                                      \
  do {                                                                           \
    ptrdiff_t i = 0;                                                             \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                       \
      const int lanes = (int) RSIMD_LANES_64;                                    \
      rsimd_vf64 a = rsimd_ew_ld(x, fl, 0, bc0, i, check), b = bc1, c = bc2, r;  \
      if ((nargs) > 1) b = rsimd_ew_ld(y, fl, 1, bc1, i, check);                 \
      if ((nargs) > 2) c = rsimd_ew_ld(z, fl, 2, bc2, i, check);                 \
      (void) b;                                                                  \
      (void) c;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu(out + i, r);                                             \
    }                                                                            \
    if (i < n) {                                                                 \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                      \
      const int lanes = rsimd_p64_count(pg);                                     \
      rsimd_vf64 a = rsimd_ew_ld_p(x, fl, 0, bc0, i, pg, check), b = bc1, c = bc2, r; \
      if ((nargs) > 1) b = rsimd_ew_ld_p(y, fl, 1, bc1, i, pg, check);           \
      if ((nargs) > 2) c = rsimd_ew_ld_p(z, fl, 2, bc2, i, pg, check);           \
      (void) b;                                                                  \
      (void) c;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu_p(pg, out + i, r);                                       \
    }                                                                            \
  } while (0)
/* As RSIMD_EW_F64_LOOP_, reading the operands through rsimd_ew_dptr(). */
#define RSIMD_EW_F64_LOOP_D_(nargs, expr)                                        \
  do {                                                                           \
    ptrdiff_t i = 0;                                                             \
    RSIMD_EW_DPTRS(x, y, z);                                                     \
    (void) p1_;                                                                  \
    (void) p2_;                                                                  \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                       \
      const int lanes = (int) RSIMD_LANES_64;                                    \
      rsimd_vf64 a = RSIMD_EW_DLD(0), b = bc1, c = bc2, r;                       \
      if ((nargs) > 1) b = RSIMD_EW_DLD(1);                                      \
      if ((nargs) > 2) c = RSIMD_EW_DLD(2);                                      \
      (void) b;                                                                  \
      (void) c;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu(out + i, r);                                             \
      RSIMD_EW_DNEXT();                                                          \
    }                                                                            \
    if (i < n) {                                                                 \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                      \
      const int lanes = rsimd_p64_count(pg);                                     \
      rsimd_vf64 a = RSIMD_EW_DLD_P(0), b = bc1, c = bc2, r;                     \
      if ((nargs) > 1) b = RSIMD_EW_DLD_P(1);                                    \
      if ((nargs) > 2) c = RSIMD_EW_DLD_P(2);                                    \
      (void) b;                                                                  \
      (void) c;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu_p(pg, out + i, r);                                       \
    }                                                                            \
  } while (0)
/* The loop with the flags as constants in the two common cases (all
   operands double vectors; a double scalar y), the loop over double
   operands, vectors or scalars, without per-operand tests otherwise, and
   the general one for int32 vectors. The first two are faster still, as
   they keep a scalar in a register and index all vectors with i. */
#define RSIMD_EW_F64_LOOP(nargs, expr)                                           \
  do {                                                                           \
    if (flags == 0) {                                                            \
      RSIMD_EW_F64_LOOP_(0, nargs, expr);                                        \
    } else if (flags == RSIMD_EW_SCALAR(1)) {                                    \
      RSIMD_EW_F64_LOOP_(RSIMD_EW_SCALAR(1), nargs, expr);                       \
    } else if (RSIMD_EW_NO_I32_VECTOR(flags)) {                                  \
      RSIMD_EW_F64_LOOP_D_(nargs, expr);                                         \
    } else {                                                                     \
      RSIMD_EW_F64_LOOP_(flags, nargs, expr);                                    \
    }                                                                            \
  } while (0)

/* x %% y: the vector form, with the lanes it leaves to the scalar form
   (huge quotients) recomputed. */
static inline rsimd_vf64 rsimd_ew_mod(rsimd_vf64 a, rsimd_vf64 b, int lanes) {
  rsimd_mf64 slow;
  rsimd_vf64 r = rsimd_vf64_mod(a, b, &slow);
  if (rsimd_mf64_any(slow)) {
    double ba[RSIMD_MAX_LANES_64], bb[RSIMD_MAX_LANES_64], br[RSIMD_MAX_LANES_64],
      bs[RSIMD_MAX_LANES_64];
    int j;
    rsimd_vf64_storeu(ba, a);
    rsimd_vf64_storeu(bb, b);
    rsimd_vf64_storeu(br, r);
    rsimd_vf64_storeu(bs, rsimd_vf64_blend(rsimd_vf64_zero(), rsimd_vf64_set1(1.0), slow));
    for (j = 0; j < lanes; j++) {
      if (bs[j] != 0) br[j] = rsimd_mod_f64(ba[j], bb[j]);
    }
    r = rsimd_vf64_loadu(br);
  }
  return r;
}

int RSIMD_KERNEL(ew1_f64)(int op, const void *x, R_xlen_t n, int flags, double *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew1_f64)(int op, const void *x, R_xlen_t n, int flags, double *out,
                          const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int check = o->na_check;
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, check), bc1 = rsimd_vf64_zero(), bc2 = bc1;
  int st = 0;
  (void) y;
  (void) z;
  switch (op) {
  case RSIMD_EW_NEG: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_neg(a)); break;
  case RSIMD_EW_ABS: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_abs(a)); break;
  case RSIMD_EW_SIGN: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_sign(a)); break;
  case RSIMD_EW_RECIP: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_div(rsimd_vf64_set1(1.0), a)); break;
  case RSIMD_EW_SQRT:
    RSIMD_EW_F64_LOOP(1, {
      if (rsimd_mf64_any(rsimd_vf64_cmp_lt(a, rsimd_vf64_zero()))) st |= RSIMD_EW_NAN_PRODUCED;
      r = rsimd_vf64_sqrt(a);
    });
    break;
  case RSIMD_EW_FLOOR: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_floor(a)); break;
  case RSIMD_EW_CEIL: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_ceil(a)); break;
  case RSIMD_EW_TRUNC: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_trunc(a)); break;
  case RSIMD_EW_ROUND: RSIMD_EW_F64_LOOP(1, r = rsimd_vf64_rint(a)); break;
  default: break;
  }
  return st;
}

/* rsimd_round_digits_f64() lane by lane. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_ew_round_digits(rsimd_vf64 x, rsimd_vf64 p10, rsimd_vf64 big) {
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_vf64 a = rsimd_vf64_abs(x), x10 = rsimd_vf64_mul(a, p10), i10 = rsimd_vf64_floor(x10),
             h = rsimd_vf64_mul(i10, rsimd_vf64_set1(0.5)), xd = rsimd_vf64_div(i10, p10),
             xu = rsimd_vf64_div(rsimd_vf64_ceil(x10), p10), du = rsimd_vf64_sub(xu, a),
             dd = rsimd_vf64_sub(a, xd), r;
  /* i10 is odd when i10 / 2 is not an integer. */
  rsimd_mf64 odd = rsimd_vf64_cmp_lt(rsimd_vf64_floor(h), h),
             up = rsimd_mf64_or(rsimd_vf64_cmp_lt(du, dd),
                                rsimd_mf64_and(odd, rsimd_vf64_cmp_eq(du, dd))),
             live = rsimd_mf64_and(rsimd_vf64_cmp_gt(a, zero), rsimd_vf64_cmp_lt(a, big));
  r = rsimd_vf64_blend(xd, xu, up);
  r = rsimd_vf64_blend(r, rsimd_vf64_neg(r), rsimd_vf64_cmp_lt(x, zero));
  return rsimd_vf64_blend(x, r, live);
}

void RSIMD_KERNEL(round_digits_f64)(const void *x, R_xlen_t n, int flags, double p10, double big,
                                    double *out, const rsimd_opts *o);
void RSIMD_KERNEL(round_digits_f64)(const void *x, R_xlen_t n, int flags, double p10, double big,
                                    double *out, const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int check = o->na_check;
  /* b and c of the loop are the broadcasts bc1 = p10 and bc2 = big. */
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, check), bc1 = rsimd_vf64_set1(p10),
                   bc2 = rsimd_vf64_set1(big);
  (void) y;
  (void) z;
  RSIMD_EW_F64_LOOP(1, r = rsimd_ew_round_digits(a, b, c));
}

int RSIMD_KERNEL(ew2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          double *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          double *out, const rsimd_opts *o) {
  const void *z = NULL;
  const int check = o->na_check;
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, check), bc1 = rsimd_ew_bcast(y, flags, 1, check),
                   bc2 = rsimd_vf64_zero();
  (void) z;
/* A NaN operand makes these results NaN, so only a vector with a NaN
   result can need NA_real_ put back. */
#define RSIMD_EW_MERGED(e)                                                       \
  do {                                                                           \
    r = (e);                                                                     \
    if (check && rsimd_mf64_any(rsimd_vf64_is_nan(r))) r = rsimd_vf64_na_merge(r, a, b); \
  } while (0)
  switch (op) {
  case RSIMD_EW_ADD: RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_vf64_add(a, b))); break;
  case RSIMD_EW_SUB: RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_vf64_sub(a, b))); break;
  case RSIMD_EW_MUL: RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_vf64_mul(a, b))); break;
  case RSIMD_EW_DIV: RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_vf64_div(a, b))); break;
  case RSIMD_EW_IDIV: {
    /* The vector form, then the elements it leaves to the scalar form
       (quotients of 2^52 and more, which are rare) recomputed in a second
       pass, so that the loop has no per-vector branch for them. */
    rsimd_mf64 big = rsimd_vf64_cmp_lt(bc2, bc2), s;
    R_xlen_t k;
    RSIMD_EW_F64_LOOP(2, {
      RSIMD_EW_MERGED(rsimd_vf64_idiv(a, b, &s));
      big = rsimd_mf64_or(big, s);
    });
    if (rsimd_mf64_any(big)) {
      for (k = 0; k < n; k++) {
        double a = rsimd_ew_get(x, flags, 0, k, check), b = rsimd_ew_get(y, flags, 1, k, check);
        if (fabs(a / b) >= 0x1p52) out[k] = rsimd_idiv_f64(a, b);
      }
    }
    break;
  }
  case RSIMD_EW_MOD: RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_ew_mod(a, b, lanes))); break;
  case RSIMD_EW_COPYSIGN:
    RSIMD_EW_F64_LOOP(2, RSIMD_EW_MERGED(rsimd_vf64_copysign(a, b)));
    break;
  case RSIMD_EW_PMIN: RSIMD_EW_F64_LOOP(2, r = rsimd_vf64_pminmax(a, b, 0, 0)); break;
  case RSIMD_EW_PMAX: RSIMD_EW_F64_LOOP(2, r = rsimd_vf64_pminmax(a, b, 1, 0)); break;
  case RSIMD_EW_PMIN_NUM: RSIMD_EW_F64_LOOP(2, r = rsimd_vf64_pminmax(a, b, 0, 1)); break;
  case RSIMD_EW_PMAX_NUM: RSIMD_EW_F64_LOOP(2, r = rsimd_vf64_pminmax(a, b, 1, 1)); break;
  default: break;
  }
#undef RSIMD_EW_MERGED
  return 0;
}

int RSIMD_KERNEL(ew3_f64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, double *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_f64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, double *out, const rsimd_opts *o) {
  const int check = o->na_check;
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, check), bc1 = rsimd_ew_bcast(y, flags, 1, check),
                   bc2 = rsimd_ew_bcast(z, flags, 2, check);
  int st = 0;
#define RSIMD_EW_MERGED(e)                                                       \
  do {                                                                           \
    r = (e);                                                                     \
    if (check && rsimd_mf64_any(rsimd_vf64_is_nan(r))) {                         \
      r = rsimd_vf64_na_merge3(r, a, b, c);                                      \
    }                                                                            \
  } while (0)
  switch (op) {
  case RSIMD_EW_FMA: RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_fma(a, b, c))); break;
  case RSIMD_EW_MUL_ADD:
    RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_add(rsimd_vf64_mul(a, b), c)));
    break;
  case RSIMD_EW_MUL_ADD_APPROX:
#if RSIMD_NATIVE_FMA
    RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_fma(a, b, c)));
#else
    RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_add(rsimd_vf64_mul(a, b), c)));
#endif
    break;
  case RSIMD_EW_ADD_MUL:
    RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_mul(rsimd_vf64_add(a, b), c)));
    break;
  case RSIMD_EW_LERP:
    RSIMD_EW_F64_LOOP(3, RSIMD_EW_MERGED(rsimd_vf64_fma(
                             c, b, rsimd_vf64_mul(rsimd_vf64_sub(rsimd_vf64_set1(1.0), c), a))));
    break;
  case RSIMD_EW_CLAMP:
    RSIMD_EW_F64_LOOP(3, {
      if (rsimd_mf64_any(rsimd_vf64_cmp_gt(b, c))) st |= RSIMD_EW_LO_GT_HI;
      r = rsimd_vf64_pminmax(rsimd_vf64_pminmax(a, b, 1, 0), c, 0, 0);
    });
    break;
  default: break;
  }
#undef RSIMD_EW_MERGED
  return st;
}

#undef RSIMD_EW_F64_LOOP
#undef RSIMD_EW_F64_LOOP_
#undef RSIMD_EW_F64_LOOP_D_

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does doubles */
#define RSIMD_SKIP_ew1_f64 1
#define RSIMD_SKIP_ew2_f64 1
#define RSIMD_SKIP_ew3_f64 1
#define RSIMD_SKIP_round_digits_f64 1
#endif

#undef RSIMD_EW_I32_LOOP
#undef RSIMD_EW_I32_CHECKED
#undef RSIMD_EW_I32_WRAP
#undef RSIMD_EW_LDI
#undef RSIMD_EW_LDI_P

#endif /* vector tiers */

#undef RSIMD_EW_IS_SCALAR
#undef RSIMD_EW_IS_I32
