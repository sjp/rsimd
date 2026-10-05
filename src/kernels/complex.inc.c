/* Complex kernels: conjugate, real and imaginary parts, sum, the
 * predicates (is.na, is.nan, is.finite, is.infinite) and the missing-value
 * scans behind any_na, count_na and which_na, on R's interleaved Rcomplex
 * {r, i} layout. n complex numbers are 2n doubles, re, im, re, im ...
 *
 * Addition, subtraction and negation have no kernels here: the entry
 * points run the f64 kernels of arith.inc.c on the 2n doubles, so every
 * part follows the double rules exactly.
 *
 * The none tier is plain C and the reference the vector tiers are tested
 * against. The vector tiers load two vectors (W complex numbers, W the
 * number of 64-bit lanes, always even) and split them into the real and
 * imaginary parts with the layer's uzp_even/uzp_odd; conj stays
 * interleaved and flips the sign bit of the odd lanes.
 *
 * Included after reduce.inc.c (sum_f64, na_f64, put_index_) and
 * predicates.inc.c (rsimd_pred_f64_1, rsimd_pred_vf64).
 */

#include "logical.inc.h"

/* The predicate of one complex number (op RSIMD_PRED_NA, NAN, FINITE or
   INFINITE), as base R: NA if either part is NA or NaN, NaN if either part
   is a NaN that is not NA, finite if both parts are finite, infinite if
   either part is. The none tier's kernels and the reference of the
   others. */
static inline int rsimd_pred_c128_1(int op, Rcomplex z) {
  int a = rsimd_pred_f64_1(op, z.r), b = rsimd_pred_f64_1(op, z.i);
  return op == RSIMD_PRED_FINITE ? a && b : a || b;
}

/* The scalar Hamming loop over complex elements [i, n) (flags
   RSIMD_EW_SCALAR(k)): the none tier's kernel and the vector tiers' tail.
   Returns 1 when it stops at a missing pair (without na_rm). */
static inline int rsimd_ham_c128_from(const Rcomplex *x, const Rcomplex *y, R_xlen_t i, R_xlen_t n,
                                      int flags, rsimd_reduce_result *r, int na_rm) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  for (; i < n; i++) {
    Rcomplex a = x[i * sx], b = y[i * sy];
    if (isnan(a.r) || isnan(a.i) || isnan(b.r) || isnan(b.i)) {
      r->saw_na = 1;
      if (!na_rm) return 1;
      continue;
    }
    r->i64 += a.r != b.r || a.i != b.i;
  }
  return 0;
}

/* Complex numbers per block of sum_c128: the parts of a block are split
   into two stack buffers and folded by sum_f64. A multiple of
   RSIMD_PAIRWISE_LEAF, so pairwise leaves line up with those of a double
   vector of the parts. */
#define RSIMD_CPLX_BLOCK 512

/* ---- Multiply, divide, product and scans: scalar helpers ------------------ */

/* The real and imaginary parts of (a + bi)(c + di) in the rounding
   variant v (RSIMD_CMUL_UNFUSED, FMA1 or FMA2; see kernel_types.h). The
   tier units are compiled without contraction, so the unfused forms round
   every product. */
static inline double rsimd_cmul_re_(int v, double a, double b, double c, double d) {
  switch (v) {
  case RSIMD_CMUL_FMA1: return rsimd_fma(a, c, -(b * d));
  case RSIMD_CMUL_FMA2: return rsimd_fma(-b, d, a * c);
  default: return a * c - b * d;
  }
}
static inline double rsimd_cmul_im_(int v, double a, double b, double c, double d) {
  switch (v) {
  case RSIMD_CMUL_FMA1: return rsimd_fma(b, c, a * d);
  case RSIMD_CMUL_FMA2: return rsimd_fma(a, d, b * c);
  default: return a * d + b * c;
  }
}

/* x * y as base R computes it: the selected variant, and base R's own
   operator (mul1) where that gives a NaN part. */
static inline Rcomplex rsimd_cmul_1_(const rsimd_c128_arith *a, Rcomplex x, Rcomplex y) {
  Rcomplex r;
  if (a->mul_re != RSIMD_CMUL_SCALAR) {
    r.r = rsimd_cmul_re_(a->mul_re, x.r, x.i, y.r, y.i);
    r.i = rsimd_cmul_im_(a->mul_im, x.r, x.i, y.r, y.i);
    if (!isnan(r.r) && !isnan(r.i)) return r;
  }
  a->mul1(&x, &y, &r);
  return r;
}

/* (a + bi) / (c + di) by libgcc's __divdc3 (GCC 12 and later) up to its
   recovery of infinities, with the multiply-adds fused or rounded twice.
   This is the branch-free form the vector tiers use, which gives the same
   results as libgcc's branches. */
static inline void rsimd_cdiv_formula_(int fused, double a, double b, double c, double d, Rcomplex *out) {
  const int sw = fabs(c) < fabs(d);
  const double abig = sw ? fabs(d) : fabs(c), rmax2 = DBL_MAX / 2 * DBL_EPSILON;
  double big = sw ? d : c, small = sw ? c : d, f = 1.0, u, v, s, t, ratio, denom, x, y;
  if (abig >= DBL_MAX / 2) f = 0.5;
  else if (abig < DBL_EPSILON ||
           (abig < rmax2 && ((fabs(a) < DBL_MIN && fabs(b) < rmax2) ||
                             (fabs(b) < DBL_MIN && fabs(a) < rmax2)))) {
    f = 1.0 / DBL_EPSILON;
  }
  a *= f;
  b *= f;
  big *= f;
  small *= f;
  u = sw ? a : b;
  v = sw ? b : a;
  s = sw ? b : -a;
  t = sw ? -a : b;
  ratio = small / big;
  denom = fused ? rsimd_fma(small, ratio, big) : small * ratio + big;
  if (fabs(ratio) > DBL_MIN) {
    x = fused ? rsimd_fma(u, ratio, v) : u * ratio + v;
    y = fused ? rsimd_fma(s, ratio, t) : s * ratio + t;
  } else {
    x = fused ? rsimd_fma(small, u / big, v) : small * (u / big) + v;
    y = fused ? rsimd_fma(small, s / big, t) : small * (s / big) + t;
  }
  out->r = x / denom;
  out->i = y / denom;
}

/* x / y as base R computes it: the selected variant of libgcc's
   algorithm, and base R's own operator (div1) where that gives a NaN
   part. */
static inline Rcomplex rsimd_cdiv_1_(const rsimd_c128_arith *a, Rcomplex x, Rcomplex y) {
  Rcomplex r;
  if (a->div != RSIMD_CDIV_SCALAR) {
    rsimd_cdiv_formula_(a->div == RSIMD_CDIV_FMA, x.r, x.i, y.r, y.i, &r);
    if (!isnan(r.r) && !isnan(r.i)) return r;
  }
  a->div1(&x, &y, &r);
  return r;
}

/* x op y for elements [i, n) one by one (the none tier's kernel and the
   vector tiers' tail). */
static inline void rsimd_ew2_c128_from(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t i,
                                       R_xlen_t n, int flags, Rcomplex *out,
                                       const rsimd_c128_arith *a) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  if (op == RSIMD_EW_MUL) {
    for (; i < n; i++) out[i] = rsimd_cmul_1_(a, x[i * sx], y[i * sy]);
  } else {
    for (; i < n; i++) out[i] = rsimd_cdiv_1_(a, x[i * sx], y[i * sy]);
  }
}

/* An element the product skips: either part NA or NaN when missing values
   are checked, which sets the flags. */
static inline int rsimd_cprod_skip_(Rcomplex e, rsimd_cprod_state *s, const rsimd_opts *o) {
  if (!(o->na_check || o->na_rm) || !(isnan(e.r) || isnan(e.i))) return 0;
  s->saw_nan = 1;
  if (rsimd_is_na_f64(e.r) || rsimd_is_na_f64(e.i)) s->saw_na = 1;
  return 1;
}

/* s->p (+ s->lo) times e in double-double arithmetic: the error-free
   products (by fma) and sums of both parts are kept in s->lo. A
   non-finite product continues in plain arithmetic with lo = 0. */
static inline void rsimd_cprod_dd_(const rsimd_c128_arith *a, rsimd_cprod_state *s, Rcomplex e) {
  const Rcomplex p = s->p, lo = s->lo, plain = rsimd_cmul_1_(a, p, e);
  double t1, t2, t3, t4, sr, si, z, lr, li, hr, hi;
  if (isfinite(plain.r) && isfinite(plain.i)) {
    t1 = p.r * e.r;
    t2 = p.i * e.i;
    sr = t1 - t2;
    z = sr - t1;
    lr = ((t1 - (sr - z)) + (-t2 - z)) + (rsimd_fma(p.r, e.r, -t1) - rsimd_fma(p.i, e.i, -t2)) +
         (lo.r * e.r - lo.i * e.i);
    t3 = p.r * e.i;
    t4 = p.i * e.r;
    si = t3 + t4;
    z = si - t3;
    li = ((t3 - (si - z)) + (t4 - z)) + (rsimd_fma(p.r, e.i, -t3) + rsimd_fma(p.i, e.r, -t4)) +
         (lo.r * e.i + lo.i * e.r);
    hr = sr + lr;
    hi = si + li;
    if (isfinite(hr) && isfinite(hi)) {
      s->p.r = hr;
      s->p.i = hi;
      s->lo.r = lr - (hr - sr);
      s->lo.i = li - (hi - si);
      return;
    }
  }
  s->p = plain;
  s->lo.r = s->lo.i = 0.0;
}

/* s->p times e in pairwise mode: a full leaf is merged into the leaf
   products like a binary counter, so the merge order depends only on the
   element's position, not on the chunks. */
static inline void rsimd_cprod_pw_(const rsimd_c128_arith *a, rsimd_cprod_state *s, Rcomplex e) {
  Rcomplex q;
  int k = 0;
  s->p = rsimd_cmul_1_(a, s->p, e);
  if (++s->fill < RSIMD_PAIRWISE_LEAF) return;
  q = s->p;
  while ((s->leaves >> k) & 1) {
    q = rsimd_cmul_1_(a, s->pw[k], q);
    s->leaves &= ~((uint64_t) 1 << k);
    k++;
  }
  s->pw[k] = q;
  s->leaves |= (uint64_t) 1 << k;
  s->p.r = 1.0;
  s->p.i = 0.0;
  s->fill = 0;
}

/* The product of elements [i, n) in order (every mode of the none tier,
   pairwise and compensated mode and the fast tail of the others). */
static inline void rsimd_cprod_seq_(const Rcomplex *x, R_xlen_t i, R_xlen_t n, rsimd_cprod_state *s,
                                    const rsimd_c128_arith *a, const rsimd_opts *o) {
  for (; i < n; i++) {
    if (rsimd_cprod_skip_(x[i], s, o)) continue;
    if (o->precision == RSIMD_PREC_COMPENSATED) rsimd_cprod_dd_(a, s, x[i]);
    else if (o->precision == RSIMD_PREC_PAIRWISE) rsimd_cprod_pw_(a, s, x[i]);
    else s->p = rsimd_cmul_1_(a, s->p, x[i]);
  }
}

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(ew2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            Rcomplex *out, const rsimd_c128_arith *a);
void RSIMD_KERNEL(ew2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            Rcomplex *out, const rsimd_c128_arith *a) {
  rsimd_ew2_c128_from(op, x, y, 0, n, flags, out, a);
}

void RSIMD_KERNEL(prod_c128)(const Rcomplex *x, R_xlen_t n, rsimd_cprod_state *s,
                             const rsimd_c128_arith *a, const rsimd_opts *o);
void RSIMD_KERNEL(prod_c128)(const Rcomplex *x, R_xlen_t n, rsimd_cprod_state *s,
                             const rsimd_c128_arith *a, const rsimd_opts *o) {
  rsimd_cprod_seq_(x, 0, n, s, a, o);
}

/* Sequential on every tier: each element depends on the one before. */
void RSIMD_KERNEL(scan_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out, Rcomplex *s,
                             const rsimd_c128_arith *a);
void RSIMD_KERNEL(scan_c128)(int op, const Rcomplex *x, R_xlen_t n, Rcomplex *out, Rcomplex *s,
                             const rsimd_c128_arith *a) {
  Rcomplex t = *s;
  R_xlen_t i;
  if (op == 0) {
    for (i = 0; i < n; i++) {
      t.r += x[i].r;
      t.i += x[i].i;
      out[i] = t;
    }
  } else if (a->cp_re == RSIMD_CMUL_SCALAR) {
    for (i = 0; i < n; i++) {
      Rcomplex u = t;
      a->cp1(x + i, &u, &t);
      out[i] = t;
    }
  } else {
    for (i = 0; i < n; i++) {
      double re = rsimd_cmul_re_(a->cp_re, x[i].r, x[i].i, t.r, t.i);
      t.i = rsimd_cmul_im_(a->cp_im, x[i].r, x[i].i, t.r, t.i);
      t.r = re;
      out[i] = t;
    }
  }
  *s = t;
}

void RSIMD_KERNEL(formula_c128)(int op, int v1, int v2, const Rcomplex *x, const Rcomplex *y,
                                R_xlen_t n, Rcomplex *out);
void RSIMD_KERNEL(formula_c128)(int op, int v1, int v2, const Rcomplex *x, const Rcomplex *y,
                                R_xlen_t n, Rcomplex *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (op == RSIMD_EW_DIV) {
      rsimd_cdiv_formula_(v1 == RSIMD_CDIV_FMA, x[i].r, x[i].i, y[i].r, y[i].i, out + i);
    } else {
      double re = rsimd_cmul_re_(v1, x[i].r, x[i].i, y[i].r, y[i].i);
      out[i].i = rsimd_cmul_im_(v2, x[i].r, x[i].i, y[i].r, y[i].i);
      out[i].r = re;
    }
  }
}

void RSIMD_KERNEL(cmp_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            int *out);
void RSIMD_KERNEL(cmp_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            int *out) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    Rcomplex u = x[i * sx], v = y[i * sy];
    if (isnan(u.r) || isnan(u.i) || isnan(v.r) || isnan(v.i)) out[i] = RSIMD_NA_I32;
    else out[i] = (u.r == v.r && u.i == v.i) == (op == RSIMD_CMP_EQ);
  }
}

void RSIMD_KERNEL(conj_c128)(const Rcomplex *x, R_xlen_t n, Rcomplex *out);
void RSIMD_KERNEL(conj_c128)(const Rcomplex *x, R_xlen_t n, Rcomplex *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    out[i].r = x[i].r;
    out[i].i = -x[i].i;
  }
}

void RSIMD_KERNEL(part_c128)(int im, const Rcomplex *x, R_xlen_t n, double *out);
void RSIMD_KERNEL(part_c128)(int im, const Rcomplex *x, R_xlen_t n, double *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) out[i] = im ? x[i].i : x[i].r;
}

/* The parts of the n complex numbers of x into re and im; with narm, both
   parts of an element with a missing part become 0. Returns the number of
   those elements. */
static ptrdiff_t RSIMD_KERNEL(split_c128_)(const Rcomplex *x, ptrdiff_t n, double *re, double *im,
                                           int narm) {
  ptrdiff_t i, removed = 0;
  for (i = 0; i < n; i++) {
    re[i] = x[i].r;
    im[i] = x[i].i;
    if (narm && (isnan(re[i]) || isnan(im[i]))) {
      re[i] = im[i] = 0.0;
      removed++;
    }
  }
  return removed;
}

int RSIMD_KERNEL(pred_c128)(int op, const Rcomplex *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_c128)(int op, const Rcomplex *x, R_xlen_t n, int mode, int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int v = rsimd_pred_c128_1(op, x[i]);
    if (mode == RSIMD_PRED_ELT) out[i] = v;
    else if (mode == RSIMD_PRED_ANY && v) return 1;
    else if (mode == RSIMD_PRED_ALL && !v) return 0;
  }
  return mode == RSIMD_PRED_ALL;
}

void RSIMD_KERNEL(hamming_c128)(const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                                rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_c128)(const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                                rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_ham_c128_from(x, y, 0, n, flags, r, o->na_rm);
}

void RSIMD_KERNEL(na_c128)(const Rcomplex *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                           rsimd_reduce_result *r);
void RSIMD_KERNEL(na_c128)(const Rcomplex *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                           rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (!(isnan(x[i].r) || isnan(x[i].i))) continue;
    if (mode == RSIMD_NAMODE_ANY) {
      r->saw_na = 1;
      return;
    }
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_KERNEL(put_index_)(mode, out, r, off + i);
  }
}

#elif !defined(RSIMD_NO_F64_SIMD) /* vector tiers */

/* The real and imaginary parts of the W complex numbers from x[i]. */
RSIMD_ALWAYS_INLINE void rsimd_c128_load_(const Rcomplex *x, ptrdiff_t i, rsimd_vf64 *re,
                                          rsimd_vf64 *im) {
  const double *d = (const double *) x + 2 * i;
  rsimd_vf64 a = rsimd_vf64_loadu(d), b = rsimd_vf64_loadu(d + RSIMD_LANES_64);
  *re = rsimd_vf64_uzp_even(a, b);
  *im = rsimd_vf64_uzp_odd(a, b);
}
/* As rsimd_c128_load_() for the last n - i < W complex numbers; both parts
   of the inactive lanes are `fill`. */
RSIMD_ALWAYS_INLINE void rsimd_c128_load_tail_(const Rcomplex *x, ptrdiff_t i, ptrdiff_t n,
                                               double fill, rsimd_vf64 *re, rsimd_vf64 *im) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const double *d = (const double *) x + 2 * i;
  rsimd_vf64 a = rsimd_vf64_loadu_p(rsimd_p64_while(2 * i, 2 * n), d, fill);
  rsimd_vf64 b = rsimd_vf64_loadu_p(rsimd_p64_while(2 * i + W, 2 * n), d + W, fill);
  *re = rsimd_vf64_uzp_even(a, b);
  *im = rsimd_vf64_uzp_odd(a, b);
}

void RSIMD_KERNEL(conj_c128)(const Rcomplex *x, R_xlen_t n, Rcomplex *out);
void RSIMD_KERNEL(conj_c128)(const Rcomplex *x, R_xlen_t n, Rcomplex *out) {
  const ptrdiff_t W = RSIMD_LANES_64, m = 2 * (ptrdiff_t) n;
  const double *d = (const double *) x;
  double *o = (double *) out, b[RSIMD_MAX_LANES_64];
  ptrdiff_t i = 0;
  rsimd_vf64 sign;
  /* -0.0 in the odd (imaginary) lanes; every vector starts at an even
     index because W is even. */
  for (i = 0; i < W; i++) b[i] = i & 1 ? -0.0 : 0.0;
  sign = rsimd_vf64_loadu(b);
  for (i = 0; i + W <= m; i += W) rsimd_vf64_storeu(o + i, rsimd_vf64_xor(rsimd_vf64_loadu(d + i), sign));
  if (i < m) {
    rsimd_p64 pg = rsimd_p64_while(i, m);
    rsimd_vf64_storeu_p(pg, o + i, rsimd_vf64_xor(rsimd_vf64_loadu_p(pg, d + i, 0.0), sign));
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(part_c128_)(const int im, const Rcomplex *x, R_xlen_t n,
                                                  double *out) {
  const ptrdiff_t W = RSIMD_LANES_64;
  ptrdiff_t i = 0;
  rsimd_vf64 re, iv;
  for (; i + W <= n; i += W) {
    rsimd_c128_load_(x, i, &re, &iv);
    rsimd_vf64_storeu(out + i, im ? iv : re);
  }
  if (i < n) {
    rsimd_c128_load_tail_(x, i, n, 0.0, &re, &iv);
    rsimd_vf64_storeu_p(rsimd_p64_while(i, n), out + i, im ? iv : re);
  }
}

void RSIMD_KERNEL(part_c128)(int im, const Rcomplex *x, R_xlen_t n, double *out);
void RSIMD_KERNEL(part_c128)(int im, const Rcomplex *x, R_xlen_t n, double *out) {
  if (im) RSIMD_KERNEL(part_c128_)(1, x, n, out);
  else RSIMD_KERNEL(part_c128_)(0, x, n, out);
}

RSIMD_ALWAYS_INLINE ptrdiff_t RSIMD_KERNEL(split_c128_v_)(const Rcomplex *x, ptrdiff_t n,
                                                          double *re, double *im,
                                                          const int narm) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 zero = rsimd_vf64_zero();
  ptrdiff_t i = 0, removed = 0;
  rsimd_vf64 vr, vi;
  rsimd_mf64 m;
#define RSIMD_CPLX_DROP_                                                                   \
  if (narm) {                                                                              \
    m = rsimd_mf64_or(rsimd_vf64_is_nan(vr), rsimd_vf64_is_nan(vi));                      \
    vr = rsimd_vf64_blend(vr, zero, m);                                                    \
    vi = rsimd_vf64_blend(vi, zero, m);                                                    \
    removed += rsimd_mf64_count(m);                                                        \
  }
  for (; i + W <= n; i += W) {
    rsimd_c128_load_(x, i, &vr, &vi);
    RSIMD_CPLX_DROP_
    rsimd_vf64_storeu(re + i, vr);
    rsimd_vf64_storeu(im + i, vi);
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    /* The inactive lanes are 0, so they are never counted as removed. */
    rsimd_c128_load_tail_(x, i, n, 0.0, &vr, &vi);
    RSIMD_CPLX_DROP_
    rsimd_vf64_storeu_p(pg, re + i, vr);
    rsimd_vf64_storeu_p(pg, im + i, vi);
  }
#undef RSIMD_CPLX_DROP_
  return removed;
}

/* As the none tier's split_c128_. */
static ptrdiff_t RSIMD_KERNEL(split_c128_)(const Rcomplex *x, ptrdiff_t n, double *re, double *im,
                                           int narm) {
  return narm ? RSIMD_KERNEL(split_c128_v_)(x, n, re, im, 1)
              : RSIMD_KERNEL(split_c128_v_)(x, n, re, im, 0);
}

RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_pred_vc128(const int op, rsimd_vf64 re, rsimd_vf64 im) {
  rsimd_mf64 a = rsimd_pred_vf64(op, re), b = rsimd_pred_vf64(op, im);
  return op == RSIMD_PRED_FINITE ? rsimd_mf64_and(a, b) : rsimd_mf64_or(a, b);
}

RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_c128_)(const int op, const Rcomplex *x, R_xlen_t n,
                                                 int mode, int32_t *out) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_mf64 none = rsimd_mf64_none();
  ptrdiff_t i = 0;
  rsimd_vf64 re, im;
  rsimd_mf64 m;
  for (; i + W <= n; i += W) {
    rsimd_c128_load_(x, i, &re, &im);
    m = rsimd_pred_vc128(op, re, im);
    if (mode == RSIMD_PRED_ELT) {
      rsimd_vi64_storeu_i32(out + i, rsimd_lgl_vi64(m, none));
    } else if (mode == RSIMD_PRED_ANY) {
      if (rsimd_mf64_any(m)) return 1;
    } else if (!rsimd_mf64_all(m)) {
      return 0;
    }
  }
  if (i < n) {
    if (mode == RSIMD_PRED_ELT) {
      rsimd_c128_load_tail_(x, i, n, 0.0, &re, &im);
      m = rsimd_pred_vc128(op, re, im);
      rsimd_vi64_storeu_i32_p(rsimd_p64_while(i, n), out + i, rsimd_lgl_vi64(m, none));
    } else {
      /* Fewer than W elements: any/all of the tail one by one, so no
         inactive lane takes part. */
      for (; i < n; i++) {
        int v = rsimd_pred_c128_1(op, x[i]);
        if (mode == RSIMD_PRED_ANY && v) return 1;
        if (mode == RSIMD_PRED_ALL && !v) return 0;
      }
    }
  }
  return mode == RSIMD_PRED_ALL;
}

int RSIMD_KERNEL(pred_c128)(int op, const Rcomplex *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_c128)(int op, const Rcomplex *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_c128_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_NAN: return RSIMD_KERNEL(pred_c128_)(RSIMD_PRED_NAN, x, n, mode, out);
  case RSIMD_PRED_FINITE: return RSIMD_KERNEL(pred_c128_)(RSIMD_PRED_FINITE, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_c128_)(RSIMD_PRED_INFINITE, x, n, mode, out);
  }
}

/* A pair differs where cmp_ne is true for either part, which includes
   every pair with a NaN part, so the differing non-missing pairs are the
   ne lanes less the missing ones. Lane counters are folded every
   RSIMD_HAMMING_BLOCK vectors. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(hamming_c128_)(const Rcomplex *x, const Rcomplex *y,
                                                     R_xlen_t n, int flags,
                                                     rsimd_reduce_result *r, int na_rm) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  const rsimd_vf64 sxr = rsimd_vf64_set1(x[0].r), sxi = rsimd_vf64_set1(x[0].i),
                   syr = rsimd_vf64_set1(y[0].r), syi = rsimd_vf64_set1(y[0].i);
  ptrdiff_t i = 0;
  while (i + W <= n) {
    const ptrdiff_t end = (n - i) / RSIMD_LANES_64 > RSIMD_HAMMING_BLOCK
                              ? i + RSIMD_HAMMING_BLOCK * RSIMD_LANES_64
                              : n;
    rsimd_vi64 ne = rsimd_vi64_zero(), m = rsimd_vi64_zero();
    int64_t missing;
    for (; i + W <= end; i += W) {
      rsimd_vf64 xr = sxr, xi = sxi, yr = syr, yi = syi;
      rsimd_mf64 nan;
      if (!sx) rsimd_c128_load_(x, i, &xr, &xi);
      if (!sy) rsimd_c128_load_(y, i, &yr, &yi);
      nan = rsimd_mf64_or(rsimd_mf64_or(rsimd_vf64_is_nan(xr), rsimd_vf64_is_nan(xi)),
                          rsimd_mf64_or(rsimd_vf64_is_nan(yr), rsimd_vf64_is_nan(yi)));
      ne = rsimd_vi64_inc(ne, rsimd_mf64_to_mi64(rsimd_mf64_or(rsimd_vf64_cmp_ne(xr, yr),
                                                               rsimd_vf64_cmp_ne(xi, yi))));
      m = rsimd_vi64_inc(m, rsimd_mf64_to_mi64(nan));
    }
    missing = rsimd_vi64_reduce_add(m);
    r->i64 += rsimd_vi64_reduce_add(ne) - missing;
    if (missing > 0) {
      r->saw_na = 1;
      if (!na_rm) return;
    }
  }
  rsimd_ham_c128_from(x, y, i, n, flags, r, na_rm);
}

void RSIMD_KERNEL(hamming_c128)(const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                                rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_c128)(const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                                rsimd_reduce_result *r, const rsimd_opts *o) {
  if (flags == 0) RSIMD_KERNEL(hamming_c128_)(x, y, n, 0, r, o->na_rm);
  else RSIMD_KERNEL(hamming_c128_)(x, y, n, flags, r, o->na_rm);
}

void RSIMD_KERNEL(na_c128)(const Rcomplex *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                           rsimd_reduce_result *r);
void RSIMD_KERNEL(na_c128)(const Rcomplex *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                           rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  ptrdiff_t i = 0, j;
  rsimd_vf64 re, im;
  rsimd_mf64 m;
  /* An element is missing when one of its doubles is NaN. */
  if (mode == RSIMD_NAMODE_ANY) {
    RSIMD_KERNEL(na_f64)((const double *) x, 2 * n, mode, 0, NULL, r);
    return;
  }
  for (; i + W <= n; i += W) {
    rsimd_c128_load_(x, i, &re, &im);
    m = rsimd_mf64_or(rsimd_vf64_is_nan(re), rsimd_vf64_is_nan(im));
    if (mode == RSIMD_NAMODE_COUNT) {
      r->i64 += rsimd_mf64_count(m);
    } else if (rsimd_mf64_any(m)) {
      for (j = i; j < i + W; j++) {
        if (isnan(x[j].r) || isnan(x[j].i)) RSIMD_KERNEL(put_index_)(mode, out, r, off + j);
      }
    }
  }
  for (; i < n; i++) {
    if (!(isnan(x[i].r) || isnan(x[i].i))) continue;
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_KERNEL(put_index_)(mode, out, r, off + i);
  }
}

/* ---- Multiply, divide, product, comparison: vector tiers ------------------- */

/* Interleaves the parts re and im of W complex numbers into out[i]. */
RSIMD_ALWAYS_INLINE void rsimd_c128_store_(Rcomplex *out, ptrdiff_t i, rsimd_vf64 re,
                                           rsimd_vf64 im) {
  double *d = (double *) out + 2 * i;
  rsimd_vf64_storeu(d, rsimd_vf64_zip_lo(re, im));
  rsimd_vf64_storeu(d + RSIMD_LANES_64, rsimd_vf64_zip_hi(re, im));
}

/* rsimd_cmul_re_/rsimd_cmul_im_ on vectors. */
RSIMD_ALWAYS_INLINE void rsimd_cmul_v_(const int vr, const int vi, rsimd_vf64 a, rsimd_vf64 b,
                                       rsimd_vf64 c, rsimd_vf64 d, rsimd_vf64 *re,
                                       rsimd_vf64 *im) {
  switch (vr) {
  case RSIMD_CMUL_FMA1:
    *re = rsimd_vf64_fma(a, c, rsimd_vf64_neg(rsimd_vf64_mul(b, d)));
    break;
  case RSIMD_CMUL_FMA2:
    *re = rsimd_vf64_fma(rsimd_vf64_neg(b), d, rsimd_vf64_mul(a, c));
    break;
  default: *re = rsimd_vf64_sub(rsimd_vf64_mul(a, c), rsimd_vf64_mul(b, d)); break;
  }
  switch (vi) {
  case RSIMD_CMUL_FMA1: *im = rsimd_vf64_fma(b, c, rsimd_vf64_mul(a, d)); break;
  case RSIMD_CMUL_FMA2: *im = rsimd_vf64_fma(a, d, rsimd_vf64_mul(b, c)); break;
  default: *im = rsimd_vf64_add(rsimd_vf64_mul(a, d), rsimd_vf64_mul(b, c)); break;
  }
}

/* (a + bi) / (c + di) as libgcc's __divdc3 computes it, without branches:
   both of its branches on |c| < |d| are one formula with the operands
   swapped by blends (big is the divisor part of larger magnitude), the
   scaling by 0.5 or 2^52 is a blended exact factor, and the form for a
   subnormal or zero ratio is blended in. Lanes where libgcc recovers
   infinities (or an operand is NaN) give a NaN part; the caller redoes
   them. `fused` contracts a*b + c as libgcc's build did. */
RSIMD_ALWAYS_INLINE void rsimd_cdiv_v_(const int fused, rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c,
                                       rsimd_vf64 d, rsimd_vf64 *re, rsimd_vf64 *im) {
  const rsimd_vf64 rmin = rsimd_vf64_set1(DBL_MIN), rmax2 = rsimd_vf64_set1(DBL_MAX / 2 * DBL_EPSILON);
  const rsimd_vf64 aa = rsimd_vf64_abs(a), ab = rsimd_vf64_abs(b), ac = rsimd_vf64_abs(c),
                   ad = rsimd_vf64_abs(d);
  const rsimd_mf64 sw = rsimd_vf64_cmp_lt(ac, ad);
  rsimd_vf64 big = rsimd_vf64_blend(c, d, sw), small = rsimd_vf64_blend(d, c, sw);
  const rsimd_vf64 abig = rsimd_vf64_blend(ac, ad, sw);
  const rsimd_mf64 halve = rsimd_vf64_cmp_ge(abig, rsimd_vf64_set1(DBL_MAX / 2));
  const rsimd_mf64 tiny_ab = rsimd_mf64_and(
    rsimd_vf64_cmp_lt(abig, rmax2),
    rsimd_mf64_or(rsimd_mf64_and(rsimd_vf64_cmp_lt(aa, rmin), rsimd_vf64_cmp_lt(ab, rmax2)),
                  rsimd_mf64_and(rsimd_vf64_cmp_lt(ab, rmin), rsimd_vf64_cmp_lt(aa, rmax2))));
  const rsimd_mf64 up = rsimd_mf64_or(rsimd_vf64_cmp_lt(abig, rsimd_vf64_set1(DBL_EPSILON)), tiny_ab);
  const rsimd_vf64 f = rsimd_vf64_blend(
    rsimd_vf64_blend(rsimd_vf64_set1(1.0), rsimd_vf64_set1(1.0 / DBL_EPSILON), up),
    rsimd_vf64_set1(0.5), halve);
  rsimd_vf64 u, v, s, t, ratio, denom, x, y, xs, ys;
  rsimd_mf64 alt;
  a = rsimd_vf64_mul(a, f);
  b = rsimd_vf64_mul(b, f);
  big = rsimd_vf64_mul(big, f);
  small = rsimd_vf64_mul(small, f);
  /* x = (u ratio + v) / denom and y = (s ratio + t) / denom: u, v = a, b
     when |c| < |d| (else b, a), s, t = b, -a (else -a, b). */
  u = rsimd_vf64_blend(b, a, sw);
  v = rsimd_vf64_blend(a, b, sw);
  s = rsimd_vf64_blend(rsimd_vf64_neg(a), b, sw);
  t = rsimd_vf64_blend(b, rsimd_vf64_neg(a), sw);
  ratio = rsimd_vf64_div(small, big);
  alt = rsimd_mf64_not(rsimd_vf64_cmp_gt(rsimd_vf64_abs(ratio), rmin));
  if (fused) {
    denom = rsimd_vf64_fma(small, ratio, big);
    x = rsimd_vf64_fma(u, ratio, v);
    y = rsimd_vf64_fma(s, ratio, t);
  } else {
    denom = rsimd_vf64_add(rsimd_vf64_mul(small, ratio), big);
    x = rsimd_vf64_add(rsimd_vf64_mul(u, ratio), v);
    y = rsimd_vf64_add(rsimd_vf64_mul(s, ratio), t);
  }
  if (rsimd_mf64_any(alt)) {
    /* x = (small (u / big) + v) / denom, and y likewise. */
    if (fused) {
      xs = rsimd_vf64_fma(small, rsimd_vf64_div(u, big), v);
      ys = rsimd_vf64_fma(small, rsimd_vf64_div(s, big), t);
    } else {
      xs = rsimd_vf64_add(rsimd_vf64_mul(small, rsimd_vf64_div(u, big)), v);
      ys = rsimd_vf64_add(rsimd_vf64_mul(small, rsimd_vf64_div(s, big)), t);
    }
    x = rsimd_vf64_blend(x, xs, alt);
    y = rsimd_vf64_blend(y, ys, alt);
  }
  *re = rsimd_vf64_div(x, denom);
  *im = rsimd_vf64_div(y, denom);
}

/* x * y (div = 0, variants vr and vi) or x / y (div = 1, `fused`) of whole
   vectors; elements with a NaN part are redone with base R's operator. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(ew2_c128_)(const int div, const int vr, const int vi,
                                                 const int fused, const Rcomplex *x,
                                                 const Rcomplex *y, R_xlen_t n, int flags,
                                                 Rcomplex *out, const rsimd_c128_arith *a) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  const rsimd_vf64 bxr = rsimd_vf64_set1(x[0].r), bxi = rsimd_vf64_set1(x[0].i),
                   byr = rsimd_vf64_set1(y[0].r), byi = rsimd_vf64_set1(y[0].i);
  const rsimd_c128_fn f = div ? a->div1 : a->mul1;
  ptrdiff_t i = 0, j;
  for (; i + W <= n; i += W) {
    rsimd_vf64 xr = bxr, xi = bxi, yr = byr, yi = byi, re, im;
    if (!sx) rsimd_c128_load_(x, i, &xr, &xi);
    if (!sy) rsimd_c128_load_(y, i, &yr, &yi);
    if (div) rsimd_cdiv_v_(fused, xr, xi, yr, yi, &re, &im);
    else rsimd_cmul_v_(vr, vi, xr, xi, yr, yi, &re, &im);
    rsimd_c128_store_(out, i, re, im);
    if (rsimd_mf64_any(rsimd_mf64_or(rsimd_vf64_is_nan(re), rsimd_vf64_is_nan(im)))) {
      for (j = i; j < i + W; j++) {
        if (isnan(out[j].r) || isnan(out[j].i)) f(x + (sx ? 0 : j), y + (sy ? 0 : j), out + j);
      }
    }
  }
  rsimd_ew2_c128_from(div ? RSIMD_EW_DIV : RSIMD_EW_MUL, x, y, i, n, flags, out, a);
}

void RSIMD_KERNEL(ew2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            Rcomplex *out, const rsimd_c128_arith *a);
void RSIMD_KERNEL(ew2_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            Rcomplex *out, const rsimd_c128_arith *a) {
  if (op == RSIMD_EW_DIV) {
    if (a->div == RSIMD_CDIV_FMA) RSIMD_KERNEL(ew2_c128_)(1, 0, 0, 1, x, y, n, flags, out, a);
    else if (a->div == RSIMD_CDIV_UNFUSED) RSIMD_KERNEL(ew2_c128_)(1, 0, 0, 0, x, y, n, flags, out, a);
    else rsimd_ew2_c128_from(op, x, y, 0, n, flags, out, a);
    return;
  }
  /* The variants of GCC builds (FMA1, FMA1 where it contracts, UNFUSED
     where it cannot) are specialised; the others pass their codes. */
  if (a->mul_re == RSIMD_CMUL_SCALAR || a->mul_im == RSIMD_CMUL_SCALAR) {
    rsimd_ew2_c128_from(op, x, y, 0, n, flags, out, a);
  } else if (a->mul_re == RSIMD_CMUL_FMA1 && a->mul_im == RSIMD_CMUL_FMA1) {
    RSIMD_KERNEL(ew2_c128_)(0, RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1, 0, x, y, n, flags, out, a);
  } else if (a->mul_re == RSIMD_CMUL_UNFUSED && a->mul_im == RSIMD_CMUL_UNFUSED) {
    RSIMD_KERNEL(ew2_c128_)(0, RSIMD_CMUL_UNFUSED, RSIMD_CMUL_UNFUSED, 0, x, y, n, flags, out, a);
  } else {
    RSIMD_KERNEL(ew2_c128_)(0, a->mul_re, a->mul_im, 0, x, y, n, flags, out, a);
  }
}

/* Fast mode: W lane products, each element multiplying into its lane (a
   skipped element as 1), combined pairwise at the end of the chunk; then
   the tail in order, and the chunk's product into s->p. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(prod_c128_)(const int vr, const int vi, const Rcomplex *x,
                                                  R_xlen_t n, rsimd_cprod_state *s,
                                                  const rsimd_c128_arith *a, const rsimd_opts *o) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int check = o->na_check || o->na_rm;
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), zero = rsimd_vf64_zero();
  rsimd_vf64 pr = one, pi = zero;
  double br[RSIMD_MAX_LANES_64], bi[RSIMD_MAX_LANES_64];
  rsimd_cprod_state t;
  ptrdiff_t i = 0, j, k;
  for (; i + W <= n; i += W) {
    rsimd_vf64 er, ei, nr, ni;
    rsimd_c128_load_(x, i, &er, &ei);
    if (check) {
      rsimd_mf64 m = rsimd_mf64_or(rsimd_vf64_is_nan(er), rsimd_vf64_is_nan(ei));
      if (rsimd_mf64_any(m)) {
        s->saw_nan = 1;
        if (rsimd_mf64_any(rsimd_mf64_or(rsimd_vf64_is_na(er), rsimd_vf64_is_na(ei)))) s->saw_na = 1;
        er = rsimd_vf64_blend(er, one, m);
        ei = rsimd_vf64_blend(ei, zero, m);
      }
    }
    rsimd_cmul_v_(vr, vi, pr, pi, er, ei, &nr, &ni);
    if (rsimd_mf64_any(rsimd_mf64_or(rsimd_vf64_is_nan(nr), rsimd_vf64_is_nan(ni)))) {
      double xr[RSIMD_MAX_LANES_64], xi[RSIMD_MAX_LANES_64], rr[RSIMD_MAX_LANES_64],
        ri[RSIMD_MAX_LANES_64];
      rsimd_vf64_storeu(br, pr);
      rsimd_vf64_storeu(bi, pi);
      rsimd_vf64_storeu(xr, er);
      rsimd_vf64_storeu(xi, ei);
      rsimd_vf64_storeu(rr, nr);
      rsimd_vf64_storeu(ri, ni);
      for (j = 0; j < W; j++) {
        if (isnan(rr[j]) || isnan(ri[j])) {
          Rcomplex p, e, r;
          p.r = br[j];
          p.i = bi[j];
          e.r = xr[j];
          e.i = xi[j];
          a->mul1(&p, &e, &r);
          rr[j] = r.r;
          ri[j] = r.i;
        }
      }
      nr = rsimd_vf64_loadu(rr);
      ni = rsimd_vf64_loadu(ri);
    }
    pr = nr;
    pi = ni;
  }
  rsimd_vf64_storeu(br, pr);
  rsimd_vf64_storeu(bi, pi);
  for (k = 1; k < W; k *= 2) {
    for (j = 0; j + k < W; j += 2 * k) {
      Rcomplex p, q;
      p.r = br[j];
      p.i = bi[j];
      q.r = br[j + k];
      q.i = bi[j + k];
      p = rsimd_cmul_1_(a, p, q);
      br[j] = p.r;
      bi[j] = p.i;
    }
  }
  t = *s;
  t.p.r = br[0];
  t.p.i = bi[0];
  rsimd_cprod_seq_(x, i, n, &t, a, o);
  s->saw_na = t.saw_na;
  s->saw_nan = t.saw_nan;
  s->p = rsimd_cmul_1_(a, s->p, t.p);
}

void RSIMD_KERNEL(prod_c128)(const Rcomplex *x, R_xlen_t n, rsimd_cprod_state *s,
                             const rsimd_c128_arith *a, const rsimd_opts *o);
void RSIMD_KERNEL(prod_c128)(const Rcomplex *x, R_xlen_t n, rsimd_cprod_state *s,
                             const rsimd_c128_arith *a, const rsimd_opts *o) {
  if (o->precision != RSIMD_PREC_FAST || a->mul_re == RSIMD_CMUL_SCALAR ||
      a->mul_im == RSIMD_CMUL_SCALAR) {
    rsimd_cprod_seq_(x, 0, n, s, a, o);
  } else if (a->mul_re == RSIMD_CMUL_FMA1 && a->mul_im == RSIMD_CMUL_FMA1) {
    RSIMD_KERNEL(prod_c128_)(RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1, x, n, s, a, o);
  } else {
    RSIMD_KERNEL(prod_c128_)(a->mul_re, a->mul_im, x, n, s, a, o);
  }
}

#define RSIMD_SKIP_scan_c128 1
#define RSIMD_SKIP_formula_c128 1

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cmp_c128_)(const int op, const Rcomplex *x,
                                                 const Rcomplex *y, R_xlen_t n, int flags,
                                                 int *out) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  const rsimd_vf64 bxr = rsimd_vf64_set1(x[0].r), bxi = rsimd_vf64_set1(x[0].i),
                   byr = rsimd_vf64_set1(y[0].r), byi = rsimd_vf64_set1(y[0].i);
  ptrdiff_t i = 0;
  for (; i + W <= n; i += W) {
    rsimd_vf64 xr = bxr, xi = bxi, yr = byr, yi = byi;
    rsimd_mf64 m, na;
    if (!sx) rsimd_c128_load_(x, i, &xr, &xi);
    if (!sy) rsimd_c128_load_(y, i, &yr, &yi);
    m = op == RSIMD_CMP_EQ
          ? rsimd_mf64_and(rsimd_vf64_cmp_eq(xr, yr), rsimd_vf64_cmp_eq(xi, yi))
          : rsimd_mf64_or(rsimd_vf64_cmp_ne(xr, yr), rsimd_vf64_cmp_ne(xi, yi));
    na = rsimd_mf64_or(rsimd_mf64_or(rsimd_vf64_is_nan(xr), rsimd_vf64_is_nan(xi)),
                       rsimd_mf64_or(rsimd_vf64_is_nan(yr), rsimd_vf64_is_nan(yi)));
    rsimd_vi64_storeu_i32(out + i, rsimd_lgl_vi64(m, na));
  }
  for (; i < n; i++) {
    Rcomplex u = x[sx ? 0 : i], v = y[sy ? 0 : i];
    if (isnan(u.r) || isnan(u.i) || isnan(v.r) || isnan(v.i)) out[i] = RSIMD_NA_I32;
    else out[i] = (u.r == v.r && u.i == v.i) == (op == RSIMD_CMP_EQ);
  }
}

void RSIMD_KERNEL(cmp_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            int *out);
void RSIMD_KERNEL(cmp_c128)(int op, const Rcomplex *x, const Rcomplex *y, R_xlen_t n, int flags,
                            int *out) {
  if (op == RSIMD_CMP_EQ) RSIMD_KERNEL(cmp_c128_)(RSIMD_CMP_EQ, x, y, n, flags, out);
  else RSIMD_KERNEL(cmp_c128_)(RSIMD_CMP_NE, x, y, n, flags, out);
}

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does complex */
#define RSIMD_SKIP_ew2_c128 1
#define RSIMD_SKIP_prod_c128 1
#define RSIMD_SKIP_scan_c128 1
#define RSIMD_SKIP_cmp_c128 1
#define RSIMD_SKIP_formula_c128 1
#define RSIMD_SKIP_conj_c128 1
#define RSIMD_SKIP_part_c128 1
#define RSIMD_SKIP_sum_c128 1
#define RSIMD_SKIP_pred_c128 1
#define RSIMD_SKIP_na_c128 1
#define RSIMD_SKIP_hamming_c128 1
#endif

#if RSIMD_TIER_IS(none) || !defined(RSIMD_NO_F64_SIMD)

/* The real parts fold into r[0] and the imaginary parts into r[1] with
   the tier's sum_f64, a block at a time. With na.rm the split has already
   zeroed every element with a missing part, so the parts are folded
   without missing-value checks and the removed elements are taken off both
   counts. */
void RSIMD_KERNEL(sum_c128)(const Rcomplex *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o);
void RSIMD_KERNEL(sum_c128)(const Rcomplex *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o) {
  double re[RSIMD_CPLX_BLOCK], im[RSIMD_CPLX_BLOCK];
  rsimd_opts po = *o;
  R_xlen_t i;
  if (o->na_rm) po.na_rm = po.na_check = 0;
  for (i = 0; i < n; i += RSIMD_CPLX_BLOCK) {
    ptrdiff_t len = n - i < RSIMD_CPLX_BLOCK ? (ptrdiff_t) (n - i) : RSIMD_CPLX_BLOCK;
    ptrdiff_t removed = RSIMD_KERNEL(split_c128_)(x + i, len, re, im, o->na_rm);
    RSIMD_KERNEL(sum_f64)(re, len, &r[0], &po);
    RSIMD_KERNEL(sum_f64)(im, len, &r[1], &po);
    r[0].count -= removed;
    r[1].count -= removed;
  }
}

#endif
