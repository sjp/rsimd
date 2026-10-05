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

#if RSIMD_TIER_IS(none)

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

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does complex */
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
