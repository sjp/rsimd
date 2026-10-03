/* Reduction kernels. Each folds one chunk of the input into the running
 * rsimd_reduce_result with the rules of na.h: precision mode, NA and NaN
 * tracking, na.rm. The none tier uses plain C loops and the scalar folds and
 * no SIMD code at all, so it is the reference the vector tiers are tested
 * against; every other tier uses the vector layer at its own width.
 *
 * Min and max: accumulators start at +-Inf and take an element only when
 * it is strictly smaller (larger), so NaN never enters them and the first
 * of equal values is kept. Equal values are identical except 0 and -0;
 * when a vector tier's extremum is zero it looks up the first zero of the
 * chunk, so the sign is the first one seen, as in base R. Chunks merge
 * with the same strict comparison, so the earlier chunk wins ties.
 */

/* Writes the 1-based index off + i + 1 to out[r->i64++] as an int or a
   double, for the which modes of the missing-value kernels. */
static inline void RSIMD_KERNEL(put_index_)(int mode, void *out, rsimd_reduce_result *r,
                                            R_xlen_t idx) {
  if (mode == RSIMD_NAMODE_WHICH_I32) ((int *) out)[r->i64] = (int) (idx + 1);
  else ((double *) out)[r->i64] = (double) (idx + 1);
  r->i64++;
}
#define RSIMD_PUT_INDEX RSIMD_KERNEL(put_index_)

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  rsimd_fold_f64(x, NULL, n, RSIMD_TERM_X, r, o);
}

void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_fold_sum_i32((const int32_t *) x, n, r, o);
}

void RSIMD_KERNEL(sum_dev_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                               const rsimd_opts *o);
void RSIMD_KERNEL(sum_dev_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                               const rsimd_opts *o) {
  rsimd_fold_f64_c(x, NULL, n, RSIMD_TERM_DEV, c, r, o);
}

void RSIMD_KERNEL(prod_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o);
void RSIMD_KERNEL(prod_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  double p = 1.0;
  R_xlen_t i, removed = 0;
  for (i = 0; i < n; i++) {
    if (check && isnan(x[i])) {
      r->saw_nan = 1;
      if (rsimd_is_na_f64(x[i])) r->saw_na = 1;
      if (o->na_rm) {
        removed++;
        continue;
      }
    }
    p *= x[i];
  }
  r->f64 *= p;
  r->count += n - removed;
}

void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  double p = 1.0;
  R_xlen_t i, removed = 0;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I32) {
      r->saw_na = 1;
      if (o->na_rm) {
        removed++;
        continue;
      }
    }
    p *= (double) x[i];
  }
  r->f64 *= p;
  r->count += n - removed;
}

void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  double lo = r->f64, hi = r->f64_hi;
  R_xlen_t i, missing = 0;
  for (i = 0; i < n; i++) {
    double v = x[i];
    if (isnan(v)) {
      if (check) {
        missing++;
        r->saw_nan = 1;
        if (rsimd_is_na_f64(v)) r->saw_na = 1;
      }
      continue;
    }
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  r->f64 = lo;
  r->f64_hi = hi;
  r->count += n - missing;
}

void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  int64_t lo = r->i64, hi = r->i64_hi;
  R_xlen_t i, missing = 0;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I32) {
      missing++;
      r->saw_na = 1;
      continue;
    }
    if (x[i] < lo) lo = x[i];
    if (x[i] > hi) hi = x[i];
  }
  r->i64 = lo;
  r->i64_hi = hi;
  r->count += n - missing;
}

R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, double v);
R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, double v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (x[i] == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int v);
R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (x[i] == v) return i;
  }
  return -1;
}

void RSIMD_KERNEL(which_u8)(const Rbyte *x, R_xlen_t n, R_xlen_t off, int max,
                            rsimd_reduce_result *r);
void RSIMD_KERNEL(which_u8)(const Rbyte *x, R_xlen_t n, R_xlen_t off, int max,
                            rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (max ? x[i] > r->i64 : x[i] < r->i64) {
      r->i64 = x[i];
      r->idx = off + i;
    }
  }
  r->count += n;
}

void RSIMD_KERNEL(anyall_f64)(const double *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_f64)(const double *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (check && isnan(x[i])) {
      r->saw_na = 1;
    } else if (x[i] == 0) {
      r->any_false = 1;
      if (stop == RSIMD_STOP_FALSE) return;
    } else {
      r->any_true = 1;
      if (stop == RSIMD_STOP_TRUE) return;
    }
  }
}

void RSIMD_KERNEL(anyall_i32)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_i32)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  rsimd_fold_lgl((const int32_t *) x, n, stop, r, o);
}

void RSIMD_KERNEL(anyall_u8)(const Rbyte *x, R_xlen_t n, int stop, rsimd_reduce_result *r);
void RSIMD_KERNEL(anyall_u8)(const Rbyte *x, R_xlen_t n, int stop, rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (x[i] == 0) {
      r->any_false = 1;
      if (stop == RSIMD_STOP_FALSE) return;
    } else {
      r->any_true = 1;
      if (stop == RSIMD_STOP_TRUE) return;
    }
  }
}

void RSIMD_KERNEL(na_f64)(const double *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_f64)(const double *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (!isnan(x[i])) continue;
    if (mode == RSIMD_NAMODE_ANY) {
      r->saw_na = 1;
      return;
    }
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_PUT_INDEX(mode, out, r, off + i);
  }
}

void RSIMD_KERNEL(na_i32)(const int *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_i32)(const int *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (x[i] != RSIMD_NA_I32) continue;
    if (mode == RSIMD_NAMODE_ANY) {
      r->saw_na = 1;
      return;
    }
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_PUT_INDEX(mode, out, r, off + i);
  }
}

/* Complex values have no vector kernels: every tier uses this one. */
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
    else RSIMD_PUT_INDEX(mode, out, r, off + i);
  }
}

#else /* vector tiers */

/* Raw which_min/which_max and complex missing values have no vector
   kernels; the none tier's are used. */
#define RSIMD_SKIP_which_u8 1
#define RSIMD_SKIP_na_c128 1

void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_vfold_sum_i32((const int32_t *) x, n, r, o);
}

/* ---- Integer and logical ---- */

/* The number of NA elements among the n. */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(count_na_i32_)(const int *x, R_xlen_t n) {
  const ptrdiff_t W = RSIMD_LANES_32;
  R_xlen_t i = 0, k = 0;
  for (; i + W <= n; i += W) k += rsimd_mi32_count(rsimd_vi32_is_na(rsimd_vi32_loadu(x + i)));
  for (; i < n; i++) k += x[i] == RSIMD_NA_I32;
  return k;
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(minmax_i32_)(const int *x, R_xlen_t n, const int check,
                                                  rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const rsimd_vi32 big = rsimd_vi32_set1(INT32_MAX);
  rsimd_vi32 lo0 = big, lo1 = big, hi0 = rsimd_vi32_set1(INT32_MIN), hi1 = hi0;
  rsimd_mi32 mna = rsimd_mi32_none();
  R_xlen_t i = 0, missing = 0;
  /* NA (INT32_MIN) lanes cannot raise the maximum; for the minimum they
     are replaced by INT32_MAX. */
#define RSIMD_MINMAX_I32_(lo, hi, v)                                                     \
  do {                                                                                   \
    rsimd_vi32 m_ = (v);                                                                 \
    if (check) {                                                                         \
      rsimd_mi32 na_ = rsimd_vi32_is_na(v);                                              \
      mna = rsimd_mi32_or(mna, na_);                                                     \
      m_ = rsimd_vi32_blend(m_, big, na_);                                               \
    }                                                                                    \
    (lo) = rsimd_vi32_min((lo), m_);                                                     \
    (hi) = rsimd_vi32_max((hi), (v));                                                    \
  } while (0)
  for (; i + 2 * W <= n; i += 2 * W) {
    rsimd_vi32 v0 = rsimd_vi32_loadu(x + i), v1 = rsimd_vi32_loadu(x + i + W);
    RSIMD_MINMAX_I32_(lo0, hi0, v0);
    RSIMD_MINMAX_I32_(lo1, hi1, v1);
  }
  for (; i < n; i += W) {
    /* The tail repeats element i in the inactive lanes. */
    rsimd_vi32 v = rsimd_vi32_loadu_p(rsimd_p32_while(i, n), x + i, x[i]);
    RSIMD_MINMAX_I32_(lo0, hi0, v);
  }
#undef RSIMD_MINMAX_I32_
  if (check && rsimd_mi32_any(mna)) {
    r->saw_na = 1;
    missing = RSIMD_KERNEL(count_na_i32_)(x, n);
  }
  if (missing < n) {
    int32_t lo = rsimd_vi32_reduce_min(rsimd_vi32_min(lo0, lo1));
    int32_t hi = rsimd_vi32_reduce_max(rsimd_vi32_max(hi0, hi1));
    if (lo < r->i64) r->i64 = lo;
    if (hi > r->i64_hi) r->i64_hi = hi;
  }
  r->count += n - missing;
}

void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  if (o->na_check || o->na_rm) RSIMD_KERNEL(minmax_i32_)(x, n, 1, r);
  else RSIMD_KERNEL(minmax_i32_)(x, n, 0, r);
}

R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int v);
R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int v) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const rsimd_vi32 vv = rsimd_vi32_set1(v);
  R_xlen_t i = 0;
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_mi32 m = rsimd_mi32_or(
        rsimd_mi32_or(rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i), vv),
                      rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + W), vv)),
        rsimd_mi32_or(rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + 2 * W), vv),
                      rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + 3 * W), vv)));
    if (rsimd_mi32_any(m)) break;
  }
  for (; i < n; i++) {
    if (x[i] == v) return i;
  }
  return -1;
}

void RSIMD_KERNEL(anyall_i32)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_i32)(const int *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  rsimd_vfold_lgl((const int32_t *) x, n, stop, r, o);
}

/* Raw bytes are read as 32-bit words: a word is non-zero when any of its
   bytes is, and (w - 0x01010101) & ~w & 0x80808080 is non-zero exactly
   when one of its bytes is zero. */
void RSIMD_KERNEL(anyall_u8)(const Rbyte *x, R_xlen_t n, int stop, rsimd_reduce_result *r);
void RSIMD_KERNEL(anyall_u8)(const Rbyte *x, R_xlen_t n, int stop, rsimd_reduce_result *r) {
  const ptrdiff_t B = 4 * RSIMD_LANES_32; /* bytes per vector */
  const rsimd_vi32 zero = rsimd_vi32_zero(), ones = rsimd_vi32_set1(0x01010101),
                   highs = rsimd_vi32_set1((int32_t) -2139062144); /* 0x80808080 */
  R_xlen_t i = 0;
  while (i + B <= n) {
    R_xlen_t end = i + 4 * B <= n ? i + 4 * B : n - (n - i) % B;
    rsimd_mi32 mt = rsimd_mi32_none(), mf = mt;
    for (; i < end; i += B) {
      rsimd_vi32 w = rsimd_vi32_loadu((const int32_t *) (const void *) (x + i));
      rsimd_vi32 z = rsimd_vi32_and(rsimd_vi32_andnot(w, rsimd_vi32_sub(w, ones)), highs);
      mt = rsimd_mi32_or(mt, rsimd_mi32_not(rsimd_vi32_cmp_eq(w, zero)));
      mf = rsimd_mi32_or(mf, rsimd_mi32_not(rsimd_vi32_cmp_eq(z, zero)));
    }
    if (rsimd_mi32_any(mt)) r->any_true = 1;
    if (rsimd_mi32_any(mf)) r->any_false = 1;
    if ((stop == RSIMD_STOP_TRUE && r->any_true) || (stop == RSIMD_STOP_FALSE && r->any_false)) {
      return;
    }
  }
  for (; i < n; i++) {
    if (x[i] == 0) {
      r->any_false = 1;
      if (stop == RSIMD_STOP_FALSE) return;
    } else {
      r->any_true = 1;
      if (stop == RSIMD_STOP_TRUE) return;
    }
  }
}

void RSIMD_KERNEL(na_i32)(const int *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_i32)(const int *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_32;
  R_xlen_t i = 0, j;
  if (mode == RSIMD_NAMODE_ANY) {
    for (; i + 4 * W <= n; i += 4 * W) {
      rsimd_mi32 m = rsimd_mi32_or(rsimd_mi32_or(rsimd_vi32_is_na(rsimd_vi32_loadu(x + i)),
                                                 rsimd_vi32_is_na(rsimd_vi32_loadu(x + i + W))),
                                   rsimd_mi32_or(rsimd_vi32_is_na(rsimd_vi32_loadu(x + i + 2 * W)),
                                                 rsimd_vi32_is_na(rsimd_vi32_loadu(x + i + 3 * W))));
      if (rsimd_mi32_any(m)) break;
    }
    for (; i < n; i++) {
      if (x[i] == RSIMD_NA_I32) {
        r->saw_na = 1;
        return;
      }
    }
  } else if (mode == RSIMD_NAMODE_COUNT) {
    r->i64 += RSIMD_KERNEL(count_na_i32_)(x, n);
  } else {
    for (; i + W <= n; i += W) {
      if (!rsimd_mi32_any(rsimd_vi32_is_na(rsimd_vi32_loadu(x + i)))) continue;
      for (j = i; j < i + W; j++) {
        if (x[j] == RSIMD_NA_I32) RSIMD_PUT_INDEX(mode, out, r, off + j);
      }
    }
    for (; i < n; i++) {
      if (x[i] == RSIMD_NA_I32) RSIMD_PUT_INDEX(mode, out, r, off + i);
    }
  }
}

/* ---- Double ---- */

#ifndef RSIMD_NO_F64_SIMD

void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  rsimd_vfold_f64(x, NULL, n, RSIMD_TERM_X, r, o);
}

void RSIMD_KERNEL(sum_dev_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                               const rsimd_opts *o);
void RSIMD_KERNEL(sum_dev_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                               const rsimd_opts *o) {
  rsimd_vfold_f64_c(x, NULL, n, RSIMD_TERM_DEV, c, r, o);
}

/* The lanes of v multiplied pairwise, as rsimd_vf64_sum_tree adds them. */
RSIMD_INLINE double RSIMD_KERNEL(prod_tree_)(rsimd_vf64 v) {
  double b[RSIMD_MAX_LANES_64];
  ptrdiff_t w = RSIMD_LANES_64, j;
  rsimd_vf64_storeu(b, v);
  while (w > 1) {
    ptrdiff_t h = w / 2;
    for (j = 0; j < h; j++) b[j] = b[2 * j] * b[2 * j + 1];
    if (w & 1) b[h++] = b[w - 1];
    w = h;
  }
  return b[0];
}

/* Four vector accumulators of products, as the fast sum. With check, NaN
   lanes are recorded (and set to 1 under narm). */
#define RSIMD_PROD_STEP_(acc, v, m)                                                      \
  do {                                                                                   \
    rsimd_vf64 v_ = (v);                                                                 \
    if (check) {                                                                         \
      mnan = rsimd_mf64_or(mnan, (m));                                                   \
      if (narm) {                                                                        \
        v_ = rsimd_vf64_blend(v_, one, (m));                                             \
        removed += rsimd_mf64_count(m);                                                  \
      }                                                                                  \
    }                                                                                    \
    (acc) = rsimd_vf64_mul((acc), v_);                                                   \
  } while (0)

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(prod_f64_)(const double *x, R_xlen_t n, const int check,
                                                const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 one = rsimd_vf64_set1(1.0);
  rsimd_vf64 a0 = one, a1 = one, a2 = one, a3 = one;
  rsimd_mf64 mnan = rsimd_mf64_none();
  R_xlen_t i = 0, removed = 0;
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_vf64 v0 = rsimd_vf64_loadu(x + i), v1 = rsimd_vf64_loadu(x + i + W);
    rsimd_vf64 v2 = rsimd_vf64_loadu(x + i + 2 * W), v3 = rsimd_vf64_loadu(x + i + 3 * W);
    RSIMD_PROD_STEP_(a0, v0, rsimd_vf64_is_nan(v0));
    RSIMD_PROD_STEP_(a1, v1, rsimd_vf64_is_nan(v1));
    RSIMD_PROD_STEP_(a2, v2, rsimd_vf64_is_nan(v2));
    RSIMD_PROD_STEP_(a3, v3, rsimd_vf64_is_nan(v3));
  }
  for (; i < n; i += W) {
    rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, n), x + i, 1.0);
    RSIMD_PROD_STEP_(a0, v, rsimd_vf64_is_nan(v));
  }
  if (check) rsimd_vfold_done_(mnan, removed, x, x, n, RSIMD_TERM_X, r);
  else r->count += n;
  r->f64 *= RSIMD_KERNEL(prod_tree_)(rsimd_vf64_mul(rsimd_vf64_mul(a0, a1),
                                                    rsimd_vf64_mul(a2, a3)));
}

void RSIMD_KERNEL(prod_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o);
void RSIMD_KERNEL(prod_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                            const rsimd_opts *o) {
  if (o->na_rm) RSIMD_KERNEL(prod_f64_)(x, n, 1, 1, r);
  else if (o->na_check) RSIMD_KERNEL(prod_f64_)(x, n, 1, 0, r);
  else RSIMD_KERNEL(prod_f64_)(x, n, 0, 0, r);
}

/* Integers are converted to double lanes (exactly) and multiplied there;
   NA converts to RSIMD_NA_I32_AS_F64. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(prod_i32_)(const int *x, R_xlen_t n, const int check,
                                                const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), na = rsimd_vf64_set1(RSIMD_NA_I32_AS_F64);
  rsimd_vf64 a0 = one, a1 = one;
  rsimd_mf64 mnan = rsimd_mf64_none(); /* the NA lanes */
  R_xlen_t i = 0, removed = 0;
  for (; i + 2 * W <= n; i += 2 * W) {
    rsimd_vf64 v0 = rsimd_vf64_loadu_i32(x + i), v1 = rsimd_vf64_loadu_i32(x + i + W);
    RSIMD_PROD_STEP_(a0, v0, rsimd_vf64_cmp_eq(v0, na));
    RSIMD_PROD_STEP_(a1, v1, rsimd_vf64_cmp_eq(v1, na));
  }
  for (; i < n; i += W) {
    rsimd_vf64 v = rsimd_vf64_loadu_i32_p(rsimd_p64_while(i, n), x + i, 1);
    RSIMD_PROD_STEP_(a0, v, rsimd_vf64_cmp_eq(v, na));
  }
  if (check && rsimd_mf64_any(mnan)) r->saw_na = 1;
  r->count += n - removed;
  r->f64 *= RSIMD_KERNEL(prod_tree_)(rsimd_vf64_mul(a0, a1));
}
#undef RSIMD_PROD_STEP_

void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  if (o->na_rm) RSIMD_KERNEL(prod_i32_)(x, n, 1, 1, r);
  else if (o->na_check) RSIMD_KERNEL(prod_i32_)(x, n, 1, 0, r);
  else RSIMD_KERNEL(prod_i32_)(x, n, 0, 0, r);
}

/* The index of the first element equal to v, or -1. */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(find_f64_)(const double *x, R_xlen_t n, double v) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 vv = rsimd_vf64_set1(v);
  R_xlen_t i = 0;
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_mf64 m = rsimd_mf64_or(
        rsimd_mf64_or(rsimd_vf64_cmp_eq(rsimd_vf64_loadu(x + i), vv),
                      rsimd_vf64_cmp_eq(rsimd_vf64_loadu(x + i + W), vv)),
        rsimd_mf64_or(rsimd_vf64_cmp_eq(rsimd_vf64_loadu(x + i + 2 * W), vv),
                      rsimd_vf64_cmp_eq(rsimd_vf64_loadu(x + i + 3 * W), vv)));
    if (rsimd_mf64_any(m)) break;
  }
  for (; i < n; i++) {
    if (x[i] == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, double v);
R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, double v) {
  return RSIMD_KERNEL(find_f64_)(x, n, v);
}

/* The number of NaN (NA included) elements among the n; sets *na if one
   of them is NA. */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(count_nan_f64_)(const double *x, R_xlen_t n, int *na) {
  const ptrdiff_t W = RSIMD_LANES_64;
  R_xlen_t i = 0, k = 0;
  rsimd_mf64 mna = rsimd_mf64_none();
  for (; i + W <= n; i += W) {
    rsimd_vf64 v = rsimd_vf64_loadu(x + i);
    k += rsimd_mf64_count(rsimd_vf64_is_nan(v));
    mna = rsimd_mf64_or(mna, rsimd_vf64_is_na(v));
  }
  if (rsimd_mf64_any(mna)) *na = 1;
  for (; i < n; i++) {
    if (isnan(x[i])) {
      k++;
      if (rsimd_is_na_f64(x[i])) *na = 1;
    }
  }
  return k;
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(minmax_f64_)(const double *x, R_xlen_t n, const int check,
                                                  rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  rsimd_vf64 lo0 = rsimd_vf64_set1(HUGE_VAL), lo1 = lo0;
  rsimd_vf64 hi0 = rsimd_vf64_set1(-HUGE_VAL), hi1 = hi0;
  rsimd_mf64 mnan = rsimd_mf64_none();
  R_xlen_t i = 0, missing = 0;
  double lo, hi;
  /* min(v, lo) is v < lo ? v : lo: strict, and lo when v is NaN. */
#define RSIMD_MINMAX_F64_(lo, hi, v)                                                     \
  do {                                                                                   \
    (lo) = rsimd_vf64_min((v), (lo));                                                    \
    (hi) = rsimd_vf64_max((v), (hi));                                                    \
    if (check) mnan = rsimd_mf64_or(mnan, rsimd_vf64_is_nan(v));                         \
  } while (0)
  for (; i + 2 * W <= n; i += 2 * W) {
    rsimd_vf64 v0 = rsimd_vf64_loadu(x + i), v1 = rsimd_vf64_loadu(x + i + W);
    RSIMD_MINMAX_F64_(lo0, hi0, v0);
    RSIMD_MINMAX_F64_(lo1, hi1, v1);
  }
  for (; i < n; i += W) {
    /* The tail repeats element i in the inactive lanes. */
    rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, n), x + i, x[i]);
    RSIMD_MINMAX_F64_(lo0, hi0, v);
  }
#undef RSIMD_MINMAX_F64_
  if (check && rsimd_mf64_any(mnan)) {
    int na = 0;
    r->saw_nan = 1;
    missing = RSIMD_KERNEL(count_nan_f64_)(x, n, &na);
    if (na) r->saw_na = 1;
  }
  r->count += n - missing;
  /* The accumulators hold no NaN, so the lane order of the reduction does
     not matter, except for the sign of zero. */
  lo = rsimd_vf64_reduce_min(rsimd_vf64_min(lo1, lo0));
  hi = rsimd_vf64_reduce_max(rsimd_vf64_max(hi1, hi0));
  if (lo == 0) lo = x[RSIMD_KERNEL(find_f64_)(x, n, 0.0)];
  if (hi == 0) hi = x[RSIMD_KERNEL(find_f64_)(x, n, 0.0)];
  if (lo < r->f64) r->f64 = lo;
  if (hi > r->f64_hi) r->f64_hi = hi;
}

void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  if (o->na_check || o->na_rm) RSIMD_KERNEL(minmax_f64_)(x, n, 1, r);
  else RSIMD_KERNEL(minmax_f64_)(x, n, 0, r);
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(anyall_f64_)(const double *x, R_xlen_t n, int stop,
                                                  const int check, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_mf64 mt = rsimd_mf64_none(), mf = mt, mna = mt;
  R_xlen_t i = 0;
  while (i < n) {
    R_xlen_t end = n - i > 4 * W ? i + 4 * W : n;
    for (; i < end; i += W) {
      /* The tail repeats element i in the inactive lanes. */
      rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, end), x + i, x[i]);
      rsimd_mf64 nan = check ? rsimd_vf64_is_nan(v) : rsimd_mf64_none();
      mna = rsimd_mf64_or(mna, nan);
      mt = rsimd_mf64_or(mt, rsimd_mf64_andnot(nan, rsimd_vf64_cmp_ne(v, zero)));
      mf = rsimd_mf64_or(mf, rsimd_vf64_cmp_eq(v, zero));
    }
    if ((stop == RSIMD_STOP_TRUE && rsimd_mf64_any(mt)) ||
        (stop == RSIMD_STOP_FALSE && rsimd_mf64_any(mf))) {
      break;
    }
  }
  if (rsimd_mf64_any(mt)) r->any_true = 1;
  if (rsimd_mf64_any(mf)) r->any_false = 1;
  if (rsimd_mf64_any(mna)) r->saw_na = 1;
}

void RSIMD_KERNEL(anyall_f64)(const double *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_f64)(const double *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  if (o->na_check || o->na_rm) RSIMD_KERNEL(anyall_f64_)(x, n, stop, 1, r);
  else RSIMD_KERNEL(anyall_f64_)(x, n, stop, 0, r);
}

void RSIMD_KERNEL(na_f64)(const double *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_f64)(const double *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  R_xlen_t i = 0, j;
  if (mode == RSIMD_NAMODE_ANY) {
    for (; i + 4 * W <= n; i += 4 * W) {
      rsimd_mf64 m = rsimd_mf64_or(
          rsimd_mf64_or(rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i)),
                        rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i + W))),
          rsimd_mf64_or(rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i + 2 * W)),
                        rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i + 3 * W))));
      if (rsimd_mf64_any(m)) break;
    }
    for (; i < n; i++) {
      if (isnan(x[i])) {
        r->saw_na = 1;
        return;
      }
    }
  } else if (mode == RSIMD_NAMODE_COUNT) {
    int na = 0;
    r->i64 += RSIMD_KERNEL(count_nan_f64_)(x, n, &na);
  } else {
    for (; i + W <= n; i += W) {
      if (!rsimd_mf64_any(rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i)))) continue;
      for (j = i; j < i + W; j++) {
        if (isnan(x[j])) RSIMD_PUT_INDEX(mode, out, r, off + j);
      }
    }
    for (; i < n; i++) {
      if (isnan(x[i])) RSIMD_PUT_INDEX(mode, out, r, off + i);
    }
  }
}

#else /* RSIMD_NO_F64_SIMD: double kernels come from the none tier */
#define RSIMD_SKIP_sum_f64 1
#define RSIMD_SKIP_sum_dev_f64 1
#define RSIMD_SKIP_prod_f64 1
#define RSIMD_SKIP_prod_i32 1
#define RSIMD_SKIP_minmax_f64 1
#define RSIMD_SKIP_find_f64 1
#define RSIMD_SKIP_anyall_f64 1
#define RSIMD_SKIP_na_f64 1
#endif /* RSIMD_NO_F64_SIMD */

#endif /* vector tiers */

#undef RSIMD_PUT_INDEX
