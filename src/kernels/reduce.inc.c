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

/* Writes the 1-based index idx + 1 to out[r->i64++] as an int or a
   double, for the which modes of the missing-value kernels (here and in
   complex.inc.c and int64.inc.c). */
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
  rsimd_fold_sum_i32((const int32_t *) x, n, 0, r, o);
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
      if (rsimd_is_na_f64(x[i])) {
        r->saw_na = 1;
        /* Without na.rm the product is NA whatever follows. */
        if (!o->na_rm) {
          r->count += n;
          return;
        }
      }
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
      if (!o->na_rm) {
        r->count += n;
        return;
      }
      removed++;
      continue;
    }
    p *= (double) x[i];
  }
  r->f64 *= p;
  r->count += n - removed;
}

void RSIMD_KERNEL(prod2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(prod2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             rsimd_reduce_result *r, const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  int xi = (flags & RSIMD_EW_I32(0)) != 0, yi = (flags & RSIMD_EW_I32(1)) != 0;
  R_xlen_t xs = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, ys = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  double p = 1.0;
  R_xlen_t i, removed = 0;
  for (i = 0; i < n; i++) {
    double a = rsimd_elt_f64(x, xi, i * xs), b = rsimd_elt_f64(y, yi, i * ys);
    double v = op == RSIMD_EW_SUB ? a - b : a + b;
    if (check && isnan(v)) {
      r->saw_nan = 1;
      if (rsimd_is_na_f64(a) || rsimd_is_na_f64(b)) {
        r->saw_na = 1;
        if (!o->na_rm) {
          r->count += n;
          return;
        }
      }
      if (o->na_rm) {
        removed++;
        continue;
      }
    }
    p *= v;
  }
  r->f64 *= p;
  r->count += n - removed;
}

void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm, na = 0;
  double lo = r->f64, hi = r->f64_hi;
  R_xlen_t i = 0, end, missing = 0;
  /* Blocks of RSIMD_FOLD_BLOCK, after each of which an NA without na.rm
     ends the scan, as the result is then NA. */
  for (end = 0; end < n;) {
    end = n - end > RSIMD_FOLD_BLOCK ? end + RSIMD_FOLD_BLOCK : n;
    for (; i < end; i++) {
      double v = absval ? fabs(x[i]) : x[i];
      if (isnan(v)) {
        if (check) {
          missing++;
          r->saw_nan = 1;
          if (rsimd_is_na_f64(v)) na = 1;
        }
        continue;
      }
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    if (na && !o->na_rm) break;
  }
  if (na) r->saw_na = 1;
  if (na && !o->na_rm) {
    r->count += n;
    return;
  }
  r->f64 = lo;
  r->f64_hi = hi;
  r->count += n - missing;
}

void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm, na = 0;
  int64_t lo = r->i64, hi = r->i64_hi;
  R_xlen_t i = 0, end, missing = 0;
  /* Blocks of RSIMD_FOLD_BLOCK, after each of which an NA without na.rm
     ends the scan, as the result is then NA. */
  for (end = 0; end < n;) {
    end = n - end > RSIMD_FOLD_BLOCK ? end + RSIMD_FOLD_BLOCK : n;
    for (; i < end; i++) {
      /* |NA| wraps to NA. */
      int32_t v = absval ? rsimd_abs_i32(x[i]) : x[i];
      if (check && v == RSIMD_NA_I32) {
        missing++;
        na = 1;
        continue;
      }
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    if (na && !o->na_rm) break;
  }
  if (na) r->saw_na = 1;
  if (na && !o->na_rm) {
    r->count += n;
    return;
  }
  r->i64 = lo;
  r->i64_hi = hi;
  r->count += n - missing;
}

R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, int absval, double v);
R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, int absval, double v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if ((absval ? fabs(x[i]) : x[i]) == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int absval, int v);
R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int absval, int v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if ((absval ? rsimd_abs_i32(x[i]) : x[i]) == v) return i;
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
  const int32_t v = (mode & RSIMD_NAMODE_TRUE) ? 1 : RSIMD_NA_I32;
  R_xlen_t i, j, k, buf[64];
  mode &= ~RSIMD_NAMODE_TRUE;
  if (v == 1 && mode != RSIMD_NAMODE_COUNT) {
    /* TRUE is common: gather each block's indices without a branch per
       element, as the vector tiers do. */
    for (i = 0; i < n; i += 64) {
      R_xlen_t len = n - i < 64 ? n - i : 64;
      for (j = 0, k = 0; j < len; j++) {
        buf[k] = off + i + j;
        k += x[i + j] == 1;
      }
      for (j = 0; j < k; j++) RSIMD_PUT_INDEX(mode, out, r, buf[j]);
    }
    return;
  }
  if (mode == RSIMD_NAMODE_COUNT) {
    for (i = 0, k = 0; i < n; i++) k += x[i] == v;
    r->i64 += k;
    return;
  }
  for (i = 0; i < n; i++) {
    if (x[i] != v) continue;
    if (mode == RSIMD_NAMODE_ANY) {
      r->saw_na = 1;
      return;
    }
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_PUT_INDEX(mode, out, r, off + i);
  }
}

#else /* vector tiers */

/* Raw which_min/which_max has no vector kernel; the none tier's is used.
   The complex missing-value kernel is in complex.inc.c. */
#define RSIMD_SKIP_which_u8 1

void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_vfold_sum_i32((const int32_t *) x, n, 0, r, o);
}

/* ---- Integer and logical ---- */

/* Loop iterations per block of count_eq_i32_(): a 32-bit lane count
   gains at most one per iteration, so it stays below INT32_MAX. The
   vector layer test sets a small value to cover the flushes. */
#ifndef RSIMD_COUNT_BLOCK
#define RSIMD_COUNT_BLOCK ((R_xlen_t) 1 << 24)
#endif

/* The number of elements equal to v (NA, or TRUE) among the n, counted
   in the lanes of four accumulators (a mask count per vector is a slow
   horizontal step) and added up once per block. */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(count_eq_i32_)(const int *x, R_xlen_t n, int32_t v) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const rsimd_vi32 zero = rsimd_vi32_zero(), vv = rsimd_vi32_set1(v);
  R_xlen_t i = 0, k = 0;
  while (i + 4 * W <= n) {
    R_xlen_t end = (n - i) / (4 * W) > RSIMD_COUNT_BLOCK ? i + RSIMD_COUNT_BLOCK * 4 * W : n;
    rsimd_vi32 c0 = zero, c1 = zero, c2 = zero, c3 = zero;
    for (; i + 4 * W <= end; i += 4 * W) {
      c0 = rsimd_vi32_inc(c0, rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i), vv));
      c1 = rsimd_vi32_inc(c1, rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + W), vv));
      c2 = rsimd_vi32_inc(c2, rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + 2 * W), vv));
      c3 = rsimd_vi32_inc(c3, rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i + 3 * W), vv));
    }
    k += rsimd_vi32_reduce_add(rsimd_vi32_add(rsimd_vi32_add(c0, c1), rsimd_vi32_add(c2, c3)));
  }
  for (; i + W <= n; i += W) k += rsimd_mi32_count(rsimd_vi32_cmp_eq(rsimd_vi32_loadu(x + i), vv));
  for (; i < n; i++) k += x[i] == v;
  return k;
}

/* Minimum and/or maximum (ext) of the n integers, or of their absolute
   values, into r, in four accumulator pairs. With stop (na_check without
   na.rm) a block of RSIMD_FOLD_BLOCK elements with an NA ends the scan,
   as the result is then NA. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(minmax_i32_)(const int *x, R_xlen_t n, const int check,
                                                  const int stop, const int absval,
                                                  const int ext, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const int want_lo = ext & RSIMD_EXT_MIN, want_hi = ext & RSIMD_EXT_MAX;
  const rsimd_vi32 big = rsimd_vi32_set1(INT32_MAX);
  rsimd_vi32 lo0 = big, lo1 = big, lo2 = big, lo3 = big;
  rsimd_vi32 hi0 = rsimd_vi32_set1(INT32_MIN), hi1 = hi0, hi2 = hi0, hi3 = hi0;
  rsimd_mi32 mna = rsimd_mi32_none();
  R_xlen_t i = 0, missing = 0;
  /* NA (INT32_MIN) lanes cannot raise the maximum; for the minimum they
     are replaced by INT32_MAX. The absolute value of NA wraps to NA. */
#define RSIMD_MINMAX_I32_(lo, hi, v0)                                                    \
  do {                                                                                   \
    rsimd_vi32 v_ = absval ? rsimd_vi32_abs_wrap(v0) : (v0), m_ = v_;                    \
    if (check) {                                                                         \
      rsimd_mi32 na_ = rsimd_vi32_is_na(v_);                                             \
      mna = rsimd_mi32_or(mna, na_);                                                     \
      if (want_lo) m_ = rsimd_vi32_blend(m_, big, na_);                                  \
    }                                                                                    \
    if (want_lo) (lo) = rsimd_vi32_min((lo), m_);                                        \
    if (want_hi) (hi) = rsimd_vi32_max((hi), v_);                                        \
  } while (0)
  while (i + 4 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (; i + 4 * W <= end; i += 4 * W) {
      rsimd_vi32 v0 = rsimd_vi32_loadu(x + i), v1 = rsimd_vi32_loadu(x + i + W);
      rsimd_vi32 v2 = rsimd_vi32_loadu(x + i + 2 * W), v3 = rsimd_vi32_loadu(x + i + 3 * W);
      RSIMD_MINMAX_I32_(lo0, hi0, v0);
      RSIMD_MINMAX_I32_(lo1, hi1, v1);
      RSIMD_MINMAX_I32_(lo2, hi2, v2);
      RSIMD_MINMAX_I32_(lo3, hi3, v3);
    }
    if (stop && rsimd_mi32_any(mna)) {
      r->saw_na = 1;
      r->count += n;
      return;
    }
  }
  for (; i < n; i += W) {
    /* The tail repeats element i in the inactive lanes. */
    rsimd_vi32 v = rsimd_vi32_loadu_p(rsimd_p32_while(i, n), x + i, x[i]);
    RSIMD_MINMAX_I32_(lo0, hi0, v);
  }
#undef RSIMD_MINMAX_I32_
  if (check && rsimd_mi32_any(mna)) {
    r->saw_na = 1;
    if (stop) {
      r->count += n;
      return;
    }
    missing = RSIMD_KERNEL(count_eq_i32_)(x, n, RSIMD_NA_I32);
  }
  if (missing < n) {
    if (want_lo) {
      int32_t lo = rsimd_vi32_reduce_min(rsimd_vi32_min(rsimd_vi32_min(lo0, lo1), rsimd_vi32_min(lo2, lo3)));
      if (lo < r->i64) r->i64 = lo;
    }
    if (want_hi) {
      int32_t hi = rsimd_vi32_reduce_max(rsimd_vi32_max(rsimd_vi32_max(hi0, hi1), rsimd_vi32_max(hi2, hi3)));
      if (hi > r->i64_hi) r->i64_hi = hi;
    }
  }
  r->count += n - missing;
}

#define RSIMD_MINMAX_I32_EXT_(c, s, a)                                                   \
  do {                                                                                   \
    if (o->extrema == RSIMD_EXT_MIN) RSIMD_KERNEL(minmax_i32_)(x, n, c, s, a, RSIMD_EXT_MIN, r); \
    else if (o->extrema == RSIMD_EXT_MAX) RSIMD_KERNEL(minmax_i32_)(x, n, c, s, a, RSIMD_EXT_MAX, r); \
    else RSIMD_KERNEL(minmax_i32_)(x, n, c, s, a, RSIMD_EXT_BOTH, r);                   \
  } while (0)

void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i32)(const int *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  if (!check) {
    if (absval) RSIMD_MINMAX_I32_EXT_(0, 0, 1);
    else RSIMD_MINMAX_I32_EXT_(0, 0, 0);
  } else if (o->na_rm) {
    if (absval) RSIMD_MINMAX_I32_EXT_(1, 0, 1);
    else RSIMD_MINMAX_I32_EXT_(1, 0, 0);
  } else {
    if (absval) RSIMD_MINMAX_I32_EXT_(1, 1, 1);
    else RSIMD_MINMAX_I32_EXT_(1, 1, 0);
  }
}
#undef RSIMD_MINMAX_I32_EXT_

RSIMD_ALWAYS_INLINE R_xlen_t RSIMD_KERNEL(find_i32_)(const int *x, R_xlen_t n, const int absval,
                                                    int v) {
  const ptrdiff_t W = RSIMD_LANES_32;
  const rsimd_vi32 vv = rsimd_vi32_set1(v);
  R_xlen_t i = 0;
#define RSIMD_FIND_I32_(j)                                                               \
  rsimd_vi32_cmp_eq(absval ? rsimd_vi32_abs_wrap(rsimd_vi32_loadu(x + i + (j) * W))       \
                           : rsimd_vi32_loadu(x + i + (j) * W),                          \
                    vv)
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_mi32 m = rsimd_mi32_or(rsimd_mi32_or(RSIMD_FIND_I32_(0), RSIMD_FIND_I32_(1)),
                                 rsimd_mi32_or(RSIMD_FIND_I32_(2), RSIMD_FIND_I32_(3)));
    if (rsimd_mi32_any(m)) break;
  }
#undef RSIMD_FIND_I32_
  for (; i < n; i++) {
    if ((absval ? rsimd_abs_i32(x[i]) : x[i]) == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int absval, int v);
R_xlen_t RSIMD_KERNEL(find_i32)(const int *x, R_xlen_t n, int absval, int v) {
  if (absval) return RSIMD_KERNEL(find_i32_)(x, n, 1, v);
  return RSIMD_KERNEL(find_i32_)(x, n, 0, v);
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
  if (mode & RSIMD_NAMODE_TRUE) {
    mode &= ~RSIMD_NAMODE_TRUE;
    if (mode == RSIMD_NAMODE_COUNT) {
      r->i64 += RSIMD_KERNEL(count_eq_i32_)(x, n, 1);
      return;
    }
    /* TRUE is common, so the indices are gathered without a branch per
       element: each block's candidates go to buf, the count advancing only
       past the TRUE ones, and are then written out. */
    {
      R_xlen_t buf[64], k;
      for (; i < n; i += 64) {
        R_xlen_t len = n - i < 64 ? n - i : 64;
        for (j = 0, k = 0; j < len; j++) {
          buf[k] = off + i + j;
          k += x[i + j] == 1;
        }
        for (j = 0; j < k; j++) RSIMD_PUT_INDEX(mode, out, r, buf[j]);
      }
    }
  } else if (mode == RSIMD_NAMODE_ANY) {
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
    r->i64 += RSIMD_KERNEL(count_eq_i32_)(x, n, RSIMD_NA_I32);
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
   lanes are recorded (and set to 1 and counted in the lanes of rcnt under
   narm). Without narm the kernels go in blocks of RSIMD_FOLD_BLOCK, as the
   folds do, and stop after a block with an NA, as the product is then
   NA. */
#define RSIMD_PROD_STEP_(acc, v, m)                                                      \
  do {                                                                                   \
    rsimd_vf64 v_ = (v);                                                                 \
    if (check) {                                                                         \
      mnan = rsimd_mf64_or(mnan, (m));                                                   \
      if (narm) {                                                                        \
        v_ = rsimd_vf64_blend(v_, one, (m));                                             \
        rcnt = rsimd_vi64_inc(rcnt, rsimd_mf64_to_mi64(m));                              \
      }                                                                                  \
    }                                                                                    \
    (acc) = rsimd_vf64_mul((acc), v_);                                                   \
  } while (0)

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(prod_f64_)(const double *x, R_xlen_t n, const int check,
                                                const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 one = rsimd_vf64_set1(1.0);
  const int stop = check && !narm;
  rsimd_vf64 a0 = one, a1 = one, a2 = one, a3 = one;
  rsimd_mf64 mnan = rsimd_mf64_none();
  rsimd_vi64 rcnt = rsimd_vi64_zero();
  R_xlen_t i = 0, from, tail;
  int nan = 0;
  while (i + 4 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (from = i; i + 4 * W <= end; i += 4 * W) {
      rsimd_vf64 v0 = rsimd_vf64_loadu(x + i), v1 = rsimd_vf64_loadu(x + i + W);
      rsimd_vf64 v2 = rsimd_vf64_loadu(x + i + 2 * W), v3 = rsimd_vf64_loadu(x + i + 3 * W);
      RSIMD_PROD_STEP_(a0, v0, rsimd_vf64_is_nan(v0));
      RSIMD_PROD_STEP_(a1, v1, rsimd_vf64_is_nan(v1));
      RSIMD_PROD_STEP_(a2, v2, rsimd_vf64_is_nan(v2));
      RSIMD_PROD_STEP_(a3, v3, rsimd_vf64_is_nan(v3));
    }
    if (stop && rsimd_mf64_any(mnan)) {
      if (rsimd_vfold_block_na_(x, 0, x, 0, 0, from, i, n, r)) return;
      mnan = rsimd_mf64_none();
      nan = 1;
    }
  }
  tail = stop ? i : 0; /* the NA scan of rsimd_vfold_done_() starts here */
  for (; i < n; i += W) {
    rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, n), x + i, 1.0);
    RSIMD_PROD_STEP_(a0, v, rsimd_vf64_is_nan(v));
  }
  if (check) {
    rsimd_vfold_done_(mnan, nan, tail, rsimd_vi64_reduce_add(rcnt), x, 0, x, 0, n, RSIMD_TERM_X,
                      r);
  } else {
    r->count += n;
  }
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
  const int stop = check && !narm;
  rsimd_vf64 a0 = one, a1 = one, a2 = one, a3 = one;
  rsimd_mf64 mnan = rsimd_mf64_none(); /* the NA lanes */
  rsimd_vi64 rcnt = rsimd_vi64_zero();
  R_xlen_t i = 0;
  while (i + 4 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (; i + 4 * W <= end; i += 4 * W) {
      rsimd_vf64 v0 = rsimd_vf64_loadu_i32(x + i), v1 = rsimd_vf64_loadu_i32(x + i + W);
      rsimd_vf64 v2 = rsimd_vf64_loadu_i32(x + i + 2 * W), v3 = rsimd_vf64_loadu_i32(x + i + 3 * W);
      RSIMD_PROD_STEP_(a0, v0, rsimd_vf64_cmp_eq(v0, na));
      RSIMD_PROD_STEP_(a1, v1, rsimd_vf64_cmp_eq(v1, na));
      RSIMD_PROD_STEP_(a2, v2, rsimd_vf64_cmp_eq(v2, na));
      RSIMD_PROD_STEP_(a3, v3, rsimd_vf64_cmp_eq(v3, na));
    }
    if (stop && rsimd_mf64_any(mnan)) {
      r->saw_na = 1;
      r->count += n;
      return;
    }
  }
  for (; i < n; i += W) {
    rsimd_vf64 v = rsimd_vf64_loadu_i32_p(rsimd_p64_while(i, n), x + i, 1);
    RSIMD_PROD_STEP_(a0, v, rsimd_vf64_cmp_eq(v, na));
  }
  if (check && rsimd_mf64_any(mnan)) r->saw_na = 1;
  r->count += n - rsimd_vi64_reduce_add(rcnt);
  r->f64 *= RSIMD_KERNEL(prod_tree_)(rsimd_vf64_mul(rsimd_vf64_mul(a0, a1),
                                                    rsimd_vf64_mul(a2, a3)));
}

void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(prod_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r, const rsimd_opts *o) {
  if (o->na_rm) RSIMD_KERNEL(prod_i32_)(x, n, 1, 1, r);
  else if (o->na_check) RSIMD_KERNEL(prod_i32_)(x, n, 1, 0, r);
  else RSIMD_KERNEL(prod_i32_)(x, n, 0, 0, r);
}

/* Operand k of prod2 at i as doubles: its broadcast value bc, or the
   doubles or int32 elements (NA becoming NA_real_) there. The operand
   flags are tested per vector, as the branches always go the same way. */
RSIMD_ALWAYS_INLINE rsimd_vf64 RSIMD_KERNEL(prod2_ld_)(const void *p, int flags, int k,
                                                      rsimd_vf64 bc, ptrdiff_t i) {
  if (flags & RSIMD_EW_SCALAR(k)) return bc;
  return rsimd_vf64_load_elt_(p, (flags & RSIMD_EW_I32(k)) != 0, i);
}

/* As rsimd_vfold_block_na_() for the pairs [from, to) of prod2, a
   broadcast operand being its one element. */
static RSIMD_NOINLINE int RSIMD_KERNEL(prod2_block_na_)(const void *x, const void *y, int flags,
                                                       R_xlen_t from, R_xlen_t to, R_xlen_t n,
                                                       rsimd_reduce_result *r) {
  const int xi = (flags & RSIMD_EW_I32(0)) != 0, yi = (flags & RSIMD_EW_I32(1)) != 0;
  const int xb = (flags & RSIMD_EW_SCALAR(0)) != 0, yb = (flags & RSIMD_EW_SCALAR(1)) != 0;
  if (rsimd_any_na_elt_(xb ? x : rsimd_elt_ptr_(x, xi, from), xi, xb ? 1 : to - from) ||
      rsimd_any_na_elt_(yb ? y : rsimd_elt_ptr_(y, yi, from), yi, yb ? 1 : to - from)) {
    r->saw_na = 1;
    r->saw_nan = 1;
    r->count += n;
    return 1;
  }
  return 0;
}

/* The product of x[i] + y[i] in four vector accumulators as prod_f64_,
   y being negated (exactly) for x - y; the last n % W pairs are
   multiplied in a scalar loop. A pair is missing when its sum is NaN;
   saw_na then comes from the operands, as an NA operand makes its pair's
   sum NaN. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(prod2_f64_)(int op, const void *x, const void *y,
                                                 R_xlen_t n, int flags, const int check,
                                                 const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int xi = (flags & RSIMD_EW_I32(0)) != 0, yi = (flags & RSIMD_EW_I32(1)) != 0;
  const R_xlen_t xs = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, ys = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), sign = rsimd_vf64_set1(op == RSIMD_EW_SUB ? -0.0 : 0.0);
  const rsimd_vf64 bx = rsimd_vf64_set1(rsimd_elt_f64(x, xi, 0));
  const rsimd_vf64 by = rsimd_vf64_set1(rsimd_elt_f64(y, yi, 0));
  const int stop = check && !narm;
  rsimd_vf64 a0 = one, a1 = one, a2 = one, a3 = one;
  rsimd_mf64 mnan = rsimd_mf64_none();
  rsimd_vi64 rcnt = rsimd_vi64_zero();
  R_xlen_t i = 0, from, removed = 0, tail;
  double p = 1.0;
  int tail_nan = 0, nan = 0;
#define RSIMD_PROD2_LD_(j)                                                               \
  rsimd_vf64_add(RSIMD_KERNEL(prod2_ld_)(x, flags, 0, bx, i + (j) * W),                  \
                 rsimd_vf64_xor(RSIMD_KERNEL(prod2_ld_)(y, flags, 1, by, i + (j) * W), sign))
  while (i + 4 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (from = i; i + 4 * W <= end; i += 4 * W) {
      rsimd_vf64 v0 = RSIMD_PROD2_LD_(0), v1 = RSIMD_PROD2_LD_(1);
      rsimd_vf64 v2 = RSIMD_PROD2_LD_(2), v3 = RSIMD_PROD2_LD_(3);
      RSIMD_PROD_STEP_(a0, v0, rsimd_vf64_is_nan(v0));
      RSIMD_PROD_STEP_(a1, v1, rsimd_vf64_is_nan(v1));
      RSIMD_PROD_STEP_(a2, v2, rsimd_vf64_is_nan(v2));
      RSIMD_PROD_STEP_(a3, v3, rsimd_vf64_is_nan(v3));
    }
    if (stop && rsimd_mf64_any(mnan)) {
      if (RSIMD_KERNEL(prod2_block_na_)(x, y, flags, from, i, n, r)) return;
      mnan = rsimd_mf64_none();
      nan = 1;
    }
  }
  tail = stop ? i : 0; /* the pairs before tail hold no NA operand */
  for (; i + W <= n; i += W) {
    rsimd_vf64 v = RSIMD_PROD2_LD_(0);
    RSIMD_PROD_STEP_(a0, v, rsimd_vf64_is_nan(v));
  }
#undef RSIMD_PROD2_LD_
  for (; i < n; i++) {
    double a = rsimd_elt_f64(x, xi, i * xs), b = rsimd_elt_f64(y, yi, i * ys);
    double v = op == RSIMD_EW_SUB ? a - b : a + b;
    if (check && isnan(v)) {
      tail_nan = 1;
      if (narm) {
        removed++;
        continue;
      }
    }
    p *= v;
  }
  if (check && (nan || tail_nan || rsimd_mf64_any(mnan))) r->saw_nan = 1;
  if (check && (tail_nan || rsimd_mf64_any(mnan)) &&
      (rsimd_any_na_elt_(xs ? rsimd_elt_ptr_(x, xi, tail) : x, xi, xs ? n - tail : 1) ||
       rsimd_any_na_elt_(ys ? rsimd_elt_ptr_(y, yi, tail) : y, yi, ys ? n - tail : 1))) {
    r->saw_na = 1;
  }
  r->count += n - removed - rsimd_vi64_reduce_add(rcnt);
  r->f64 *= RSIMD_KERNEL(prod_tree_)(rsimd_vf64_mul(rsimd_vf64_mul(a0, a1),
                                                    rsimd_vf64_mul(a2, a3))) * p;
}

void RSIMD_KERNEL(prod2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(prod2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             rsimd_reduce_result *r, const rsimd_opts *o) {
  if (o->na_rm) RSIMD_KERNEL(prod2_f64_)(op, x, y, n, flags, 1, 1, r);
  else if (o->na_check) RSIMD_KERNEL(prod2_f64_)(op, x, y, n, flags, 1, 0, r);
  else RSIMD_KERNEL(prod2_f64_)(op, x, y, n, flags, 0, 0, r);
}
#undef RSIMD_PROD_STEP_

/* The index of the first element equal to v (whose absolute value is,
   with absval), or -1. */
RSIMD_ALWAYS_INLINE R_xlen_t RSIMD_KERNEL(find_f64_)(const double *x, R_xlen_t n,
                                                    const int absval, double v) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vf64 vv = rsimd_vf64_set1(v);
  R_xlen_t i = 0;
#define RSIMD_FIND_F64_(j)                                                               \
  rsimd_vf64_cmp_eq(absval ? rsimd_vf64_abs(rsimd_vf64_loadu(x + i + (j) * W))           \
                           : rsimd_vf64_loadu(x + i + (j) * W),                          \
                    vv)
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_mf64 m = rsimd_mf64_or(rsimd_mf64_or(RSIMD_FIND_F64_(0), RSIMD_FIND_F64_(1)),
                                 rsimd_mf64_or(RSIMD_FIND_F64_(2), RSIMD_FIND_F64_(3)));
    if (rsimd_mf64_any(m)) break;
  }
#undef RSIMD_FIND_F64_
  for (; i < n; i++) {
    if ((absval ? fabs(x[i]) : x[i]) == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, int absval, double v);
R_xlen_t RSIMD_KERNEL(find_f64)(const double *x, R_xlen_t n, int absval, double v) {
  if (absval) return RSIMD_KERNEL(find_f64_)(x, n, 1, v);
  return RSIMD_KERNEL(find_f64_)(x, n, 0, v);
}

/* The number of NaN (NA included) elements among the n, counted in the
   64-bit lanes of four accumulators as count_na_i32_(). */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(count_nan_f64_)(const double *x, R_xlen_t n) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vi64 zero = rsimd_vi64_zero();
  rsimd_vi64 c0 = zero, c1 = zero, c2 = zero, c3 = zero;
  R_xlen_t i = 0, k;
#define RSIMD_NAN_MASK_(j) rsimd_mf64_to_mi64(rsimd_vf64_is_nan(rsimd_vf64_loadu(x + i + (j) * W)))
  for (; i + 4 * W <= n; i += 4 * W) {
    c0 = rsimd_vi64_inc(c0, RSIMD_NAN_MASK_(0));
    c1 = rsimd_vi64_inc(c1, RSIMD_NAN_MASK_(1));
    c2 = rsimd_vi64_inc(c2, RSIMD_NAN_MASK_(2));
    c3 = rsimd_vi64_inc(c3, RSIMD_NAN_MASK_(3));
  }
  for (; i + W <= n; i += W) c0 = rsimd_vi64_inc(c0, RSIMD_NAN_MASK_(0));
#undef RSIMD_NAN_MASK_
  k = rsimd_vi64_reduce_add(rsimd_vi64_add(rsimd_vi64_add(c0, c1), rsimd_vi64_add(c2, c3)));
  for (; i < n; i++) k += isnan(x[i]) != 0;
  return k;
}

/* The lane minimum and maximum of v and acc, which holds no NaN: acc when v
   is NaN (the layer's min and max take their second operand then). Not
   FMINNM/FMAXNM: R's NA is a signalling NaN, for which they give NaN. */
#define RSIMD_MIN_ACC_(v, acc) rsimd_vf64_min((v), (acc))
#define RSIMD_MAX_ACC_(v, acc) rsimd_vf64_max((v), (acc))

/* Minimum and/or maximum (ext) of the n doubles, or of their absolute
   values, into r, in four accumulator pairs. NaN lanes cannot win
   (RSIMD_MIN_ACC_); without stop they are counted afterwards when the NaN
   mask is set. The mask tests the sum of each four vectors, through which
   a NaN carries, so it can also be set by Inf + -Inf: whether there is a
   NaN is then settled by an exact scan. With stop (na_check without na.rm)
   a block of RSIMD_FOLD_BLOCK elements with an NA ends the scan, as the
   result is then NA, and the count does not matter. The sign of a zero extremum is
   that of the first zero in x, as in base R, found by a second scan. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(minmax_f64_)(const double *x, R_xlen_t n, const int check,
                                                  const int stop, const int absval,
                                                  const int ext, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int want_lo = ext & RSIMD_EXT_MIN, want_hi = ext & RSIMD_EXT_MAX;
  rsimd_vf64 lo0 = rsimd_vf64_set1(HUGE_VAL), lo1 = lo0, lo2 = lo0, lo3 = lo0;
  rsimd_vf64 hi0 = rsimd_vf64_set1(-HUGE_VAL), hi1 = hi0, hi2 = hi0, hi3 = hi0;
  rsimd_mf64 mnan = rsimd_mf64_none();
  R_xlen_t i = 0, missing = 0, from, tail;
  int nan = 0;
  double lo, hi;
#define RSIMD_MINMAX_F64_(lo, hi, v0)                                                    \
  do {                                                                                   \
    rsimd_vf64 v_ = absval ? rsimd_vf64_abs(v0) : (v0);                                  \
    if (want_lo) (lo) = RSIMD_MIN_ACC_(v_, (lo));                                        \
    if (want_hi) (hi) = RSIMD_MAX_ACC_(v_, (hi));                                        \
  } while (0)
  while (i + 4 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (from = i; i + 4 * W <= end; i += 4 * W) {
      rsimd_vf64 v0 = rsimd_vf64_loadu(x + i), v1 = rsimd_vf64_loadu(x + i + W);
      rsimd_vf64 v2 = rsimd_vf64_loadu(x + i + 2 * W), v3 = rsimd_vf64_loadu(x + i + 3 * W);
      RSIMD_MINMAX_F64_(lo0, hi0, v0);
      RSIMD_MINMAX_F64_(lo1, hi1, v1);
      RSIMD_MINMAX_F64_(lo2, hi2, v2);
      RSIMD_MINMAX_F64_(lo3, hi3, v3);
      if (check) {
        rsimd_vf64 s = rsimd_vf64_add(rsimd_vf64_add(v0, v1), rsimd_vf64_add(v2, v3));
        mnan = rsimd_mf64_or(mnan, rsimd_vf64_is_nan(s));
      }
    }
    if (stop && rsimd_mf64_any(mnan)) {
      if (rsimd_vf64_any_na(x + from, i - from)) {
        r->saw_na = 1;
        r->saw_nan = 1;
        r->count += n;
        return;
      }
      mnan = rsimd_mf64_none();
      if (!nan) nan = RSIMD_KERNEL(count_nan_f64_)(x + from, i - from) > 0;
    }
  }
  tail = stop ? i : 0; /* the blocks before tail hold no NA */
  for (; i < n; i += W) {
    /* The tail repeats element i in the inactive lanes. */
    rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, n), x + i, x[i]);
    RSIMD_MINMAX_F64_(lo0, hi0, v);
    if (check) mnan = rsimd_mf64_or(mnan, rsimd_vf64_is_nan(v));
  }
#undef RSIMD_MINMAX_F64_
  if (check && (nan || rsimd_mf64_any(mnan))) {
    /* Under stop the result is NA or NaN, whatever the count, and the
       mask is exact, as it only holds the vectors after tail. */
    if (!stop) missing = RSIMD_KERNEL(count_nan_f64_)(x, n);
    if (stop || missing) {
      r->saw_nan = 1;
      if (rsimd_mf64_any(mnan) && rsimd_vf64_any_na(x + tail, n - tail)) r->saw_na = 1;
    }
  }
  r->count += n - missing;
  /* The accumulators hold no NaN, so the lane order of the reduction does
     not matter, except for the sign of zero (always +0 for absolute
     values). */
  if (want_lo) {
    lo = rsimd_vf64_reduce_min(rsimd_vf64_min(rsimd_vf64_min(lo1, lo0), rsimd_vf64_min(lo3, lo2)));
    if (lo == 0 && !absval) lo = x[RSIMD_KERNEL(find_f64_)(x, n, 0, 0.0)];
    if (lo < r->f64) r->f64 = lo;
  }
  if (want_hi) {
    hi = rsimd_vf64_reduce_max(rsimd_vf64_max(rsimd_vf64_max(hi1, hi0), rsimd_vf64_max(hi3, hi2)));
    if (hi == 0 && !absval) hi = x[RSIMD_KERNEL(find_f64_)(x, n, 0, 0.0)];
    if (hi > r->f64_hi) r->f64_hi = hi;
  }
}
#undef RSIMD_MIN_ACC_
#undef RSIMD_MAX_ACC_

/* The kernel copy for the constant arguments. */
#define RSIMD_MINMAX_F64_EXT_(c, s, a)                                                   \
  do {                                                                                   \
    if (o->extrema == RSIMD_EXT_MIN) RSIMD_KERNEL(minmax_f64_)(x, n, c, s, a, RSIMD_EXT_MIN, r); \
    else if (o->extrema == RSIMD_EXT_MAX) RSIMD_KERNEL(minmax_f64_)(x, n, c, s, a, RSIMD_EXT_MAX, r); \
    else RSIMD_KERNEL(minmax_f64_)(x, n, c, s, a, RSIMD_EXT_BOTH, r);                   \
  } while (0)

void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_f64)(const double *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  if (!check) {
    if (absval) RSIMD_MINMAX_F64_EXT_(0, 0, 1);
    else RSIMD_MINMAX_F64_EXT_(0, 0, 0);
  } else if (o->na_rm) {
    if (absval) RSIMD_MINMAX_F64_EXT_(1, 0, 1);
    else RSIMD_MINMAX_F64_EXT_(1, 0, 0);
  } else {
    if (absval) RSIMD_MINMAX_F64_EXT_(1, 1, 1);
    else RSIMD_MINMAX_F64_EXT_(1, 1, 0);
  }
}
#undef RSIMD_MINMAX_F64_EXT_

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
    r->i64 += RSIMD_KERNEL(count_nan_f64_)(x, n);
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
#define RSIMD_SKIP_prod2_f64 1
#define RSIMD_SKIP_minmax_f64 1
#define RSIMD_SKIP_find_f64 1
#define RSIMD_SKIP_anyall_f64 1
#define RSIMD_SKIP_na_f64 1
#endif /* RSIMD_NO_F64_SIMD */

#endif /* vector tiers */

#undef RSIMD_PUT_INDEX

/* ---- Fused reductions: sum_sq, sum_abs, dot, dist, cosine, var ---------- */

/* The none tier folds with the scalar helpers of na.h, the other tiers
   with their vector forms; both read int32 operands as doubles. */
#if RSIMD_TIER_IS(none)
#define RSIMD_FOLD_ rsimd_fold_gen
#elif !defined(RSIMD_NO_F64_SIMD)
#define RSIMD_FOLD_ rsimd_vfold_gen
#endif

void RSIMD_KERNEL(sumabs_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(sumabs_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
#if RSIMD_TIER_IS(none)
  rsimd_fold_sum_i32((const int32_t *) x, n, 1, r, o);
#else
  rsimd_vfold_sum_i32((const int32_t *) x, n, 1, r, o);
#endif
}

#ifdef RSIMD_FOLD_

/* cosine folds blocks of this many elements three times while they are
   in cache; a multiple of RSIMD_PAIRWISE_LEAF. */
#define RSIMD_COSINE_BLOCK 1024

/* Element i of an operand of doubles or int32 elements, as a pointer. */
static inline const void *RSIMD_KERNEL(elt_ptr_)(const void *p, int i32, R_xlen_t i) {
  return i32 ? (const void *) ((const int *) p + i) : (const void *) ((const double *) p + i);
}

void RSIMD_KERNEL(sumsq_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                             const rsimd_opts *o);
void RSIMD_KERNEL(sumsq_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                             const rsimd_opts *o) {
  RSIMD_FOLD_(x, 0, x, 0, n, RSIMD_TERM_SQ, 0.0, r, o);
}

void RSIMD_KERNEL(sumsq_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                             const rsimd_opts *o);
void RSIMD_KERNEL(sumsq_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                             const rsimd_opts *o) {
  RSIMD_FOLD_(x, 1, x, 1, n, RSIMD_TERM_SQ, 0.0, r, o);
}

void RSIMD_KERNEL(sumabs_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(sumabs_f64)(const double *x, R_xlen_t n, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  RSIMD_FOLD_(x, 0, x, 0, n, RSIMD_TERM_ABS, 0.0, r, o);
}

void RSIMD_KERNEL(var_pass2_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                                 const rsimd_opts *o);
void RSIMD_KERNEL(var_pass2_f64)(const double *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                                 const rsimd_opts *o) {
  RSIMD_FOLD_(x, 0, x, 0, n, RSIMD_TERM_DEVSQ, c, r, o);
}

void RSIMD_KERNEL(var_pass2_i32)(const int *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                                 const rsimd_opts *o);
void RSIMD_KERNEL(var_pass2_i32)(const int *x, R_xlen_t n, double c, rsimd_reduce_result *r,
                                 const rsimd_opts *o) {
  RSIMD_FOLD_(x, 1, x, 1, n, RSIMD_TERM_DEVSQ, c, r, o);
}

/* A pair term over operands of the element types `types`. */
#define RSIMD_PAIR_FOLD_(term)                                                           \
  do {                                                                                   \
    switch (types) {                                                                     \
    case RSIMD_PAIR_I32_I32: RSIMD_FOLD_(x, 1, y, 1, n, (term), 0.0, r, o); break;       \
    case RSIMD_PAIR_F64_I32: RSIMD_FOLD_(x, 0, y, 1, n, (term), 0.0, r, o); break;       \
    default: RSIMD_FOLD_(x, 0, y, 0, n, (term), 0.0, r, o); break;                       \
    }                                                                                    \
  } while (0)

void RSIMD_KERNEL(dot_f64)(const void *x, const void *y, R_xlen_t n, int types,
                           rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(dot_f64)(const void *x, const void *y, R_xlen_t n, int types,
                           rsimd_reduce_result *r, const rsimd_opts *o) {
  RSIMD_PAIR_FOLD_(RSIMD_TERM_XY);
}

void RSIMD_KERNEL(dist_f64)(const void *x, const void *y, R_xlen_t n, int types,
                            rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(dist_f64)(const void *x, const void *y, R_xlen_t n, int types,
                            rsimd_reduce_result *r, const rsimd_opts *o) {
  RSIMD_PAIR_FOLD_(RSIMD_TERM_SQDIFF);
}
#undef RSIMD_PAIR_FOLD_

static inline void RSIMD_KERNEL(cosine_)(const void *x, const int xi32, const void *y,
                                         const int yi32, R_xlen_t n, rsimd_reduce_result *r,
                                         const rsimd_opts *o) {
  R_xlen_t i;
  for (i = 0; i < n; i += RSIMD_COSINE_BLOCK) {
    R_xlen_t len = n - i < RSIMD_COSINE_BLOCK ? n - i : RSIMD_COSINE_BLOCK;
    const void *xp = RSIMD_KERNEL(elt_ptr_)(x, xi32, i), *yp = RSIMD_KERNEL(elt_ptr_)(y, yi32, i);
    RSIMD_FOLD_(xp, xi32, yp, yi32, len, RSIMD_TERM_XY, 0.0, &r[0], o);
    RSIMD_FOLD_(xp, xi32, xp, xi32, len, RSIMD_TERM_SQ, 0.0, &r[1], o);
    RSIMD_FOLD_(yp, yi32, yp, yi32, len, RSIMD_TERM_SQ, 0.0, &r[2], o);
    /* A missing pair makes the result NA. */
    if (r[0].saw_na && !o->na_rm) break;
  }
}

void RSIMD_KERNEL(cosine_f64)(const void *x, const void *y, R_xlen_t n, int types,
                              rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(cosine_f64)(const void *x, const void *y, R_xlen_t n, int types,
                              rsimd_reduce_result *r, const rsimd_opts *o) {
  switch (types) {
  case RSIMD_PAIR_I32_I32: RSIMD_KERNEL(cosine_)(x, 1, y, 1, n, r, o); break;
  case RSIMD_PAIR_F64_I32: RSIMD_KERNEL(cosine_)(x, 0, y, 1, n, r, o); break;
  default: RSIMD_KERNEL(cosine_)(x, 0, y, 0, n, r, o); break;
  }
}

#undef RSIMD_FOLD_

#else /* RSIMD_NO_F64_SIMD: the none tier's kernels are used */
#define RSIMD_SKIP_sumsq_f64 1
#define RSIMD_SKIP_sumsq_i32 1
#define RSIMD_SKIP_sumabs_f64 1
#define RSIMD_SKIP_var_pass2_f64 1
#define RSIMD_SKIP_var_pass2_i32 1
#define RSIMD_SKIP_dot_f64 1
#define RSIMD_SKIP_dist_f64 1
#define RSIMD_SKIP_cosine_f64 1
#endif
