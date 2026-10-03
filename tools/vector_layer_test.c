/*
 * Behavioural test of the vector layer in src/kernels/common.inc.h.
 *
 * Compiled once per tier by tools/check_vector_layer.sh with
 * -DRSIMD_TIER=<tier> and that tier's flags, then run (natively or under
 * qemu). Every operation is applied through the predicated loop that
 * kernels use, for every length from 0 to N, and compared lane by lane with
 * a plain C reference. Prints a summary and exits non-zero on any mismatch.
 */

#include <float.h>
#include <stdio.h>
#include "kernels/common.inc.h"

#define N 200 /* > 3 vectors at 2048-bit SVE for 32-bit lanes */
#define SENTINEL_F64 -12345.678
#define SENTINEL_I32 0x5A5A5A5A
#define SENTINEL_I64 INT64_C(0x5A5A5A5A5A5A5A5A)

static long n_checks = 0, n_fail = 0;

static void fail(const char *what, ptrdiff_t n, ptrdiff_t i, const char *detail) {
  if (n_fail < 40) printf("  FAIL %s (n=%ld, i=%ld): %s\n", what, (long) n, (long) i, detail);
  n_fail++;
}

static double from_bits(uint64_t u) {
  double x;
  memcpy(&x, &u, sizeof x);
  return x;
}

#ifndef RSIMD_NO_F64_SIMD
static uint64_t bits(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof u);
  return u;
}

/* exact: bit-identical; otherwise NaN matches any NaN. */
static void check_f64(const char *what, ptrdiff_t n, const double *got, const double *want,
                      int exact) {
  ptrdiff_t i;
  char buf[160];
  for (i = 0; i < n; i++) {
    int ok = bits(got[i]) == bits(want[i]) || (!exact && isnan(got[i]) && isnan(want[i]));
    n_checks++;
    if (!ok) {
      snprintf(buf, sizeof buf, "got %a (%016llx) want %a (%016llx)", got[i],
               (unsigned long long) bits(got[i]), want[i], (unsigned long long) bits(want[i]));
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (bits(got[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}
#endif

static void check_i32(const char *what, ptrdiff_t n, const int32_t *got, const int32_t *want) {
  ptrdiff_t i;
  char buf[96];
  for (i = 0; i < n; i++) {
    n_checks++;
    if (got[i] != want[i]) {
      snprintf(buf, sizeof buf, "got %ld want %ld", (long) got[i], (long) want[i]);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (got[n] != SENTINEL_I32) fail(what, n, n, "wrote past the end");
}

static void check_i64(const char *what, ptrdiff_t n, const int64_t *got, const int64_t *want) {
  ptrdiff_t i;
  char buf[96];
  for (i = 0; i < n; i++) {
    n_checks++;
    if (got[i] != want[i]) {
      snprintf(buf, sizeof buf, "got %lld want %lld", (long long) got[i], (long long) want[i]);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (got[n] != SENTINEL_I64) fail(what, n, n, "wrote past the end");
}

static void check_int(const char *what, ptrdiff_t n, ptrdiff_t i, long got, long want) {
  char buf[96];
  n_checks++;
  if (got != want) {
    snprintf(buf, sizeof buf, "got %ld want %ld", got, want);
    fail(what, n, i, buf);
  }
}

/* Inputs */
static double fa[N + 1], fb[N + 1], fc[N + 1], fint[N + 1], fout[N + 1];
#ifndef RSIMD_NO_F64_SIMD
static double fref[N + 1];
#endif
static int32_t ia[N + 1], ib[N + 1], iout[N + 1], iref[N + 1];
static int64_t la[N + 1], lb[N + 1], lout[N + 1], lref[N + 1];

static uint64_t rng = UINT64_C(0x9E3779B97F4A7C15);
static uint64_t next_rand(void) {
  rng ^= rng << 13;
  rng ^= rng >> 7;
  rng ^= rng << 17;
  return rng;
}

static void init_inputs(void) {
  const double fvals[] = {
    0.0, -0.0, 1.0, -1.0, 1.5, -2.25, 1e308, -1e308, DBL_MIN, 4.9e-324, -4.9e-324,
    HUGE_VAL, -HUGE_VAL, from_bits(UINT64_C(0x7FF8000000000000)), /* NaN */
    from_bits(UINT64_C(0x7FF00000000007A2)),                       /* NA_real_ */
    from_bits(UINT64_C(0x7FF80000000007A2)),                       /* quiet NA */
    from_bits(UINT64_C(0xFFF80000000007A2)),                       /* negative NA */
    from_bits(UINT64_C(0x3FF00000000007A2)),                       /* finite, low word 1954 */
    from_bits(UINT64_C(0x7FF80000000007A3)),                       /* NaN, low word 1955 */
    3.0, 1e-300, 2.0
  };
  const int32_t ivals[] = {0, 1, -1, 2, INT32_MAX, INT32_MIN, INT32_MIN + 1, 12345, -98765,
                           0x55555555, (int32_t) 0xAAAAAAAAu, 65536};
  const int64_t lvals[] = {0, 1, -1, 2, INT64_MAX, INT64_MIN, INT64_MIN + 1,
                           INT64_C(4294967296), -INT64_C(4294967297), INT32_MIN, INT32_MAX,
                           INT64_C(0x5555555555555555)};
  const int nf = (int) (sizeof fvals / sizeof fvals[0]);
  const int ni = (int) (sizeof ivals / sizeof ivals[0]);
  const int nl = (int) (sizeof lvals / sizeof lvals[0]);
  int i;
  for (i = 0; i < N; i++) {
    if (i < 3 * nf) {
      fa[i] = fvals[i % nf];
      fb[i] = fvals[(i / nf + i * 7 + 3) % nf];
      fc[i] = fvals[(i * 5 + 1) % nf];
    } else {
      fa[i] = (double) (int64_t) (next_rand() >> 11) / 9007199254740992.0 * 200.0 - 100.0;
      fb[i] = (double) (int64_t) (next_rand() >> 11) / 9007199254740992.0 * 200.0 - 100.0;
      fc[i] = i % 5 == 0 ? fa[i] : (double) (int64_t) (next_rand() >> 11) / 1e12;
    }
    fint[i] = (double) ((int) (next_rand() % 2001) - 1000);
    if (i < 3 * ni) {
      ia[i] = ivals[i % ni];
      ib[i] = ivals[(i / ni + i * 5 + 1) % ni];
    } else {
      ia[i] = (int32_t) (uint32_t) next_rand();
      ib[i] = i % 4 == 0 ? ia[i] : (int32_t) (uint32_t) next_rand();
    }
    if (i < 3 * nl) {
      la[i] = lvals[i % nl];
      lb[i] = lvals[(i / nl + i * 5 + 1) % nl];
    } else {
      la[i] = (int64_t) next_rand();
      lb[i] = i % 4 == 0 ? la[i] : (int64_t) next_rand();
    }
  }
}

static void reset_out(void) {
  int i;
  for (i = 0; i <= N; i++) {
    fout[i] = SENTINEL_F64;
    iout[i] = SENTINEL_I32;
    lout[i] = SENTINEL_I64;
  }
}

/* The predicated loop every kernel uses: full vectors, then one predicated
   step. EXPR computes the output vector from the loaded inputs x, y, z. */
#define LOOP(L, P, V, LOADS, EXPR, OUTPTR)                                       \
  do {                                                                         \
    ptrdiff_t i = 0;                                                           \
    for (; i + (L) <= n; i += (L)) {                                           \
      rsimd_##P pg = rsimd_##P##_true();                                       \
      LOADS;                                                                   \
      rsimd_##V##_storeu((OUTPTR) + i, (EXPR));                                \
      (void) pg;                                                               \
    }                                                                          \
    if (i < n) {                                                               \
      rsimd_##P pg = rsimd_##P##_while(i, n);                                  \
      LOADS;                                                                   \
      rsimd_##V##_storeu_p(pg, (OUTPTR) + i, (EXPR));                          \
    }                                                                          \
  } while (0)

#ifndef RSIMD_NO_F64_SIMD
#define F64_LOADS                                                              \
  rsimd_vf64 x = rsimd_vf64_loadu_p(pg, fa + i, 7.0);                          \
  rsimd_vf64 y = rsimd_vf64_loadu_p(pg, fb + i, 7.0);                          \
  rsimd_vf64 z = rsimd_vf64_loadu_p(pg, fc + i, 7.0);                          \
  (void) x; (void) y; (void) z
#define F64_OP(name, EXPR, REF, exact)                                         \
  do {                                                                         \
    reset_out();                                                               \
    LOOP(RSIMD_LANES_64, p64, vf64, F64_LOADS, EXPR, fout);                    \
    for (j = 0; j < n; j++) {                                                  \
      double a = fa[j], b = fb[j], c = fc[j];                                  \
      (void) a; (void) b; (void) c;                                            \
      fref[j] = (REF);                                                         \
    }                                                                          \
    check_f64(name, n, fout, fref, exact);                                     \
  } while (0)
/* A mask is checked through blend(0, 1, m). */
#define F64_MASK(name, MEXPR, REF)                                             \
  F64_OP(name, rsimd_vf64_blend(rsimd_vf64_zero(), rsimd_vf64_set1(1.0), (MEXPR)), \
         (REF) ? 1.0 : 0.0, 1)

static int is_na_ref(double x) {
  return isnan(x) && (bits(x) & 0xFFFFFFFFu) == 1954;
}

static void test_f64(ptrdiff_t n) {
  ptrdiff_t j;
  F64_OP("f64 copy", x, a, 1);
  F64_OP("f64 add", rsimd_vf64_add(x, y), a + b, 0);
  F64_OP("f64 sub", rsimd_vf64_sub(x, y), a - b, 0);
  F64_OP("f64 mul", rsimd_vf64_mul(x, y), a * b, 0);
  F64_OP("f64 div", rsimd_vf64_div(x, y), a / b, 0);
  F64_OP("f64 fma", rsimd_vf64_fma(x, y, z), fma(a, b, c), 0);
  F64_OP("f64 min", rsimd_vf64_min(x, y), a < b ? a : b, 1);
  F64_OP("f64 max", rsimd_vf64_max(x, y), a > b ? a : b, 1);
  F64_OP("f64 abs", rsimd_vf64_abs(x), from_bits(bits(a) & ~(UINT64_C(1) << 63)), 1);
  F64_OP("f64 neg", rsimd_vf64_neg(x), from_bits(bits(a) ^ (UINT64_C(1) << 63)), 1);
  F64_OP("f64 sqrt", rsimd_vf64_sqrt(x), sqrt(a), 0);
  F64_OP("f64 and", rsimd_vf64_and(x, y), from_bits(bits(a) & bits(b)), 1);
  F64_OP("f64 or", rsimd_vf64_or(x, y), from_bits(bits(a) | bits(b)), 1);
  F64_OP("f64 xor", rsimd_vf64_xor(x, y), from_bits(bits(a) ^ bits(b)), 1);
  F64_OP("f64 andnot", rsimd_vf64_andnot(x, y), from_bits(~bits(a) & bits(b)), 1);
  F64_OP("f64 set1", rsimd_vf64_set1(-3.5), -3.5, 1);
  F64_OP("f64 blend", rsimd_vf64_blend(x, y, rsimd_vf64_cmp_lt(z, x)), c < a ? b : a, 1);
  F64_OP("f64 as_vi64", rsimd_vi64_as_vf64(rsimd_vi64_add(rsimd_vf64_as_vi64(x),
                                                         rsimd_vi64_set1(1))),
         from_bits(bits(a) + 1), 1);
  F64_MASK("f64 cmp_eq", rsimd_vf64_cmp_eq(x, z), a == c);
  F64_MASK("f64 cmp_ne", rsimd_vf64_cmp_ne(x, z), !(a == c));
  F64_MASK("f64 cmp_lt", rsimd_vf64_cmp_lt(x, y), a < b);
  F64_MASK("f64 cmp_le", rsimd_vf64_cmp_le(x, y), a <= b);
  F64_MASK("f64 cmp_gt", rsimd_vf64_cmp_gt(x, y), a > b);
  F64_MASK("f64 cmp_ge", rsimd_vf64_cmp_ge(x, y), a >= b);
  F64_MASK("f64 is_nan", rsimd_vf64_is_nan(x), isnan(a));
  F64_MASK("f64 is_na", rsimd_vf64_is_na(x), is_na_ref(a));
  F64_MASK("mf64 and", rsimd_mf64_and(rsimd_vf64_cmp_lt(x, y), rsimd_vf64_cmp_lt(z, x)),
           a < b && c < a);
  F64_MASK("mf64 or", rsimd_mf64_or(rsimd_vf64_cmp_lt(x, y), rsimd_vf64_cmp_lt(z, x)),
           a < b || c < a);
  F64_MASK("mf64 not", rsimd_mf64_not(rsimd_vf64_cmp_lt(x, y)), !(a < b));
  F64_MASK("mf64 andnot", rsimd_mf64_andnot(rsimd_vf64_cmp_lt(x, y), rsimd_vf64_cmp_lt(z, x)),
           !(a < b) && c < a);
  F64_MASK("mf64 via mi64",
           rsimd_mi64_to_mf64(rsimd_mf64_to_mi64(rsimd_vf64_is_nan(x))), isnan(a));
}

/* Horizontal operations, on full vectors only. */
static void test_f64_horizontal(void) {
  ptrdiff_t L = RSIMD_LANES_64, i, j;
  for (i = 0; i + L <= N; i += L) {
    rsimd_vf64 x = rsimd_vf64_loadu(fa + i), y = rsimd_vf64_loadu(fb + i);
    rsimd_mf64 m = rsimd_vf64_cmp_lt(x, y);
    double s = 0, lo = fint[i], hi = fint[i];
    long cnt = 0;
    for (j = 0; j < L; j++) cnt += fa[i + j] < fb[i + j];
    check_int("mf64 count", L, i, rsimd_mf64_count(m), cnt);
    check_int("mf64 any", L, i, rsimd_mf64_any(m), cnt > 0);
    check_int("mf64 all", L, i, rsimd_mf64_all(m), cnt == L);
    for (j = 0; j < L; j++) {
      s += fint[i + j];
      lo = fint[i + j] < lo ? fint[i + j] : lo;
      hi = fint[i + j] > hi ? fint[i + j] : hi;
    }
    x = rsimd_vf64_loadu(fint + i);
    check_int("f64 reduce_add", L, i, (long) rsimd_vf64_reduce_add(x), (long) s);
    check_int("f64 reduce_min", L, i, (long) rsimd_vf64_reduce_min(x), (long) lo);
    check_int("f64 reduce_max", L, i, (long) rsimd_vf64_reduce_max(x), (long) hi);
  }
}
#endif /* RSIMD_NO_F64_SIMD */

/* Integers: the same pattern for 32- and 64-bit lanes. */
#define INT_LOADS(V, A, B, FILL)                                               \
  rsimd_##V x = rsimd_##V##_loadu_p(pg, A + i, FILL);                          \
  rsimd_##V y = rsimd_##V##_loadu_p(pg, B + i, FILL);                          \
  (void) x; (void) y

#define INT_OP(V, P, L, T, U, A, B, OUT, REF_ARR, CHECK, name, EXPR, REF)       \
  do {                                                                         \
    reset_out();                                                               \
    LOOP(L, P, V, INT_LOADS(V, A, B, 7), EXPR, OUT);                           \
    for (j = 0; j < n; j++) {                                                  \
      T a = A[j], b = B[j];                                                    \
      (void) a; (void) b;                                                      \
      REF_ARR[j] = (REF);                                                      \
    }                                                                          \
    CHECK(name, n, OUT, REF_ARR);                                              \
  } while (0)

#define I32_OP(name, EXPR, REF)                                                \
  INT_OP(vi32, p32, RSIMD_LANES_32, int32_t, uint32_t, ia, ib, iout, iref, check_i32, \
         name, EXPR, REF)
#define I64_OP(name, EXPR, REF)                                                \
  INT_OP(vi64, p64, RSIMD_LANES_64, int64_t, uint64_t, la, lb, lout, lref, check_i64, \
         name, EXPR, REF)
#define I32_MASK(name, MEXPR, REF)                                             \
  I32_OP(name, rsimd_vi32_blend(rsimd_vi32_zero(), rsimd_vi32_set1(1), (MEXPR)), (REF) ? 1 : 0)
#define I64_MASK(name, MEXPR, REF)                                             \
  I64_OP(name, rsimd_vi64_blend(rsimd_vi64_zero(), rsimd_vi64_set1(1), (MEXPR)), (REF) ? 1 : 0)

static void test_int(ptrdiff_t n) {
  ptrdiff_t j;
  I32_OP("i32 copy", x, a);
  I32_OP("i32 add", rsimd_vi32_add(x, y), (int32_t) ((uint32_t) a + (uint32_t) b));
  I32_OP("i32 sub", rsimd_vi32_sub(x, y), (int32_t) ((uint32_t) a - (uint32_t) b));
  I32_OP("i32 mul", rsimd_vi32_mul(x, y), (int32_t) ((uint32_t) a * (uint32_t) b));
  I32_OP("i32 min", rsimd_vi32_min(x, y), a < b ? a : b);
  I32_OP("i32 max", rsimd_vi32_max(x, y), a > b ? a : b);
  I32_OP("i32 and", rsimd_vi32_and(x, y), a & b);
  I32_OP("i32 or", rsimd_vi32_or(x, y), a | b);
  I32_OP("i32 xor", rsimd_vi32_xor(x, y), a ^ b);
  I32_OP("i32 andnot", rsimd_vi32_andnot(x, y), ~a & b);
  I32_OP("i32 set1", rsimd_vi32_set1(-42), -42);
  I32_OP("i32 blend", rsimd_vi32_blend(x, y, rsimd_vi32_cmp_gt(x, y)), a > b ? b : a);
  I32_MASK("i32 cmp_eq", rsimd_vi32_cmp_eq(x, y), a == b);
  I32_MASK("i32 cmp_gt", rsimd_vi32_cmp_gt(x, y), a > b);
  I32_MASK("i32 cmp_lt", rsimd_vi32_cmp_lt(x, y), a < b);
  I32_MASK("i32 is_na", rsimd_vi32_is_na(x), a == INT32_MIN);
  I32_MASK("mi32 and", rsimd_mi32_and(rsimd_vi32_cmp_gt(x, y), rsimd_vi32_is_na(y)),
           a > b && b == INT32_MIN);
  I32_MASK("mi32 or", rsimd_mi32_or(rsimd_vi32_cmp_gt(x, y), rsimd_vi32_is_na(y)),
           a > b || b == INT32_MIN);
  I32_MASK("mi32 not", rsimd_mi32_not(rsimd_vi32_cmp_gt(x, y)), !(a > b));
  I32_MASK("mi32 andnot", rsimd_mi32_andnot(rsimd_vi32_cmp_gt(x, y), rsimd_vi32_cmp_lt(x, y)),
           !(a > b) && a < b);

  I64_OP("i64 copy", x, a);
  I64_OP("i64 add", rsimd_vi64_add(x, y), (int64_t) ((uint64_t) a + (uint64_t) b));
  I64_OP("i64 sub", rsimd_vi64_sub(x, y), (int64_t) ((uint64_t) a - (uint64_t) b));
  I64_OP("i64 mul", rsimd_vi64_mul(x, y), (int64_t) ((uint64_t) a * (uint64_t) b));
  I64_OP("i64 min", rsimd_vi64_min(x, y), a < b ? a : b);
  I64_OP("i64 max", rsimd_vi64_max(x, y), a > b ? a : b);
  I64_OP("i64 and", rsimd_vi64_and(x, y), a & b);
  I64_OP("i64 or", rsimd_vi64_or(x, y), a | b);
  I64_OP("i64 xor", rsimd_vi64_xor(x, y), a ^ b);
  I64_OP("i64 andnot", rsimd_vi64_andnot(x, y), ~a & b);
  I64_OP("i64 set1", rsimd_vi64_set1(INT64_MIN + 3), INT64_MIN + 3);
  I64_OP("i64 blend", rsimd_vi64_blend(x, y, rsimd_vi64_cmp_gt(x, y)), a > b ? b : a);
  I64_MASK("i64 cmp_eq", rsimd_vi64_cmp_eq(x, y), a == b);
  I64_MASK("i64 cmp_gt", rsimd_vi64_cmp_gt(x, y), a > b);
  I64_MASK("i64 cmp_lt", rsimd_vi64_cmp_lt(x, y), a < b);
  I64_MASK("i64 is_na", rsimd_vi64_is_na(x), a == INT64_MIN);
  I64_MASK("mi64 and", rsimd_mi64_and(rsimd_vi64_cmp_gt(x, y), rsimd_vi64_is_na(y)),
           a > b && b == INT64_MIN);
  I64_MASK("mi64 or", rsimd_mi64_or(rsimd_vi64_cmp_gt(x, y), rsimd_vi64_is_na(y)),
           a > b || b == INT64_MIN);
  I64_MASK("mi64 not", rsimd_mi64_not(rsimd_vi64_cmp_gt(x, y)), !(a > b));
  I64_MASK("mi64 andnot", rsimd_mi64_andnot(rsimd_vi64_cmp_gt(x, y), rsimd_vi64_cmp_lt(x, y)),
           !(a > b) && a < b);
}

static void test_int_horizontal(void) {
  ptrdiff_t L = RSIMD_LANES_32, i, j;
  for (i = 0; i + L <= N; i += L) {
    rsimd_vi32 x = rsimd_vi32_loadu(ia + i), y = rsimd_vi32_loadu(ib + i);
    rsimd_mi32 m = rsimd_vi32_cmp_gt(x, y);
    int64_t s = 0;
    int32_t lo = ia[i], hi = ia[i];
    long cnt = 0;
    for (j = 0; j < L; j++) {
      cnt += ia[i + j] > ib[i + j];
      s += ia[i + j];
      lo = ia[i + j] < lo ? ia[i + j] : lo;
      hi = ia[i + j] > hi ? ia[i + j] : hi;
    }
    check_int("mi32 count", L, i, rsimd_mi32_count(m), cnt);
    check_int("mi32 any", L, i, rsimd_mi32_any(m), cnt > 0);
    check_int("mi32 all", L, i, rsimd_mi32_all(m), cnt == L);
    check_int("mi32 all (true)", L, i, rsimd_mi32_all(rsimd_vi32_cmp_eq(x, x)), 1);
    check_int("mi32 any (false)", L, i, rsimd_mi32_any(rsimd_vi32_cmp_gt(x, x)), 0);
    n_checks++;
    if (rsimd_vi32_reduce_add(x) != s) fail("i32 reduce_add", L, i, "mismatch");
    check_int("i32 reduce_min", L, i, rsimd_vi32_reduce_min(x), lo);
    check_int("i32 reduce_max", L, i, rsimd_vi32_reduce_max(x), hi);
  }
  L = RSIMD_LANES_64;
  for (i = 0; i + L <= N; i += L) {
    rsimd_vi64 x = rsimd_vi64_loadu(la + i), y = rsimd_vi64_loadu(lb + i);
    rsimd_mi64 m = rsimd_vi64_cmp_gt(x, y);
    uint64_t s = 0;
    int64_t lo = la[i], hi = la[i];
    long cnt = 0;
    for (j = 0; j < L; j++) {
      cnt += la[i + j] > lb[i + j];
      s += (uint64_t) la[i + j];
      lo = la[i + j] < lo ? la[i + j] : lo;
      hi = la[i + j] > hi ? la[i + j] : hi;
    }
    check_int("mi64 count", L, i, rsimd_mi64_count(m), cnt);
    check_int("mi64 any", L, i, rsimd_mi64_any(m), cnt > 0);
    check_int("mi64 all", L, i, rsimd_mi64_all(m), cnt == L);
    check_int("mi64 all (true)", L, i, rsimd_mi64_all(rsimd_vi64_cmp_eq(x, x)), 1);
    n_checks += 3;
    if (rsimd_vi64_reduce_add(x) != (int64_t) s) fail("i64 reduce_add", L, i, "mismatch");
    if (rsimd_vi64_reduce_min(x) != lo) fail("i64 reduce_min", L, i, "mismatch");
    if (rsimd_vi64_reduce_max(x) != hi) fail("i64 reduce_max", L, i, "mismatch");
  }
}

/* Predicates and fill values. */
static void test_predicates(ptrdiff_t n) {
  ptrdiff_t L64 = RSIMD_LANES_64, L32 = RSIMD_LANES_32, i;
  for (i = 0; i < n; i += L64) {
    ptrdiff_t want = n - i < L64 ? n - i : L64;
    check_int("p64 count", n, i, rsimd_p64_count(rsimd_p64_while(i, n)), (long) want);
#ifndef RSIMD_NO_F64_SIMD
    {
      /* Inactive lanes take the fill value. */
      rsimd_vf64 v = rsimd_vf64_loadu_p(rsimd_p64_while(i, n), fint + i, 0.0);
      double s = 0;
      ptrdiff_t j;
      for (j = 0; j < want; j++) s += fint[i + j];
      check_int("f64 loadu_p fill", n, i, (long) rsimd_vf64_reduce_add(v), (long) s);
    }
#endif
  }
  for (i = 0; i < n; i += L32) {
    ptrdiff_t want = n - i < L32 ? n - i : L32;
    rsimd_vi32 v = rsimd_vi32_loadu_p(rsimd_p32_while(i, n), ia + i, 0);
    int64_t s = 0;
    ptrdiff_t j;
    for (j = 0; j < want; j++) s += ia[i + j];
    check_int("p32 count", n, i, rsimd_p32_count(rsimd_p32_while(i, n)), (long) want);
    n_checks++;
    if (rsimd_vi32_reduce_add(v) != s) fail("i32 loadu_p fill", n, i, "mismatch");
  }
  check_int("p64 true", n, 0, rsimd_p64_count(rsimd_p64_true()), (long) L64);
  check_int("p32 true", n, 0, rsimd_p32_count(rsimd_p32_true()), (long) L32);
  check_int("p64 empty", n, n, rsimd_p64_count(rsimd_p64_while(n, n)), 0);
}

int main(void) {
  ptrdiff_t n;
  init_inputs();
  for (n = 0; n <= N; n++) {
#ifndef RSIMD_NO_F64_SIMD
    test_f64(n);
#endif
    test_int(n);
    test_predicates(n);
  }
#ifndef RSIMD_NO_F64_SIMD
  test_f64_horizontal();
#endif
  test_int_horizontal();
  printf("tier %s: lanes64=%ld lanes32=%ld, %ld checks, %ld failures\n", RSIMD_TIER_STRING,
         (long) RSIMD_LANES_64, (long) RSIMD_LANES_32, n_checks, n_fail);
  return n_fail != 0;
}
