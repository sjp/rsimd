/*
 * Behavioural test of the vector layer in src/kernels/common.inc.h and of
 * the vector forms of the NA, precision and overflow helpers in src/na.h,
 * which must agree with their scalar forms.
 *
 * Compiled once per tier by tools/check_vector_layer.sh with
 * -DRSIMD_TIER=<tier> and that tier's flags, then run (natively or under
 * qemu). Every operation is applied through the predicated loop that
 * kernels use, for every length from 0 to N, and compared lane by lane with
 * a plain C reference. The reduction and scan kernels
 * (src/kernels/reduce.inc.c, src/kernels/scan.inc.c) are run on inputs full
 * of ties, signed zeros and missing values, whole and split into two
 * chunks, and compared with plain C references too. Prints a summary and
 * exits non-zero on any mismatch.
 */

#include <float.h>
#include <stdio.h>
#include "kernels/common.inc.h"
#include "kernels/reduce.inc.c"
#include "kernels/scan.inc.c"

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
static double fconv[N + 1]; /* within the int32 range, for storeu_i32 */
static int32_t ia[N + 1], ib[N + 1], iout[N + 1], iref[N + 1];
static int32_t ismall[N + 1], ismall2[N + 1]; /* near the int32 multiply limits */
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
    fconv[i] = fint[i] * 2147483.0 + (double) (i % 7) * (i % 2 ? -0.37 : 0.37);
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
  {
    const int32_t edges[] = {46340, 46341, -46341, 65536, -65536, 32768, -32768, 0, 1, -1,
                             INT32_MIN, INT32_MAX, -INT32_MAX};
    const int ne = (int) (sizeof edges / sizeof edges[0]);
    for (i = 0; i < N; i++) {
      ismall[i] = i < 2 * ne ? edges[i % ne] : (int32_t) (next_rand() % 131073) - 65536;
      ismall2[i] =
        i < 2 * ne ? edges[(i * 5 + 2) % ne] : (int32_t) (next_rand() % 131073) - 65536;
    }
    /* Results of exactly INT32_MIN, which is NA: add, sub, mul. */
    ismall[0] = -INT32_MAX;
    ismall2[0] = -1;
    ismall[1] = -INT32_MAX;
    ismall2[1] = 1;
    ismall[2] = -65536;
    ismall2[2] = 32768;
    ismall[3] = INT32_MAX;
    ismall2[3] = 1;
  }
  fconv[0] = 2147483647.0;
  fconv[1] = -2147483648.0;
  fconv[2] = -2147483647.9;
  fconv[3] = -0.5;
  fconv[4] = 0.999;
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
  F64_OP("f64 floor", rsimd_vf64_floor(x), floor(a), 0);
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
  I32_OP("i32 mulhi", rsimd_vi32_mulhi(x, y),
         (int32_t) (uint32_t) ((uint64_t) ((int64_t) a * b) >> 32));
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

/* int32 elements <-> 64-bit lanes: RSIMD_LANES_64 elements per step. */
#define CONV_LOOP(STEP_FULL, STEP_PART)                                        \
  do {                                                                         \
    ptrdiff_t i = 0;                                                           \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                     \
      STEP_FULL;                                                               \
    }                                                                          \
    if (i < n) {                                                               \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                    \
      STEP_PART;                                                               \
    }                                                                          \
  } while (0)

static void test_convert(ptrdiff_t n) {
  ptrdiff_t j;
  reset_out();
  CONV_LOOP(rsimd_vi64_storeu(lout + i, rsimd_vi64_loadu_i32(ia + i)),
            rsimd_vi64_storeu_p(pg, lout + i, rsimd_vi64_loadu_i32_p(pg, ia + i, 7)));
  for (j = 0; j < n; j++) lref[j] = ia[j];
  check_i64("i64 loadu_i32", n, lout, lref);

  reset_out();
  CONV_LOOP(rsimd_vi64_storeu_i32(iout + i, rsimd_vi64_loadu(la + i)),
            rsimd_vi64_storeu_i32_p(pg, iout + i, rsimd_vi64_loadu_p(pg, la + i, 7)));
  for (j = 0; j < n; j++) iref[j] = (int32_t) (uint32_t) (uint64_t) la[j];
  check_i32("i64 storeu_i32", n, iout, iref);

  /* Fill values of inactive lanes. */
  if (n % RSIMD_LANES_64 != 0) {
    ptrdiff_t i = n - n % RSIMD_LANES_64;
    rsimd_vi64 v = rsimd_vi64_loadu_i32_p(rsimd_p64_while(i, n), ia + i, -5);
    int64_t s = rsimd_vi64_reduce_add(v), want = 0;
    for (j = i; j < n; j++) want += ia[j];
    want += -5 * (RSIMD_LANES_64 - (n - i));
    n_checks++;
    if (s != want) fail("i64 loadu_i32_p fill", n, i, "mismatch");
  }

#ifndef RSIMD_NO_F64_SIMD
  reset_out();
  CONV_LOOP(rsimd_vf64_storeu(fout + i, rsimd_vf64_loadu_i32(ia + i)),
            rsimd_vf64_storeu_p(pg, fout + i, rsimd_vf64_loadu_i32_p(pg, ia + i, 7)));
  for (j = 0; j < n; j++) fref[j] = (double) ia[j];
  check_f64("f64 loadu_i32", n, fout, fref, 1);

  reset_out();
  CONV_LOOP(rsimd_vf64_storeu_i32(iout + i, rsimd_vf64_loadu(fconv + i)),
            rsimd_vf64_storeu_i32_p(pg, iout + i, rsimd_vf64_loadu_p(pg, fconv + i, 0.0)));
  for (j = 0; j < n; j++) iref[j] = (int32_t) fconv[j];
  check_i32("f64 storeu_i32", n, iout, iref);
#endif
}

/* ---- na.h: vector forms against scalar forms ---------------------------- */

static int same_result(const rsimd_reduce_result *a, const rsimd_reduce_result *b) {
  return a->saw_na == b->saw_na && a->saw_nan == b->saw_nan && a->count == b->count &&
         a->overflow == b->overflow && a->i64 == b->i64;
}

/* Checked (op 0-2: add, sub, mul) or wrapping (op 3-5) arithmetic through
   the vector helpers of na.h. */
static rsimd_vi32 arith_lanes(int op, int check, rsimd_vi32 x, rsimd_vi32 y) {
  rsimd_vi32 r = op % 3 == 0 ? rsimd_vi32_add(x, y)
                 : op % 3 == 1 ? rsimd_vi32_sub(x, y)
                               : rsimd_vi32_mul(x, y);
  rsimd_mi32 bad = rsimd_mi32_none();
  if (op == 0) bad = rsimd_vi32_add_ovf(x, y, r);
  else if (op == 1) bad = rsimd_vi32_sub_ovf(x, y, r);
  else if (op == 2) bad = rsimd_vi32_mul_ovf(x, y, r);
  if (check) bad = rsimd_mi32_or(bad, rsimd_vi32_na2(x, y));
  return rsimd_vi32_set_na(r, bad);
}

/* 1 in the lanes where checked op 0-2 overflows (NA operands excluded). */
static rsimd_vi32 ovf_lanes(int op, rsimd_vi32 x, rsimd_vi32 y) {
  rsimd_vi32 r = op == 0 ? rsimd_vi32_add(x, y) : op == 1 ? rsimd_vi32_sub(x, y)
                                                          : rsimd_vi32_mul(x, y);
  rsimd_mi32 bad = op == 0 ? rsimd_vi32_add_ovf(x, y, r)
                   : op == 1 ? rsimd_vi32_sub_ovf(x, y, r)
                             : rsimd_vi32_mul_ovf(x, y, r);
  bad = rsimd_mi32_andnot(rsimd_vi32_na2(x, y), bad);
  return rsimd_vi32_blend(rsimd_vi32_zero(), rsimd_vi32_set1(1), bad);
}

static void test_na_int(ptrdiff_t n) {
  static const rsimd_opts opts[3] = {{0, 1, 0}, {1, 1, 0}, {0, 0, 0}};
  ptrdiff_t j;
  int k;
  char buf[96];
  for (k = 0; k < 3; k++) {
    rsimd_reduce_result a, b;
    memset(&a, 0, sizeof a);
    memset(&b, 0, sizeof b);
    rsimd_fold_sum_i32(ia, n, 0, &a, &opts[k]);
    rsimd_vfold_sum_i32(ia, n, 0, &b, &opts[k]);
    n_checks++;
    if (!same_result(&a, &b)) {
      snprintf(buf, sizeof buf, "opts %d: sum %lld/%lld count %ld/%ld na %d/%d", k,
               (long long) a.i64, (long long) b.i64, (long) a.count, (long) b.count, a.saw_na,
               b.saw_na);
      fail("na sum_i32", n, 0, buf);
    }
  }

  /* Logical flags: identical without early exit; with it, the same answer. */
  for (k = 0; k < 3; k++) {
    int stop;
    for (stop = RSIMD_STOP_NONE; stop <= RSIMD_STOP_FALSE; stop++) {
      rsimd_reduce_result a, b;
      memset(&a, 0, sizeof a);
      memset(&b, 0, sizeof b);
      rsimd_fold_lgl(ib, n, stop, &a, &opts[k]);
      rsimd_vfold_lgl(ib, n, stop, &b, &opts[k]);
      n_checks++;
      if (stop == RSIMD_STOP_NONE
            ? (a.any_true != b.any_true || a.any_false != b.any_false || a.saw_na != b.saw_na)
            : stop == RSIMD_STOP_TRUE ? a.any_true != b.any_true
                                      : a.any_false != b.any_false) {
        snprintf(buf, sizeof buf, "opts %d stop %d", k, stop);
        fail("na lgl flags", n, 0, buf);
      }
    }
  }

  /* Checked and wrapping arithmetic, lane by lane. */
  for (k = 0; k < 2; k++) {
    const int check = k == 0;
    const int32_t *xs = n % 2 ? ia : ismall, *ys = n % 2 ? ib : ismall2;
    int op;
    for (op = 0; op < 6; op++) {
      int want_ovf = 0;
      reset_out();
      LOOP(RSIMD_LANES_32, p32, vi32, INT_LOADS(vi32, xs, ys, 0), arith_lanes(op, check, x, y),
           iout);
      for (j = 0; j < n; j++) {
        int32_t a = xs[j], b = ys[j];
        switch (op) {
        case 0: iref[j] = rsimd_add_i32(a, b, check, &want_ovf); break;
        case 1: iref[j] = rsimd_sub_i32(a, b, check, &want_ovf); break;
        case 2: iref[j] = rsimd_mul_i32(a, b, check, &want_ovf); break;
        case 3: iref[j] = rsimd_add_wrap_i32(a, b, check); break;
        case 4: iref[j] = rsimd_sub_wrap_i32(a, b, check); break;
        default: iref[j] = rsimd_mul_wrap_i32(a, b, check); break;
        }
      }
      snprintf(buf, sizeof buf, "na arith op %d check %d", op, check);
      check_i32(buf, n, iout, iref);
      if (op < 3) {
        reset_out();
        LOOP(RSIMD_LANES_32, p32, vi32, INT_LOADS(vi32, xs, ys, 0), ovf_lanes(op, x, y), iout);
        for (j = 0; j < n; j++) {
          int ovf = 0;
          if (op == 0) rsimd_add_i32(xs[j], ys[j], 1, &ovf);
          else if (op == 1) rsimd_sub_i32(xs[j], ys[j], 1, &ovf);
          else rsimd_mul_i32(xs[j], ys[j], 1, &ovf);
          iref[j] = ovf;
        }
        snprintf(buf, sizeof buf, "na overflow flag op %d", op);
        check_i32(buf, n, iout, iref);
      }
    }
  }
  I32_OP("na neg_wrap", rsimd_vi32_neg_wrap(x), rsimd_neg_i32(a));
  I32_OP("na abs_wrap", rsimd_vi32_abs_wrap(x), rsimd_abs_i32(a));

#ifndef RSIMD_NO_F64_SIMD
  /* %/% and %% through double lanes. */
  for (k = 0; k < 2; k++) {
    int mod;
    for (mod = 0; mod < 2; mod++) {
      const int32_t *ys = k ? ib : ismall2;
      reset_out();
      CONV_LOOP(rsimd_vf64_storeu_i32(iout + i, rsimd_vf64_intdiv(rsimd_vf64_loadu_i32(ia + i),
                                                                  rsimd_vf64_loadu_i32(ys + i),
                                                                  mod, 1)),
                rsimd_vf64_storeu_i32_p(pg, iout + i,
                                        rsimd_vf64_intdiv(rsimd_vf64_loadu_i32_p(pg, ia + i, 0),
                                                          rsimd_vf64_loadu_i32_p(pg, ys + i, 1),
                                                          mod, 1)));
      for (j = 0; j < n; j++) iref[j] = rsimd_intdiv_i32(ia[j], ys[j], mod, 1);
      check_i32(mod ? "na mod" : "na idiv", n, iout, iref);
    }
  }
#endif
}

#ifndef RSIMD_NO_F64_SIMD
/* Missing values are found in every lane, and integer-valued data sums
   exactly in every mode, so vector and scalar folds must agree exactly. */
static void test_na_fold(ptrdiff_t n) {
  static const rsimd_opts base[3] = {{0, 1, 0}, {1, 1, 0}, {0, 0, 0}};
  static double xs[N + 1], ys[N + 1];
  ptrdiff_t j;
  int k, prec, term;
  char buf[160];
  for (j = 0; j < n; j++) {
    xs[j] = fint[j];
    ys[j] = fint[(j * 7 + 3) % N];
    if (j % 13 == 5) xs[j] = from_bits(UINT64_C(0x7FF00000000007A2)); /* NA */
    if (j % 17 == 9) ys[j] = from_bits(UINT64_C(0x7FF8000000000000)); /* NaN */
  }
  for (k = 0; k < 3; k++) {
    for (prec = 0; prec < 3; prec++) {
      for (term = 0; term < 4; term++) {
        rsimd_opts o = base[k];
        rsimd_reduce_result a, b;
        double va, vb;
        o.precision = prec;
        memset(&a, 0, sizeof a);
        memset(&b, 0, sizeof b);
        rsimd_fold_f64(xs, ys, n, term, &a, &o);
        switch (term) {
        case RSIMD_TERM_X: rsimd_vfold_f64(xs, ys, n, RSIMD_TERM_X, &b, &o); break;
        case RSIMD_TERM_SQ: rsimd_vfold_f64(xs, ys, n, RSIMD_TERM_SQ, &b, &o); break;
        case RSIMD_TERM_ABS: rsimd_vfold_f64(xs, ys, n, RSIMD_TERM_ABS, &b, &o); break;
        default: rsimd_vfold_f64(xs, ys, n, RSIMD_TERM_XY, &b, &o); break;
        }
        va = rsimd_reduce_value(&a, prec);
        vb = rsimd_reduce_value(&b, prec);
        n_checks++;
        if (!same_result(&a, &b) || !(bits(va) == bits(vb) || (isnan(va) && isnan(vb)))) {
          snprintf(buf, sizeof buf, "opts %d prec %d term %d: %g/%g count %ld/%ld na %d/%d",
                   k, prec, term, va, vb, (long) a.count, (long) b.count, a.saw_na, b.saw_na);
          fail("na fold_f64", n, 0, buf);
        }
      }
    }
  }
  /* Elementwise NA merge. */
  F64_OP("na merge", rsimd_vf64_na_merge(rsimd_vf64_add(x, y), x, y),
         rsimd_na_merge_f64(a + b, a, b), 0);
}

/* Cancellation that every lane width leaves to accumulator 0, lane 0. */
static void test_na_cancel(void) {
  static double x[257];
  static const double want[3] = {0.0, 0.0, 1.0};
  rsimd_opts o = {0, 1, 0};
  int prec;
  memset(x, 0, sizeof x);
  x[0] = 1e16;
  x[128] = 1.0;
  x[256] = -1e16;
  for (prec = 0; prec < 3; prec++) {
    rsimd_reduce_result r;
    memset(&r, 0, sizeof r);
    o.precision = prec;
    rsimd_vfold_f64(x, x, 257, RSIMD_TERM_X, &r, &o);
    check_int("na cancel", 257, prec, (long) rsimd_reduce_value(&r, prec), (long) want[prec]);
  }
}
#endif

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


/* ---- Reduction kernels ---------------------------------------------------- */

static double rd[N + 1];
static int32_t ri[N + 1];
static Rbyte rb[N + 1];
static double ridx[N + 1];
static int iidx[N + 1];

/* variant 0: ties and signed zeros, no missing values; 1: some NaN and NA;
   2: positive values and zeros; 3: all missing. Doubles are exact in
   products (powers of two) and in sums of deviations from 0.5. */
static void fill_reduce(int variant) {
  const double dv[] = {0.0, -0.0, 1.0, -1.0, 2.0, 0.5, -0.5, 4.0, -2.0, HUGE_VAL, -HUGE_VAL};
  const double pos[] = {0.0, -0.0, 1.0, 2.0, 0.5};
  const int32_t iv[] = {0, 1, -1, 2, -2, 5, INT32_MAX, -INT32_MAX};
  const double nan = from_bits(UINT64_C(0x7FF8000000000000));
  const double na = from_bits(UINT64_C(0x7FF00000000007A2));
  int i;
  for (i = 0; i < N; i++) {
    uint64_t u = next_rand();
    if (variant == 3) {
      rd[i] = u % 2 ? nan : na;
      ri[i] = RSIMD_NA_I32;
    } else if (variant == 2) {
      rd[i] = pos[u % 5];
      ri[i] = (int32_t) (u % 7);
    } else {
      rd[i] = dv[u % 11];
      ri[i] = iv[(u >> 8) % 8];
      if (variant == 1 && (u >> 20) % 16 == 0) rd[i] = (u >> 30) % 2 ? nan : na;
      if (variant == 1 && (u >> 24) % 16 == 0) ri[i] = RSIMD_NA_I32;
    }
    rb[i] = (Rbyte) (variant == 3 ? 0 : (variant == 2 ? 1 + u % 255 : ((u >> 12) % 9 ? u >> 40 : 0)));
  }
}

static void reduce_init(rsimd_reduce_result *r, double f64) {
  memset(r, 0, sizeof *r);
  r->idx = -1;
  r->f64 = f64;
  r->f64_hi = -HUGE_VAL;
  r->i64 = INT64_MAX;
  r->i64_hi = INT64_MIN;
}

static void check_flags(const char *what, ptrdiff_t n, const rsimd_reduce_result *got,
                        const rsimd_reduce_result *want) {
  check_int(what, n, 0, (long) got->count, (long) want->count);
  check_int(what, n, 1, got->saw_na, want->saw_na);
  check_int(what, n, 2, got->saw_nan, want->saw_nan);
}

/* Three-valued any/all from the flags. */
static int lgl3(const rsimd_reduce_result *r, int all, int narm) {
  if (all) return r->any_false ? 0 : (r->saw_na && !narm ? 2 : 1);
  return r->any_true ? 1 : (r->saw_na && !narm ? 2 : 0);
}

/* Runs `call` (a statement using px, len, off and &r) on [0, n) whole and
   on [0, k) then [k, n). */
#define TWO_WAYS(r, init, call)                                                  \
  for (way = 0; way < 2; way++) {                                                \
    ptrdiff_t k = way ? n / 3 : n, off, len;                                     \
    init;                                                                        \
    for (off = 0; off < n; off += len) {                                         \
      len = off < k ? k - off : n - off;                                         \
      call;                                                                      \
    }                                                                            \
    REF_CHECK;                                                                   \
  }

static void test_reduce(ptrdiff_t n) {
  int variant, way, narm, check, stop, mode;
  char what[64];
  for (variant = 0; variant < 4; variant++) {
    fill_reduce(variant);
    for (narm = 0; narm < 2; narm++) {
      for (check = 0; check < 2; check++) {
        rsimd_opts o = {narm, check, RSIMD_PREC_FAST};
        int chk = check || narm;
        rsimd_reduce_result r, want;
        ptrdiff_t i;
        if (!chk && (variant == 1 || variant == 3)) continue;

        /* min/max of int32 */
        reduce_init(&want, HUGE_VAL);
        for (i = 0; i < n; i++) {
          if (chk && ri[i] == RSIMD_NA_I32) {
            want.saw_na = 1;
            continue;
          }
          want.count++;
          if (ri[i] < want.i64) want.i64 = ri[i];
          if (ri[i] > want.i64_hi) want.i64_hi = ri[i];
        }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "minmax_i32 v%d narm%d way%d", variant, narm, way); \
  check_flags(what, n, &r, &want);                                               \
  if (want.count > 0) {                                                          \
    check_int(what, n, 3, (long) r.i64, (long) want.i64);                        \
    check_int(what, n, 4, (long) r.i64_hi, (long) want.i64_hi);                  \
  }
        TWO_WAYS(r, reduce_init(&r, HUGE_VAL), RSIMD_KERNEL(minmax_i32)(ri + off, len, &r, &o))
#undef REF_CHECK

        /* prod of int32 (exact) */
        reduce_init(&want, 1.0);
        for (i = 0; i < n; i++) {
          if (chk && ri[i] == RSIMD_NA_I32) {
            want.saw_na = 1;
            if (narm) continue;
          }
          want.count++;
          want.f64 *= variant == 2 ? 1.0 : (double) (ri[i] % 3);
        }
#ifndef RSIMD_NO_F64_SIMD
        {
          static int32_t small[N + 1];
          for (i = 0; i < n; i++) small[i] = variant == 2 ? 1 : (ri[i] == RSIMD_NA_I32 ? ri[i] : ri[i] % 3);
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "prod_i32 v%d narm%d way%d", variant, narm, way);  \
  if (!want.saw_na || narm) {                                                    \
    check_int(what, n, 0, (long) r.count, (long) want.count);                    \
    n_checks++;                                                                  \
    if (bits(r.f64) != bits(want.f64)) fail(what, n, 0, "product differs");      \
  }                                                                              \
  check_int(what, n, 1, r.saw_na, want.saw_na);
          TWO_WAYS(r, reduce_init(&r, 1.0), RSIMD_KERNEL(prod_i32)(small + off, len, &r, &o))
#undef REF_CHECK
        }
#endif

        /* find_i32 */
        for (i = 0; i < 3; i++) {
          int v = i == 0 ? (n > 0 ? ri[n - 1] : 0) : (i == 1 ? 5 : 77);
          ptrdiff_t j, wantj = -1;
          for (j = 0; j < n; j++) {
            if (ri[j] == v) {
              wantj = j;
              break;
            }
          }
          check_int("find_i32", n, i, (long) RSIMD_KERNEL(find_i32)(ri, n, v), (long) wantj);
        }

        /* any/all of int32, logical and raw */
        for (stop = 0; stop < 3; stop++) {
          int all = stop == RSIMD_STOP_FALSE;
          reduce_init(&want, 0.0);
          for (i = 0; i < n; i++) {
            if (chk && ri[i] == RSIMD_NA_I32) want.saw_na = 1;
            else if (ri[i] == 0) want.any_false = 1;
            else want.any_true = 1;
          }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "anyall_i32 v%d stop%d way%d", variant, stop, way); \
  if (stop == RSIMD_STOP_NONE) {                                                 \
    check_int(what, n, 0, r.any_true, want.any_true);                            \
    check_int(what, n, 1, r.any_false, want.any_false);                          \
    check_int(what, n, 2, r.saw_na, want.saw_na);                                \
  } else {                                                                       \
    check_int(what, n, 3, lgl3(&r, all, narm), lgl3(&want, all, narm));          \
  }
          TWO_WAYS(r, reduce_init(&r, 0.0),
                   if (!(stop && (all ? r.any_false : r.any_true)))
                     RSIMD_KERNEL(anyall_i32)(ri + off, len, stop, &r, &o))
#undef REF_CHECK
          reduce_init(&want, 0.0);
          for (i = 0; i < n; i++) {
            if (rb[i] == 0) want.any_false = 1;
            else want.any_true = 1;
          }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "anyall_u8 v%d stop%d way%d", variant, stop, way); \
  if (stop == RSIMD_STOP_NONE) {                                                 \
    check_int(what, n, 0, r.any_true, want.any_true);                            \
    check_int(what, n, 1, r.any_false, want.any_false);                          \
  } else {                                                                       \
    check_int(what, n, 3, lgl3(&r, all, 0), lgl3(&want, all, 0));                \
  }
          TWO_WAYS(r, reduce_init(&r, 0.0),
                   if (!(stop && (all ? r.any_false : r.any_true)))
                     RSIMD_KERNEL(anyall_u8)(rb + off, len, stop, &r))
#undef REF_CHECK
        }

#ifndef RSIMD_NO_F64_SIMD
        /* min/max of double: exact, including the sign of zero */
        reduce_init(&want, HUGE_VAL);
        for (i = 0; i < n; i++) {
          if (isnan(rd[i])) {
            if (chk) {
              want.saw_nan = 1;
              if (is_na_ref(rd[i])) want.saw_na = 1;
            } else {
              want.count++;
            }
            continue;
          }
          want.count++;
          if (rd[i] < want.f64) want.f64 = rd[i];
          if (rd[i] > want.f64_hi) want.f64_hi = rd[i];
        }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "minmax_f64 v%d narm%d way%d", variant, narm, way); \
  check_flags(what, n, &r, &want);                                               \
  n_checks += 2;                                                                 \
  if (bits(r.f64) != bits(want.f64)) fail(what, n, 3, "minimum differs");        \
  if (bits(r.f64_hi) != bits(want.f64_hi)) fail(what, n, 4, "maximum differs");
        TWO_WAYS(r, reduce_init(&r, HUGE_VAL), RSIMD_KERNEL(minmax_f64)(rd + off, len, &r, &o))
#undef REF_CHECK

        /* prod of double: exact (powers of two), NaN matches any NaN */
        reduce_init(&want, 1.0);
        for (i = 0; i < n; i++) {
          if (chk && isnan(rd[i])) {
            want.saw_nan = 1;
            if (is_na_ref(rd[i])) want.saw_na = 1;
            if (narm) continue;
          }
          want.count++;
          want.f64 *= rd[i];
        }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "prod_f64 v%d narm%d way%d", variant, narm, way);  \
  check_flags(what, n, &r, &want);                                               \
  n_checks++;                                                                    \
  if (!(bits(r.f64) == bits(want.f64) || (isnan(r.f64) && isnan(want.f64))))     \
    fail(what, n, 3, "product differs");
        TWO_WAYS(r, reduce_init(&r, 1.0), RSIMD_KERNEL(prod_f64)(rd + off, len, &r, &o))
#undef REF_CHECK

        /* sum of deviations from 0.5 (exact in every mode) */
        for (mode = 0; mode < 3; mode++) {
          rsimd_opts om = {narm, check, mode};
          double s = 0;
          reduce_init(&want, 0.0);
          for (i = 0; i < n; i++) {
            if (chk && isnan(rd[i])) {
              want.saw_nan = 1;
              if (is_na_ref(rd[i])) want.saw_na = 1;
              if (narm) continue;
            }
            want.count++;
            s += isinf(rd[i]) ? 0.0 : rd[i] - 0.5;
          }
          if (variant != 2) continue; /* finite values only for the value */
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "sum_dev_f64 v%d mode%d way%d", variant, mode, way); \
  check_flags(what, n, &r, &want);                                               \
  n_checks++;                                                                    \
  if (rsimd_reduce_value(&r, mode) != s) fail(what, n, 3, "sum differs");
          TWO_WAYS(r, reduce_init(&r, 0.0),
                   RSIMD_KERNEL(sum_dev_f64)(rd + off, len, 0.5, &r, &om))
#undef REF_CHECK
        }

        /* find_f64: numeric equality, so 0 finds -0 */
        for (i = 0; i < 4; i++) {
          double v = i == 0 ? 0.0 : (i == 1 ? -2.0 : (i == 2 ? 7.0 : -HUGE_VAL));
          ptrdiff_t j, wantj = -1;
          for (j = 0; j < n; j++) {
            if (rd[j] == v) {
              wantj = j;
              break;
            }
          }
          check_int("find_f64", n, i, (long) RSIMD_KERNEL(find_f64)(rd, n, v), (long) wantj);
        }

        /* any/all of double */
        for (stop = 0; stop < 3; stop++) {
          int all = stop == RSIMD_STOP_FALSE;
          reduce_init(&want, 0.0);
          for (i = 0; i < n; i++) {
            if (chk && isnan(rd[i])) want.saw_na = 1;
            else if (rd[i] == 0) want.any_false = 1;
            else want.any_true = 1;
          }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "anyall_f64 v%d stop%d way%d", variant, stop, way); \
  if (stop == RSIMD_STOP_NONE) {                                                 \
    check_int(what, n, 0, r.any_true, want.any_true);                            \
    check_int(what, n, 1, r.any_false, want.any_false);                          \
    check_int(what, n, 2, r.saw_na, want.saw_na);                                \
  } else {                                                                       \
    check_int(what, n, 3, lgl3(&r, all, narm), lgl3(&want, all, narm));          \
  }
          TWO_WAYS(r, reduce_init(&r, 0.0),
                   if (!(stop && (all ? r.any_false : r.any_true)))
                     RSIMD_KERNEL(anyall_f64)(rd + off, len, stop, &r, &o))
#undef REF_CHECK
        }
#endif
      }
    }

    /* Missing values: any, count, which (as int and as double). */
    for (mode = 0; mode < 4; mode++) {
      rsimd_reduce_result r, want;
      ptrdiff_t i;
      reduce_init(&want, 0.0);
      want.i64 = 0;
      for (i = 0; i < n; i++) {
        if (ri[i] == RSIMD_NA_I32) want.i64++;
      }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "na_i32 v%d mode%d way%d", variant, mode, way);    \
  if (mode == RSIMD_NAMODE_ANY) {                                                \
    check_int(what, n, 0, r.saw_na, want.i64 > 0);                               \
  } else {                                                                       \
    check_int(what, n, 1, (long) r.i64, (long) want.i64);                        \
    for (i = 0, j = 0; mode != RSIMD_NAMODE_COUNT && i < n; i++) {               \
      if (ri[i] != RSIMD_NA_I32) continue;                                       \
      check_int(what, n, i, mode == RSIMD_NAMODE_WHICH_I32 ? iidx[j] : (long) ridx[j], \
                (long) (i + 1001));                                              \
      j++;                                                                       \
    }                                                                            \
  }
      {
        ptrdiff_t j;
        TWO_WAYS(r, (reduce_init(&r, 0.0), r.i64 = 0),
                 if (!r.saw_na) RSIMD_KERNEL(na_i32)(ri + off, len, mode, 1000 + off,
                                                    mode == RSIMD_NAMODE_WHICH_I32
                                                        ? (void *) iidx
                                                        : (void *) ridx,
                                                    &r))
      }
#undef REF_CHECK
#ifndef RSIMD_NO_F64_SIMD
      reduce_init(&want, 0.0);
      want.i64 = 0;
      for (i = 0; i < n; i++) {
        if (isnan(rd[i])) want.i64++;
      }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "na_f64 v%d mode%d way%d", variant, mode, way);    \
  if (mode == RSIMD_NAMODE_ANY) {                                                \
    check_int(what, n, 0, r.saw_na, want.i64 > 0);                               \
  } else {                                                                       \
    check_int(what, n, 1, (long) r.i64, (long) want.i64);                        \
    for (i = 0, j = 0; mode != RSIMD_NAMODE_COUNT && i < n; i++) {               \
      if (!isnan(rd[i])) continue;                                               \
      check_int(what, n, i, mode == RSIMD_NAMODE_WHICH_I32 ? iidx[j] : (long) ridx[j], \
                (long) (i + 1001));                                              \
      j++;                                                                       \
    }                                                                            \
  }
      {
        ptrdiff_t j;
        TWO_WAYS(r, (reduce_init(&r, 0.0), r.i64 = 0),
                 if (!r.saw_na) RSIMD_KERNEL(na_f64)(rd + off, len, mode, 1000 + off,
                                                    mode == RSIMD_NAMODE_WHICH_I32
                                                        ? (void *) iidx
                                                        : (void *) ridx,
                                                    &r))
      }
#undef REF_CHECK
#endif
    }
  }
}

/* ---- Fused reductions and scans ----------------------------------------- */

static double ld[N + 1], ld2[N + 1];
static int32_t li[N + 1], li2[N + 1];

/* Small integer-valued data, so that every sum, product and square below
   is exact in every precision mode and with or without fused
   multiply-adds; variant 1 adds NaN, NA and signed zeros, variant 2 large
   integers (int32 cumsum overflow). */
static void fill_linalg(int variant) {
  const double nan = from_bits(UINT64_C(0x7FF8000000000000));
  const double na = from_bits(UINT64_C(0x7FF00000000007A2));
  int i;
  for (i = 0; i < N; i++) {
    uint64_t u = next_rand();
    int32_t a = (int32_t) (u % 41) - 20, b = (int32_t) ((u >> 16) % 41) - 20;
    if (variant == 2) a = (int32_t) ((u >> 8) % 2 ? INT32_MAX / 3 : -(INT32_MAX / 3)) + a;
    li[i] = a;
    li2[i] = b;
    ld[i] = a == 0 && (u >> 40) % 2 ? -0.0 : (double) a;
    ld2[i] = (double) b;
    if (variant == 1) {
      if ((u >> 44) % 23 == 0) ld[i] = (u >> 50) % 2 ? nan : na;
      if ((u >> 52) % 29 == 0) ld2[i] = nan;
      if ((u >> 56) % 31 == 0) li[i] = RSIMD_NA_I32;
      if ((u >> 58) % 37 == 0) li2[i] = RSIMD_NA_I32;
    }
  }
}

#if !defined(RSIMD_SKIP_sumsq_f64) || !defined(RSIMD_SKIP_cumsum_f64)
/* Element i of an operand as the kernels read it. */
static double elt(const void *p, int i32, ptrdiff_t i) {
  if (i32) {
    int32_t v = ((const int32_t *) p)[i];
    return v == RSIMD_NA_I32 ? from_bits(UINT64_C(0x7FF00000000007A2)) : (double) v;
  }
  return ((const double *) p)[i];
}
#endif

#if !defined(RSIMD_SKIP_sumsq_f64)
/* The reference of a fold of term (sq 0, abs 1, xy 2, (x-y)^2 3, (x-c)^2 4)
   over x and y. */
static void ref_fold(const void *x, int xi, const void *y, int yi, ptrdiff_t n, int term,
                     double c, int check, int narm, rsimd_reduce_result *want) {
  ptrdiff_t i;
  reduce_init(want, 0.0);
  for (i = 0; i < n; i++) {
    double a = elt(x, xi, i), b = y ? elt(y, yi, i) : 0.0, t;
    int pair = term == 2 || term == 3;
    if (check && (isnan(a) || (pair && isnan(b)))) {
      want->saw_nan = 1;
      if (is_na_ref(a) || (pair && is_na_ref(b))) want->saw_na = 1;
      if (narm) continue;
    }
    want->count++;
    t = term == 0 ? a * a : term == 1 ? fabs(a) : term == 2 ? a * b
        : term == 3 ? (a - b) * (a - b) : (a - c) * (a - c);
    want->f64 += t;
  }
}

static void check_fold(const char *what, ptrdiff_t n, const rsimd_reduce_result *r, int mode,
                       const rsimd_reduce_result *want, int narm) {
  double v = rsimd_reduce_value(r, mode);
  check_flags(what, n, r, want);
  n_checks++;
  if ((!want->saw_nan || narm) && v != want->f64) fail(what, n, 3, "sum differs");
}
#endif

static void test_linalg(ptrdiff_t n) {
  int variant, narm, check, mode, way;
  char what[80];
  for (variant = 0; variant < 2; variant++) {
    fill_linalg(variant);
    for (narm = 0; narm < 2; narm++) {
      for (check = 0; check < 2; check++) {
        int chk = check || narm;
        if (!chk && variant == 1) continue;
        for (mode = 0; mode < 3; mode++) {
          rsimd_opts o = {narm, check, mode};
          rsimd_reduce_result r;
          int64_t isum = 0;
          ptrdiff_t i;
          int ina = 0;
          /* sumabs_i32: exact int64 */
          for (i = 0; i < n; i++) {
            if (chk && li[i] == RSIMD_NA_I32) {
              ina = 1;
              if (narm) continue;
            }
            isum += li[i] < 0 ? -(int64_t) li[i] : li[i];
          }
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "sumabs_i32 v%d narm%d chk%d way%d", variant, narm, check, way); \
  check_int(what, n, 1, r.saw_na, ina);                                          \
  if (!ina || narm) check_int(what, n, 2, (long) r.i64, (long) isum);
          TWO_WAYS(r, reduce_init(&r, 0.0); r.i64 = 0,
                   RSIMD_KERNEL(sumabs_i32)(li + off, len, &r, &o))
#undef REF_CHECK
#if !defined(RSIMD_SKIP_sumsq_f64)
          rsimd_reduce_result want;
          int t;
          ref_fold(ld, 0, NULL, 0, n, 0, 0.0, chk, narm, &want);
#define REF_CHECK                                                                \
  snprintf(what, sizeof what, "%s v%d narm%d chk%d mode%d way%d", name, variant, narm, check, \
           mode, way);                                                           \
  check_fold(what, n, &r, mode, &want, narm);
          {
            const char *name = "sumsq_f64";
            TWO_WAYS(r, reduce_init(&r, 0.0), RSIMD_KERNEL(sumsq_f64)(ld + off, len, &r, &o))
            name = "sumsq_i32";
            ref_fold(li, 1, NULL, 0, n, 0, 0.0, chk, narm, &want);
            TWO_WAYS(r, reduce_init(&r, 0.0), RSIMD_KERNEL(sumsq_i32)(li + off, len, &r, &o))
            name = "sumabs_f64";
            ref_fold(ld, 0, NULL, 0, n, 1, 0.0, chk, narm, &want);
            TWO_WAYS(r, reduce_init(&r, 0.0), RSIMD_KERNEL(sumabs_f64)(ld + off, len, &r, &o))
            name = "var_pass2_f64";
            ref_fold(ld, 0, NULL, 0, n, 4, 0.5, chk, narm, &want);
            TWO_WAYS(r, reduce_init(&r, 0.0),
                     RSIMD_KERNEL(var_pass2_f64)(ld + off, len, 0.5, &r, &o))
            name = "var_pass2_i32";
            ref_fold(li, 1, NULL, 0, n, 4, -1.5, chk, narm, &want);
            TWO_WAYS(r, reduce_init(&r, 0.0),
                     RSIMD_KERNEL(var_pass2_i32)(li + off, len, -1.5, &r, &o))
            for (t = 0; t < 3; t++) {
              const void *x = t == RSIMD_PAIR_I32_I32 ? (const void *) li : (const void *) ld;
              const void *y = t == RSIMD_PAIR_F64_F64 ? (const void *) ld2 : (const void *) li2;
              int xi = t == RSIMD_PAIR_I32_I32, yi = t != RSIMD_PAIR_F64_F64;
              size_t xs = xi ? sizeof(int32_t) : sizeof(double);
              size_t ys = yi ? sizeof(int32_t) : sizeof(double);
              name = t == 0 ? "dot ff" : t == 1 ? "dot ii" : "dot fi";
              ref_fold(x, xi, y, yi, n, 2, 0.0, chk, narm, &want);
              TWO_WAYS(r, reduce_init(&r, 0.0),
                       RSIMD_KERNEL(dot_f64)((const char *) x + off * xs,
                                             (const char *) y + off * ys, len, t, &r, &o))
              name = t == 0 ? "dist ff" : t == 1 ? "dist ii" : "dist fi";
              ref_fold(x, xi, y, yi, n, 3, 0.0, chk, narm, &want);
              TWO_WAYS(r, reduce_init(&r, 0.0),
                       RSIMD_KERNEL(dist_f64)((const char *) x + off * xs,
                                              (const char *) y + off * ys, len, t, &r, &o))
              if (!narm) {
                rsimd_reduce_result rc[3], w3[3];
                int k;
                ref_fold(x, xi, y, yi, n, 2, 0.0, chk, 0, &w3[0]);
                ref_fold(x, xi, NULL, 0, n, 0, 0.0, chk, 0, &w3[1]);
                ref_fold(y, yi, NULL, 0, n, 0, 0.0, chk, 0, &w3[2]);
                for (way = 0; way < 2; way++) {
                  ptrdiff_t kk = way ? n / 3 : n, off, len;
                  for (k = 0; k < 3; k++) reduce_init(&rc[k], 0.0);
                  for (off = 0; off < n; off += len) {
                    len = off < kk ? kk - off : n - off;
                    RSIMD_KERNEL(cosine_f64)((const char *) x + off * xs,
                                             (const char *) y + off * ys, len, t, rc, &o);
                  }
                  for (k = 0; k < 3; k++) {
                    snprintf(what, sizeof what, "cosine t%d sum%d v%d chk%d mode%d way%d", t, k,
                             variant, check, mode, way);
                    check_fold(what, n, &rc[k], mode, &w3[k], 0);
                  }
                }
              }
            }
          }
#undef REF_CHECK
#endif
        }
      }
    }
  }
}

/* Scans: exact on small integer-valued data, including the stop index and
   the sign of zero. */
#if !defined(RSIMD_SKIP_cumsum_f64)
static double lout_f[N + 1], lref_f[N + 1];
static int32_t lout_i[N + 1], lref_i[N + 1];

static void ref_scan(const void *x, int in_i32, ptrdiff_t n, int op, int out_i32, double init,
                     ptrdiff_t *stop, int *ovf) {
  double acc = init;
  ptrdiff_t i;
  *stop = -1;
  *ovf = 0;
  for (i = 0; i < n; i++) {
    double v;
    if (in_i32 ? ((const int32_t *) x)[i] == RSIMD_NA_I32 : isnan(((const double *) x)[i])) {
      *stop = i;
      return;
    }
    v = elt(x, in_i32, i);
    acc = op == 0   ? acc + v
          : op == 1 ? acc * v
          : op == 2 ? (acc < v ? acc : v)
                    : (acc > v ? acc : v);
    if (out_i32) {
      if (op == 0 && (acc > INT32_MAX || acc < -INT32_MAX)) {
        *stop = i;
        *ovf = 1;
        return;
      }
      lref_i[i] = (int32_t) acc;
    } else {
      lref_f[i] = acc;
    }
  }
}

static void test_scan(ptrdiff_t n) {
  static const double inits[4] = {0.0, 1.0, HUGE_VAL, -HUGE_VAL};
  int variant, op, in_i32, way;
  char what[80];
  for (variant = 0; variant < 3; variant++) {
    fill_linalg(variant);
    for (op = 0; op < 4; op++) {
      for (in_i32 = 0; in_i32 < 2; in_i32++) {
        int out_i32 = in_i32 && op != 1, ovf;
        ptrdiff_t stop, i;
        const void *x = in_i32 ? (const void *) li : (const void *) ld;
        size_t xs = in_i32 ? sizeof(int32_t) : sizeof(double);
        if (op == 1) {
          /* products of +-1, +-2 and +-0.5 stay exact */
          static double pd[N + 1];
          static int32_t pi[N + 1];
          for (i = 0; i < n; i++) {
            int32_t a = li[i] == RSIMD_NA_I32 ? 0 : li[i];
            pi[i] = li[i] == RSIMD_NA_I32 ? li[i] : (a % 2 ? -1 : 1) * (a % 3 ? 1 : 2);
            pd[i] = isnan(ld[i]) ? ld[i] : (a % 5 == 0 ? 0.5 : (double) pi[i]);
          }
          x = in_i32 ? (const void *) pi : (const void *) pd;
        }
        if (in_i32 && variant == 2 && op != 0) continue;
        ref_scan(x, in_i32, n, op, out_i32, inits[op], &stop, &ovf);
        for (way = 0; way < 2; way++) {
          ptrdiff_t k = way ? n / 3 : n, off, len, got = -1;
          rsimd_scan_state s = {inits[op], 0.0, 0};
          rsimd_opts o = {0, 1, RSIMD_PREC_FAST};
          for (i = 0; i <= n; i++) {
            lout_f[i] = SENTINEL_F64;
            lout_i[i] = SENTINEL_I32;
          }
          for (off = 0; off < n && got < 0; off += len) {
            const void *px = (const char *) x + off * xs;
            ptrdiff_t r;
            len = off < k ? k - off : n - off;
            if (in_i32) {
              r = op == 0 ? RSIMD_KERNEL(cumsum_i32)((const int *) px, len, lout_i + off, &s)
                  : op == 1 ? RSIMD_KERNEL(cumprod_i32)((const int *) px, len, lout_f + off, &s)
                            : RSIMD_KERNEL(cumminmax_i32)((const int *) px, len, op == 3,
                                                          lout_i + off, &s);
            } else {
              r = op == 0 ? RSIMD_KERNEL(cumsum_f64)((const double *) px, len, lout_f + off, &s, &o)
                  : op == 1 ? RSIMD_KERNEL(cumprod_f64)((const double *) px, len, lout_f + off, &s)
                            : RSIMD_KERNEL(cumminmax_f64)((const double *) px, len, op == 3,
                                                          lout_f + off, &s);
            }
            if (r >= 0) got = off + r;
          }
          snprintf(what, sizeof what, "scan op%d i32 %d v%d way%d", op, in_i32, variant, way);
          check_int(what, n, 0, (long) got, (long) stop);
          check_int(what, n, 1, s.overflow, ovf);
          if (out_i32) check_i32(what, stop < 0 ? n : stop, lout_i, lref_i);
          else check_f64(what, stop < 0 ? n : stop, lout_f, lref_f, 0);
        }
      }
    }
    /* compensated cumsum: the scalar Neumaier sum on every tier */
    {
      rsimd_scan_state s = {0.0, 0.0, 0};
      rsimd_opts o = {0, 1, RSIMD_PREC_COMPENSATED};
      double S = 0.0, C = 0.0;
      ptrdiff_t i, stop = -1;
      for (i = 0; i < n; i++) {
        if (isnan(ld[i])) {
          stop = i;
          break;
        }
        rsimd_neumaier_add(&S, &C, ld[i] * 1e15 + 0.25);
        lref_f[i] = rsimd_neumaier_value(S, C);
        fc[i] = ld[i] * 1e15 + 0.25;
      }
      for (; i < n; i++) fc[i] = ld[i];
      check_int("cumsum comp stop", n, 0, (long) RSIMD_KERNEL(cumsum_f64)(fc, n, lout_f, &s, &o),
                (long) stop);
      check_f64("cumsum comp", stop < 0 ? n : stop, lout_f, lref_f, 0);
    }
  }
}
#else
static void test_scan(ptrdiff_t n) { (void) n; }
#endif

#undef TWO_WAYS

int main(void) {
  ptrdiff_t n;
  init_inputs();
  for (n = 0; n <= N; n++) {
#ifndef RSIMD_NO_F64_SIMD
    test_f64(n);
#endif
    test_int(n);
    test_convert(n);
    test_predicates(n);
    test_na_int(n);
#ifndef RSIMD_NO_F64_SIMD
    test_na_fold(n);
#endif
    test_reduce(n);
    test_linalg(n);
    test_scan(n);
  }
#ifndef RSIMD_NO_F64_SIMD
  test_na_cancel();
#endif
#ifndef RSIMD_NO_F64_SIMD
  test_f64_horizontal();
#endif
  test_int_horizontal();
  printf("tier %s: lanes64=%ld lanes32=%ld, %ld checks, %ld failures\n", RSIMD_TIER_STRING,
         (long) RSIMD_LANES_64, (long) RSIMD_LANES_32, n_checks, n_fail);
  return n_fail != 0;
}
