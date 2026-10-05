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
 * chunks, and compared with plain C references too, and so are the
 * elementwise kernels (src/kernels/arith.inc.c) with every combination of
 * broadcast and int32 operands, and so are the predicate, comparison,
 * logical, bitwise and conversion kernels (src/kernels/predicates.inc.c,
 * compare.inc.c, bitwise.inc.c, convert.inc.c), against the scalar forms
 * those files define for the none tier, and so are the complex kernels
 * (src/kernels/complex.inc.c: multiply and divide in every rounding
 * variant, products, comparisons, scans; and cmath.inc.c's Mod and Arg)
 * and the Hamming kernels (compare.inc.c,
 * bitwise.inc.c, int64.inc.c, complex.inc.c, against plain loops). On tiers built with SLEEF the
 * elementary-function wrappers (accurate and fast) are compared with long
 * double libm, and the elementary-function kernels (src/kernels/math.inc.c,
 * in both accuracy modes) with the scalar forms that file defines for the
 * none tier, and the softmax passes
 * (src/kernels/ml.inc.c) with libm and the scalar folds. Prints a
 * summary and exits non-zero on any mismatch.
 */

#include <float.h>
#include <stdio.h>
/* Flush the lane counts of count_na_i32_() every two iterations, and
   those of the Hamming kernels every two vectors. */
#define RSIMD_COUNT_BLOCK 2
#define RSIMD_HAMMING_BLOCK 2
/* Model the float estimates of sse2 and avx2 (RCPPS, RSQRTPS) with the
   largest relative error they are allowed, 1.5 * 2^-12, so that the
   approximation bound is checked against the worst case. */
#define RSIMD_EST_SKEW 0x1.8p-12
#include "kernels/common.inc.h"
#include "kernels/reduce.inc.c"
#include "kernels/scan.inc.c"
#include "kernels/arith.inc.c"
#include "kernels/predicates.inc.c"
#include "kernels/complex.inc.c"
#include "kernels/compare.inc.c"
#include "kernels/bitwise.inc.c"
#include "kernels/int64.inc.c"
#include "kernels/convert.inc.c"
#include "kernels/math.inc.c"
#include "kernels/cmath.inc.c"
#include "kernels/ml.inc.c"

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
  F64_OP("f64 ceil", rsimd_vf64_ceil(x), ceil(a), 0);
  F64_OP("f64 trunc", rsimd_vf64_trunc(x), trunc(a), 0);
  F64_OP("f64 rint", rsimd_vf64_rint(x), nearbyint(a), 0);
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
  {
    int k;
    for (k = 0; k < 32; k++) {
      I32_OP("i32 sll", rsimd_vi32_sll(x, k), (int32_t) ((uint32_t) a << k));
      I32_OP("i32 srl", rsimd_vi32_srl(x, k), (int32_t) ((uint32_t) a >> k));
      I32_OP("i32 sra", rsimd_vi32_sra(x, k),
             a < 0 ? (int32_t) ~(~(uint32_t) a >> k) : (int32_t) ((uint32_t) a >> k));
    }
  }
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
  {
    int k;
    for (k = 0; k < 64; k += k < 62 ? 9 : 1) {
      I64_OP("i64 sll", rsimd_vi64_sll(x, k), (int64_t) ((uint64_t) a << k));
      I64_OP("i64 srl", rsimd_vi64_srl(x, k), (int64_t) ((uint64_t) a >> k));
      I64_OP("i64 sra", rsimd_vi64_sra(x, k),
             a < 0 ? (int64_t) ~(~(uint64_t) a >> k) : (int64_t) ((uint64_t) a >> k));
    }
  }
  I64_OP("i64 sign", rsimd_vi64_sign(x), a < 0 ? -1 : 0);
  I64_OP("i64 mulu32", rsimd_vi64_mulu32(x, y),
         (int64_t) (((uint64_t) a & 0xFFFFFFFFu) * ((uint64_t) b & 0xFFFFFFFFu)));
  I64_OP("i64 mul32", rsimd_vi64_mul32(x, y),
         (int64_t) (int32_t) (uint32_t) a * (int64_t) (int32_t) (uint32_t) b);
#ifndef RSIMD_NO_F64_SIMD
  reset_out();
  LOOP(RSIMD_LANES_64, p64, vf64, rsimd_vi64 x = rsimd_vi64_loadu_p(pg, la + i, 7),
       rsimd_vi64_to_vf64(x), fout);
  for (j = 0; j < n; j++) fref[j] = (double) la[j];
  check_f64("i64 to_vf64", n, fout, fref, 1);
#endif
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

#ifdef RSIMD_HAVE_SLEEF
/* The SLEEF wrappers of common.inc.h, through the predicated loop, against
   SLEEF's documented error bounds: 1.0 ULP for the _u10 functions, 0.5
   ULP for the _u05 ones and 3.5 ULP for the _u35 ones (the _fast
   wrappers that have one), measured against the long double libm function
   (quad precision on aarch64, x87 extended on x86), with a little slack for
   the reference's own error. Results whose reference rounds to NaN, an
   infinity or zero must match exactly, sign included. Inputs are random
   over each function's domain of accuracy, plus special values. SLEEF's
   bounds only hold for sinh and cosh on [-709, 709] and asinh and acosh on
   [-1.34e154, 1.34e154] (beyond, they may give infinities); sinpi, cospi
   and sincospi return 0 and 1 beyond 2.5e8 (SLEEF documents 1e9) and give
   sinpi(n) = -0 for odd n > 0, so for them a zero matches either zero. */
static long double ldref[N + 1];

static double ulp_of(double x) {
  int e;
  if (x == 0 || !isfinite(x)) return DBL_TRUE_MIN;
  frexp(x, &e);
  return ldexp(1.0, e - 53) > DBL_TRUE_MIN ? ldexp(1.0, e - 53) : DBL_TRUE_MIN;
}

static int sleef_any_zero = 0; /* 1: a zero result matches either zero */
/* The worst errors seen at bounds 1, 0.5 and 3.5. */
static double sleef_worst[3] = {0, 0, 0};
static const char *sleef_worst_name[3] = {"", "", ""};

static void check_ulp(const char *what, ptrdiff_t n, const double *got, double bound) {
  ptrdiff_t i;
  char buf[200];
  const int k = bound < 1 ? 1 : bound > 1 ? 2 : 0;
  for (i = 0; i < n; i++) {
    double g = got[i], w = (double) ldref[i];
    int ok;
    n_checks++;
    if (isnan(w) || isinf(w) || w == 0) {
      ok = bits(g) == bits(w) || (isnan(g) && isnan(w)) || (sleef_any_zero && g == 0 && w == 0);
    } else {
      double err = (double) (fabsl((long double) g - ldref[i]) / ulp_of(w));
      /* SLEEF's bounds are max(bound, DBL_MIN) for atan2, sinpi and cospi:
         subnormal results may be off by up to DBL_MIN. */
      ok = err <= bound + 0.01 ||
           (fabs(w) < DBL_MIN && fabsl((long double) g - ldref[i]) <= DBL_MIN);
      if (ok && fabs(w) >= DBL_MIN && err > sleef_worst[k]) {
        sleef_worst[k] = err;
        sleef_worst_name[k] = what;
      }
    }
    if (!ok) {
      snprintf(buf, sizeof buf, "input %a %a: got %a want %a (%.4Lg)", fa[i], fb[i], g, w,
               ldref[i]);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (bits(got[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}

/* Fills v with n values: specials first, then uniform on [lo, hi] or,
   when log_scale, signs times 10^u for u uniform on [lo, hi]. */
static void sleef_inputs(double *v, ptrdiff_t n, double lo, double hi, int log_scale, int sign) {
  static const double specials[] = {0.0, -0.0, 1.0, -1.0, 0.5, -0.5, INFINITY, -INFINITY, NAN};
  ptrdiff_t i;
  for (i = 0; i < n; i++) {
    double u = (double) (next_rand() >> 11) * 0x1.0p-53;
    if (i < (ptrdiff_t) (sizeof specials / sizeof specials[0])) {
      v[i] = specials[i];
    } else if (log_scale) {
      v[i] = pow(10.0, lo + u * (hi - lo));
      if (sign && (next_rand() & 1)) v[i] = -v[i];
    } else {
      v[i] = lo + u * (hi - lo);
    }
  }
}

static rsimd_vf64 sleef_sincos_s(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincos(x, &s, &c);
  return s;
}
static rsimd_vf64 sleef_sincos_c(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincos(x, &s, &c);
  return c;
}
static rsimd_vf64 sleef_sincospi_s(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincospi(x, &s, &c);
  return s;
}
static rsimd_vf64 sleef_sincospi_c(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincospi(x, &s, &c);
  return c;
}
static rsimd_vf64 sleef_sincos_fast_s(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincos_fast(x, &s, &c);
  return s;
}
static rsimd_vf64 sleef_sincos_fast_c(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincos_fast(x, &s, &c);
  return c;
}
static rsimd_vf64 sleef_sincospi_fast_s(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincospi_fast(x, &s, &c);
  return s;
}
static rsimd_vf64 sleef_sincospi_fast_c(rsimd_vf64 x) {
  rsimd_vf64 s, c;
  rsimd_sleef_sincospi_fast(x, &s, &c);
  return c;
}

static long double ref_exp10(double x) { return powl(10.0L, (long double) x); }
/* i686 glibc's acoshl(-0) is -Inf. */
static long double ref_acosh(double x) {
  return x < 1 ? (long double) NAN : acoshl((long double) x);
}

/* sin(pi x) and cos(pi x) from an exact reduction to pi t, |t| <= 1/4,
   so that x87 long double is accurate enough near the zeros; exact at
   multiples of 1/2 (sinpi keeps the sign of a zero argument, as C23 does). */
/* sin(pi (r + shift / 2)) with q the nearest integer to 2 r. */
static long double sinpi_quadrant(long double r, int q, int shift) {
  const long double pi = 3.141592653589793238462643383279502884L;
  long double t = r - 0.5L * (long double) q; /* exact, |t| <= 1/4 */
  switch ((q + shift) & 3) {
  case 0: return sinl(pi * t);
  case 1: return cosl(pi * t);
  case 2: return -sinl(pi * t);
  default: return -cosl(pi * t);
  }
}
static long double ref_sinpi(double x) {
  long double r;
  if (!isfinite(x)) return NAN;
  r = fmodl((long double) x, 2.0L); /* exact */
  if (r == 0 || fabsl(r) == 1) return copysignl(0.0L, (long double) x);
  return sinpi_quadrant(r, (int) nearbyintl(2 * r), 0);
}
static long double ref_cospi(double x) {
  long double r;
  if (!isfinite(x)) return NAN;
  r = fmodl(fabsl((long double) x), 2.0L);
  if (r == 0.5L || r == 1.5L) return 0.0L;
  return sinpi_quadrant(r, (int) nearbyintl(2 * r), 1);
}

#define SLEEF1(name, FN, REF, bound, lo, hi, log_scale, sign)                  \
  do {                                                                         \
    sleef_inputs(fa, n, lo, hi, log_scale, sign);                              \
    reset_out();                                                               \
    LOOP(RSIMD_LANES_64, p64, vf64, F64_LOADS, FN(x), fout);                   \
    for (j = 0; j < n; j++) ldref[j] = REF(fa[j]);                             \
    check_ulp(name, n, fout, bound);                                           \
  } while (0)
#define SLEEF2(name, FN, REF, bound, lo, hi, log_scale, sign)                  \
  do {                                                                         \
    sleef_inputs(fa, n, lo, hi, log_scale, sign);                              \
    sleef_inputs(fb, n, lo, hi, log_scale, sign);                              \
    reset_out();                                                               \
    LOOP(RSIMD_LANES_64, p64, vf64, F64_LOADS, FN(x, y), fout);                \
    for (j = 0; j < n; j++) ldref[j] = REF(fa[j], fb[j]);                      \
    check_ulp(name, n, fout, bound);                                           \
  } while (0)

static void test_sleef(void) {
  const ptrdiff_t n = N;
  ptrdiff_t j;
  int rep;
  for (rep = 0; rep < 20; rep++) {
    SLEEF1("sleef exp", rsimd_sleef_exp, expl, 1, -745, 709.7, 0, 0);
    SLEEF1("sleef exp2", rsimd_sleef_exp2, exp2l, 1, -1075, 1023.9, 0, 0);
    SLEEF1("sleef exp10", rsimd_sleef_exp10, ref_exp10, 1, -323, 308.2, 0, 0);
    SLEEF1("sleef expm1", rsimd_sleef_expm1, expm1l, 1, -40, 709.7, 0, 0);
    SLEEF1("sleef expm1 small", rsimd_sleef_expm1, expm1l, 1, -12, -1, 1, 1);
    SLEEF1("sleef log", rsimd_sleef_log, logl, 1, -320, 308, 1, 0);
    SLEEF1("sleef log negative", rsimd_sleef_log, logl, 1, -10, 10, 1, 1);
    SLEEF1("sleef log2", rsimd_sleef_log2, log2l, 1, -320, 308, 1, 0);
    SLEEF1("sleef log10", rsimd_sleef_log10, log10l, 1, -320, 308, 1, 0);
    SLEEF1("sleef log1p", rsimd_sleef_log1p, log1pl, 1, -0.999, 1e3, 0, 0);
    SLEEF1("sleef log1p small", rsimd_sleef_log1p, log1pl, 1, -300, -1, 1, 1);
    SLEEF1("sleef cbrt", rsimd_sleef_cbrt, cbrtl, 1, -320, 308, 1, 1);
    SLEEF1("sleef sin", rsimd_sleef_sin, sinl, 1, -1e4, 1e4, 0, 0);
    SLEEF1("sleef sin wide", rsimd_sleef_sin, sinl, 1, -300, 300, 1, 1);
    SLEEF1("sleef cos", rsimd_sleef_cos, cosl, 1, -1e4, 1e4, 0, 0);
    SLEEF1("sleef cos wide", rsimd_sleef_cos, cosl, 1, -300, 300, 1, 1);
    SLEEF1("sleef tan", rsimd_sleef_tan, tanl, 1, -1e4, 1e4, 0, 0);
    SLEEF1("sleef tan wide", rsimd_sleef_tan, tanl, 1, -300, 300, 1, 1);
    SLEEF1("sleef sincos sin", sleef_sincos_s, sinl, 1, -300, 300, 1, 1);
    SLEEF1("sleef sincos cos", sleef_sincos_c, cosl, 1, -300, 300, 1, 1);
    SLEEF1("sleef asin", rsimd_sleef_asin, asinl, 1, -1.01, 1.01, 0, 0);
    SLEEF1("sleef acos", rsimd_sleef_acos, acosl, 1, -1.01, 1.01, 0, 0);
    SLEEF1("sleef atan", rsimd_sleef_atan, atanl, 1, -300, 300, 1, 1);
    SLEEF1("sleef sinh", rsimd_sleef_sinh, sinhl, 1, -709, 709, 0, 0);
    SLEEF1("sleef cosh", rsimd_sleef_cosh, coshl, 1, -709, 709, 0, 0);
    SLEEF1("sleef tanh", rsimd_sleef_tanh, tanhl, 1, -20, 20, 0, 0);
    SLEEF1("sleef asinh", rsimd_sleef_asinh, asinhl, 1, -300, 154.12, 1, 1);
    SLEEF1("sleef acosh", rsimd_sleef_acosh, ref_acosh, 1, 0, 154.12, 1, 0);
    SLEEF1("sleef atanh", rsimd_sleef_atanh, atanhl, 1, -1.01, 1.01, 0, 0);
    sleef_any_zero = 1;
    SLEEF1("sleef sinpi", rsimd_sleef_sinpi, ref_sinpi, 0.5, -2.5e8, 2.5e8, 0, 0);
    SLEEF1("sleef sinpi small", rsimd_sleef_sinpi, ref_sinpi, 0.5, -4, 4, 0, 0);
    SLEEF1("sleef cospi", rsimd_sleef_cospi, ref_cospi, 0.5, -2.5e8, 2.5e8, 0, 0);
    SLEEF1("sleef cospi small", rsimd_sleef_cospi, ref_cospi, 0.5, -4, 4, 0, 0);
    SLEEF1("sleef sincospi sin", sleef_sincospi_s, ref_sinpi, 0.5, -100, 100, 0, 0);
    SLEEF1("sleef sincospi cos", sleef_sincospi_c, ref_cospi, 0.5, -100, 100, 0, 0);
    sleef_any_zero = 0;
    SLEEF2("sleef pow", rsimd_sleef_pow, powl, 1, -3, 3, 1, 1);
    SLEEF2("sleef atan2", rsimd_sleef_atan2, atan2l, 1, -300, 300, 1, 1);
    SLEEF2("sleef hypot", rsimd_sleef_hypot, hypotl, 0.5, -300, 300, 1, 1);
    /* Fast mode: the _u35 functions over the same domains. */
    SLEEF1("fast log", rsimd_sleef_log_fast, logl, 3.5, -320, 308, 1, 0);
    SLEEF1("fast log negative", rsimd_sleef_log_fast, logl, 3.5, -10, 10, 1, 1);
    SLEEF1("fast log2", rsimd_sleef_log2_fast, log2l, 3.5, -320, 308, 1, 0);
    SLEEF1("fast cbrt", rsimd_sleef_cbrt_fast, cbrtl, 3.5, -320, 308, 1, 1);
    SLEEF1("fast sin", rsimd_sleef_sin_fast, sinl, 3.5, -1e4, 1e4, 0, 0);
    SLEEF1("fast sin wide", rsimd_sleef_sin_fast, sinl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast cos", rsimd_sleef_cos_fast, cosl, 3.5, -1e4, 1e4, 0, 0);
    SLEEF1("fast cos wide", rsimd_sleef_cos_fast, cosl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast tan", rsimd_sleef_tan_fast, tanl, 3.5, -1e4, 1e4, 0, 0);
    SLEEF1("fast tan wide", rsimd_sleef_tan_fast, tanl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast sincos sin", sleef_sincos_fast_s, sinl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast sincos cos", sleef_sincos_fast_c, cosl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast asin", rsimd_sleef_asin_fast, asinl, 3.5, -1.01, 1.01, 0, 0);
    SLEEF1("fast acos", rsimd_sleef_acos_fast, acosl, 3.5, -1.01, 1.01, 0, 0);
    SLEEF1("fast atan", rsimd_sleef_atan_fast, atanl, 3.5, -300, 300, 1, 1);
    SLEEF1("fast sinh", rsimd_sleef_sinh_fast, sinhl, 3.5, -709, 709, 0, 0);
    SLEEF1("fast cosh", rsimd_sleef_cosh_fast, coshl, 3.5, -709, 709, 0, 0);
    SLEEF1("fast tanh", rsimd_sleef_tanh_fast, tanhl, 3.5, -20, 20, 0, 0);
    sleef_any_zero = 1;
    SLEEF1("fast sinpi", rsimd_sleef_sinpi_fast, ref_sinpi, 3.5, -2.5e8, 2.5e8, 0, 0);
    SLEEF1("fast sinpi small", rsimd_sleef_sinpi_fast, ref_sinpi, 3.5, -4, 4, 0, 0);
    SLEEF1("fast cospi", rsimd_sleef_cospi_fast, ref_cospi, 3.5, -2.5e8, 2.5e8, 0, 0);
    SLEEF1("fast cospi small", rsimd_sleef_cospi_fast, ref_cospi, 3.5, -4, 4, 0, 0);
    SLEEF1("fast sincospi sin", sleef_sincospi_fast_s, ref_sinpi, 3.5, -100, 100, 0, 0);
    SLEEF1("fast sincospi cos", sleef_sincospi_fast_c, ref_cospi, 3.5, -100, 100, 0, 0);
    sleef_any_zero = 0;
    SLEEF2("fast atan2", rsimd_sleef_atan2_fast, atan2l, 3.5, -300, 300, 1, 1);
    SLEEF2("fast hypot", rsimd_sleef_hypot_fast, hypotl, 3.5, -300, 300, 1, 1);
  }
}

/* The none tier's scalar sinpi, cospi and tanpi (the reference for the
   kernels) against long double, over random arguments in [-4, 4] and
   multiples of 2^-15: within 1.1 ULP. */
static long double ref_tanpi(double x) {
  const long double pi = 3.141592653589793238462643383279502884L;
  long double r = fmodl((long double) x, 1.0L), a;
  if (r <= -0.5L) r += 1;
  if (r > 0.5L) r -= 1;
  a = fabsl(r);
  a = a <= 0.25L ? tanl(pi * a) : 1 / tanl(pi * (0.5L - a));
  return r < 0 ? -a : a;
}
static double pi_oracle_worst[3] = {0, 0, 0};
static void test_pi_oracle(void) {
  static const char *const names[] = {"none sinpi", "none cospi", "none tanpi"};
  const double bounds[] = {1.1, 1.1, 1.1};
  long i;
  int k;
  char buf[160];
  for (i = 0; i < 400000; i++) {
    double x = i % 2 ? ((double) (next_rand() >> 11) * 0x1.0p-53 - 0.5) * 8
                     : (double) (i / 2 - 100000) * 0x1.0p-15;
    for (k = 0; k < 3; k++) {
      double g = k == 0 ? rsimd_sinpi_f64(x) : k == 1 ? rsimd_cospi_f64(x) : rsimd_tanpi_f64(x);
      long double w = k == 0 ? ref_sinpi(x) : k == 1 ? ref_cospi(x) : ref_tanpi(x);
      double err;
      n_checks++;
      if (w == 0 || isnan(g)) {
        if (g != 0 && !(k == 2 && isnan(g) && fabsl(w) > 1e15L)) fail(names[k], 0, i, "zero");
        continue;
      }
      err = (double) (fabsl((long double) g - w) / ulp_of((double) w));
      if (err > pi_oracle_worst[k]) pi_oracle_worst[k] = err;
      if (err > bounds[k] + 0.001) {
        snprintf(buf, sizeof buf, "input %a: got %a want %.21Lg (%.3f ULP)", x, g, w, err);
        fail(names[k], 0, i, buf);
      }
    }
  }
}

/* The elementary-function kernels against the scalar forms of math.inc.c
   (libm, as the none tier): within 2 ULP (4 for cbrt), or 2.5 ULP more in
   fast mode (RSIMD_MATH_FAST), with NaN kinds (NA or NaN),
   infinities, signed zeros and the status bit exactly equal; pow exactly
   equal at base R's special cases. Each op runs over every short length
   and N, with double and int32 operands and every broadcast. */
static int math_ref_status;
static double math_ref1(int op, double a, double p) {
  double r;
  if (isnan(a)) return a;
  r = rsimd_math1_f64(op, a, p);
  if (isnan(r)) math_ref_status = RSIMD_EW_NAN_PRODUCED;
  return r;
}
static double math_ref2(int op, double a, double b) {
  double r;
  switch (op) {
  case RSIMD_MATH_POW: return rsimd_pow_f64(a, b);
  case RSIMD_MATH_ATAN2: r = rsimd_math2_na_f64(atan2(a, b), a, b); break;
  case RSIMD_MATH_NEXTAFTER: r = rsimd_math2_na_f64(nextafter(a, b), a, b); break;
  case RSIMD_MATH_REMAINDER: r = rsimd_math2_na_f64(remainder(a, b), a, b); break;
  case RSIMD_MATH_SCALEB:
    r = rsimd_math2_na_f64(isnan(b) ? b : rsimd_scaleb_f64(a, b), a, b);
    break;
  case RSIMD_MATH_ROOTN:
    r = rsimd_math2_na_f64(isnan(b) ? b : rsimd_rootn_f64(a, b), a, b);
    break;
  default: r = rsimd_math2_na_f64(hypot(a, b), a, b); break;
  }
  if (isnan(r) && !isnan(a) && !isnan(b)) math_ref_status = RSIMD_EW_NAN_PRODUCED;
  return r;
}
/* Base R's exact cases of x ^ y (rsimd_rpow_f64): y = 3 and 4 are
   products only for |x| <= 11. */
static int math_exact_pow(double a, double b) {
  return b == 0 || b == 2 || ((b == 3 || b == 4) && fabs(a) <= 11) || a == 0 || a == 1 ||
         !isfinite(a) || !isfinite(b);
}
static void check_math(const char *what, ptrdiff_t n, const double *got, const double *want,
                       const double *a, const double *b, int exact_pow, double bound) {
  ptrdiff_t i;
  char buf[200];
  for (i = 0; i < n; i++) {
    double g = got[i], w = want[i];
    int ok;
    n_checks++;
    if (isnan(w) || isnan(g)) {
      ok = isnan(g) && isnan(w) && is_na_ref(g) == is_na_ref(w);
    } else if (isinf(w) || w == 0 || (exact_pow && math_exact_pow(a[i], b[i]))) {
      ok = bits(g) == bits(w);
    } else {
      ok = fabs(g - w) <= bound * ulp_of(w) || (fabs(w) < DBL_MIN && fabs(g - w) <= DBL_MIN);
    }
    if (!ok) {
      snprintf(buf, sizeof buf, "input %a %a: got %a want %a", a[i], b ? b[i] : 0.0, g, w);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (bits(got[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}

static double ma[N + 1], mb[N + 1], mref[N + 1], mref2[N + 1], mout2[N + 1];
static int32_t mia[N + 1], mib[N + 1];

/* Inputs of op: specials and NA first, then values over its domain. */
static void math_inputs(int op, double *v, int32_t *iv, ptrdiff_t n) {
  double lo = -20, hi = 20;
  int log_scale = 0, j;
  switch (op) {
  case RSIMD_MATH_EXP: case RSIMD_MATH_EXPM1: lo = -746; hi = 710; break;
  case RSIMD_MATH_EXP2: lo = -1076; hi = 1025; break;
  case RSIMD_MATH_EXP10: lo = -324; hi = 309; break;
  case RSIMD_MATH_LOG: case RSIMD_MATH_LOG2: case RSIMD_MATH_LOG10: case RSIMD_MATH_LOGB:
  case RSIMD_MATH_CBRT: case RSIMD_MATH_ASINH: case RSIMD_MATH_ACOSH:
    lo = -320; hi = 308; log_scale = 1; break;
  case RSIMD_MATH_LOG1P: lo = -1.5; hi = 1e3; break;
  case RSIMD_MATH_ASIN: case RSIMD_MATH_ACOS: case RSIMD_MATH_ATANH: lo = -1.1; hi = 1.1; break;
  case RSIMD_MATH_SINPI: case RSIMD_MATH_COSPI: case RSIMD_MATH_TANPI:
    lo = -2; hi = 2; break;
  case RSIMD_MATH_SINH: case RSIMD_MATH_COSH: lo = -712; hi = 712; break;
  case RSIMD_MATH_SIGMOID: lo = -800; hi = 800; break;
  default: break;
  }
  sleef_inputs(v, n, lo, hi, log_scale, 1);
  if (n > 9) v[9] = rsimd_na_real();
  /* Multiples of 1/4, huge and tiny values for the pi functions. */
  if (op >= RSIMD_MATH_SINPI && op <= RSIMD_MATH_TANPI) {
    static const double extra[] = {0.25, -0.25, 0.75, 1.5, -1.5, 2, -2, 3, 1e300, -1e17,
                                   9007199254740993.0, 2.5e8 + 0.5, 1e-310, -5e-324};
    for (j = 0; j < (int) (sizeof extra / sizeof extra[0]) && 10 + j < n; j++) v[10 + j] = extra[j];
    for (j = 24; j < n && j < 120; j++) v[j] = (j - 72) * 0.25;
  }
  if (op == RSIMD_MATH_SINH || op == RSIMD_MATH_COSH) {
    static const double extra[] = {709.5, -709.9, 710.4, -710.47, 710.5, 750};
    for (j = 0; j < 6 && 10 + j < n; j++) v[10 + j] = extra[j];
  }
  for (j = 0; j < n; j++) iv[j] = (int32_t) (next_rand() % 41) - 20;
  if (n > 3) iv[3] = RSIMD_NA_I32;
}

/* The softmax passes (src/kernels/ml.inc.c) against libm and the scalar
   folds of na.h: exp(x - m) within 2 ULP and its sum within 1e-14
   (relative; the sum is at least 1) in every precision mode, SUM giving
   exactly the sum of EXP_SUM, DIV and LOG exactly (one or two IEEE
   operations), and every op the same in place as into another buffer. */
static double sx[N + 1], sref[N + 1], sout[N + 1];
static void test_softmax(void) {
  static const ptrdiff_t lens[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 65, 127, 128, 129, N};
  size_t li;
  ptrdiff_t j, n;
  int prec;
  char what[96], buf[160];
  for (li = 0; li < sizeof lens / sizeof lens[0]; li++) {
    double m = -INFINITY;
    n = lens[li];
    for (j = 0; j < n; j++) sx[j] = (double) (next_rand() % 1000000) * -6e-5 + 3;
    if (n > 3) sx[3] = -INFINITY;
    if (n > 5) sx[5] = -800; /* exp(x - m) subnormal or 0 */
    if (n > 6) sx[6] = -743;
    for (j = 0; j < n; j++) m = sx[j] > m ? sx[j] : m;
    for (j = 0; j < n; j++) sref[j] = exp(sx[j] - m);
    for (prec = RSIMD_PREC_FAST; prec <= RSIMD_PREC_COMPENSATED; prec++) {
      const rsimd_opts o = {0, 0, prec};
      rsimd_reduce_result r, r2, rr;
      double got, want;
      memset(&r, 0, sizeof r);
      memset(&r2, 0, sizeof r2);
      memset(&rr, 0, sizeof rr);
      snprintf(what, sizeof what, "softmax exp_sum precision %d", prec);
      reset_out();
      RSIMD_KERNEL(softmax_f64)(RSIMD_SOFTMAX_EXP_SUM, sx, n, m, 0.0, fout, &r, &o);
      check_math(what, n, fout, sref, sx, NULL, 0, 2);
      rsimd_fold_f64(sref, NULL, n, RSIMD_TERM_X, &rr, &o);
      got = rsimd_reduce_value(&r, prec);
      want = rsimd_reduce_value(&rr, prec);
      n_checks++;
      if (!(fabs(got - want) <= 1e-14 * want)) {
        snprintf(buf, sizeof buf, "sum got %a want %a", got, want);
        fail(what, n, 0, buf);
      }
      /* SUM writes nothing and gives the same sum. */
      snprintf(what, sizeof what, "softmax sum precision %d", prec);
      memcpy(sout, fout, sizeof sout);
      RSIMD_KERNEL(softmax_f64)(RSIMD_SOFTMAX_SUM, sx, n, m, 0.0, NULL, &r2, &o);
      check_int(what, n, 0, bits(rsimd_reduce_value(&r2, prec)) == bits(got), 1);
      /* In place. */
      snprintf(what, sizeof what, "softmax exp_sum in place precision %d", prec);
      memcpy(fout, sx, sizeof(double) * (size_t) n);
      memset(&r2, 0, sizeof r2);
      RSIMD_KERNEL(softmax_f64)(RSIMD_SOFTMAX_EXP_SUM, fout, n, m, 0.0, fout, &r2, &o);
      for (j = 0; j < n; j++) check_int(what, n, j, bits(fout[j]) == bits(sout[j]), 1);
      check_int(what, n, 0, bits(rsimd_reduce_value(&r2, prec)) == bits(got), 1);
    }
    {
      const rsimd_opts o = {0, 0, RSIMD_PREC_FAST};
      const double c = 1.2345678901234567;
      int op;
      for (op = RSIMD_SOFTMAX_DIV; op <= RSIMD_SOFTMAX_LOG; op++) {
        snprintf(what, sizeof what, "softmax op %d", op);
        reset_out();
        RSIMD_KERNEL(softmax_f64)(op, sx, n, m, c, fout, NULL, &o);
        for (j = 0; j < n; j++) {
          double w = op == RSIMD_SOFTMAX_DIV ? sx[j] / c : (sx[j] - m) - c;
          check_int(what, n, j, bits(fout[j]) == bits(w), 1);
        }
        n_checks++;
        if (bits(fout[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
        memcpy(sout, sx, sizeof(double) * (size_t) n);
        RSIMD_KERNEL(softmax_f64)(op, sout, n, m, c, sout, NULL, &o);
        for (j = 0; j < n; j++) check_int(what, n, j, bits(sout[j]) == bits(fout[j]), 1);
      }
    }
  }
}

/* log2 and exp2 of the powers of 2 and integers, log10 and exp10 of
   10^(0:22) and 0:22: exact, as in base R, in both accuracy modes. */
static void test_math_exact(void) {
  static double xs[2100], got[2100], want[2100];
  static const struct { int op, lo, hi, pow; } cases[] = {
    {RSIMD_MATH_LOG2, -1074, 1023, 2}, {RSIMD_MATH_EXP2, -1074, 1023, 0},
    {RSIMD_MATH_LOG10, 0, 22, 10}, {RSIMD_MATH_EXP10, 0, 22, 0}};
  size_t c;
  int k, m, fast;
  char buf[96];
  for (fast = 0; fast < 2; fast++) {
    for (c = 0; c < sizeof cases / sizeof cases[0]; c++) {
      for (m = 0, k = cases[c].lo; k <= cases[c].hi; k++, m++) {
        double p = cases[c].op == RSIMD_MATH_LOG10 || cases[c].op == RSIMD_MATH_EXP10
                     ? pow(10.0, k) : ldexp(1.0, k);
        xs[m] = cases[c].pow ? p : (double) k;
        want[m] = cases[c].pow ? (double) k : p;
      }
      RSIMD_KERNEL(math1_f64)(cases[c].op | (fast ? RSIMD_MATH_FAST : 0), xs, m, 0, 1.0, got);
      for (k = 0; k < m; k++) {
        n_checks++;
        if (bits(got[k]) != bits(want[k])) {
          snprintf(buf, sizeof buf, "input %a: got %a want %a", xs[k], got[k], want[k]);
          fail(fast ? "math exact powers fast" : "math exact powers", m, k, buf);
        }
      }
    }
  }
}

static void test_math_mode(int fast) {
  static const ptrdiff_t lens[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 65, N};
  const int bit = fast ? RSIMD_MATH_FAST : 0;
  const double extra = fast ? 2.5 : 0;
  int op, f, st;
  size_t li;
  ptrdiff_t j, n;
  char what[96];
  for (li = 0; li < sizeof lens / sizeof lens[0]; li++) {
    n = lens[li];
    for (op = RSIMD_MATH_EXP; op <= RSIMD_MATH_SIGMOID; op++) {
      const double p = op == RSIMD_MATH_LOGB ? log(3.0) : 1.0;
      math_inputs(op, ma, mia, n);
      for (f = 0; f < 2; f++) {
        reset_out();
        st = RSIMD_KERNEL(math1_f64)(op | bit, f ? (const void *) mia : ma, n,
                                     f ? RSIMD_EW_I32(0) : 0, p, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          mref[j] = math_ref1(op, f ? (mia[j] == RSIMD_NA_I32 ? rsimd_na_real() : mia[j]) : ma[j], p);
        }
        snprintf(what, sizeof what, "math1 op %d int32 %d fast %d", op, f, fast);
        /* glibc's cbrt is up to 3 ULP from the exact value on aarch64. */
        check_math(what, n, fout, mref, ma, NULL, 0, (op == RSIMD_MATH_CBRT ? 4 : 2) + extra);
        check_int(what, n, 0, st, math_ref_status);
      }
    }
    for (op = RSIMD_MATH_POW; op <= RSIMD_MATH_HYPOT; op++) {
      static const double pa[] = {-8, -2, -0.0, 0, 1, 2, 10, 11, -11, 12, INFINITY, -INFINITY, NAN};
      static const double pb[] = {2, 3, 4, 0.5, 1.0 / 3, 0, -1, -2, INFINITY, -INFINITY, NAN};
      sleef_inputs(ma, n, -300, 300, 1, 1);
      sleef_inputs(mb, n, -300, 300, 1, 1);
      if (op == RSIMD_MATH_POW) {
        for (j = 0; j < n; j++) {
          if (j % 3 == 0) ma[j] = pa[(j / 3) % 13];
          if (j % 2 == 0) mb[j] = pb[(j / 2) % 11];
          if (j % 5 == 4) ma[j] = ma[j] < 0 ? -fmod(-ma[j], 30) : fmod(ma[j], 30);
          if (j % 7 == 6) mb[j] = fmod(mb[j], 8);
        }
      }
      if (n > 9) ma[9] = rsimd_na_real();
      if (n > 12) mb[12] = rsimd_na_real();
      for (j = 0; j < n; j++) {
        mia[j] = (int32_t) (next_rand() % 21) - 10;
        mib[j] = (int32_t) (next_rand() % 9) - 4;
      }
      if (n > 5) mia[5] = RSIMD_NA_I32;
      /* f: bits 0-1 broadcast x, y; bits 2-3 int32 x, y. */
      for (f = 0; f < 16; f++) {
        const void *x = (f & 4) ? (const void *) mia : ma, *y = (f & 8) ? (const void *) mib : mb;
        int flags = ((f & 1) ? RSIMD_EW_SCALAR(0) : 0) | ((f & 2) ? RSIMD_EW_SCALAR(1) : 0) |
                    ((f & 4) ? RSIMD_EW_I32(0) : 0) | ((f & 8) ? RSIMD_EW_I32(1) : 0);
        double *xa = mref2, *yb = mout2;
        reset_out();
        st = RSIMD_KERNEL(math2_f64)(op | bit, x, y, n, flags, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          ptrdiff_t jx = (f & 1) ? 0 : j, jy = (f & 2) ? 0 : j;
          xa[j] = (f & 4) ? (mia[jx] == RSIMD_NA_I32 ? rsimd_na_real() : mia[jx]) : ma[jx];
          yb[j] = (f & 8) ? (mib[jy] == RSIMD_NA_I32 ? rsimd_na_real() : mib[jy]) : mb[jy];
          mref[j] = math_ref2(op, xa[j], yb[j]);
        }
        snprintf(what, sizeof what, "math2 op %d flags %d fast %d", op, f, fast);
        check_math(what, n, fout, mref, xa, yb, op == RSIMD_MATH_POW, 2 + extra);
        check_int(what, n, 0, st, math_ref_status);
      }
    }
    /* sincos: both outputs. */
    sleef_inputs(ma, n, -1e3, 1e3, 0, 1);
    if (n > 9) ma[9] = rsimd_na_real();
    reset_out();
    for (j = 0; j <= n; j++) mout2[j] = SENTINEL_F64;
    st = RSIMD_KERNEL(sincos_f64)(RSIMD_MATH_SIN | bit, ma, n, 0, fout, mout2);
    math_ref_status = 0;
    for (j = 0; j < n; j++) {
      mref[j] = math_ref1(RSIMD_MATH_SIN, ma[j], 1);
      mref2[j] = math_ref1(RSIMD_MATH_COS, ma[j], 1);
    }
    check_math("sincos sin", n, fout, mref, ma, NULL, 0, 2 + extra);
    check_math("sincos cos", n, mout2, mref2, ma, NULL, 0, 2 + extra);
    check_int("sincos status", n, 0, st, math_ref_status);
  }
}

static void test_math(void) {
  test_math_mode(0);
  test_math_mode(1);
}
#endif /* RSIMD_HAVE_SLEEF */

/* ---- Math extras ------------------------------------------------------------
   ilogb_f64 against C's ilogb, with every broadcast and int32 operand; the
   layer's recip_approx and rsqrt_approx within 2^-22 over a sweep of
   their range; on the tiers with SLEEF the math1 ops NEXT_UP ..
   RSQRT_APPROX and the math2 ops NEXTAFTER .. ROOTN against the scalar
   forms of math.inc.c (exactly, except the approximations, within 2^-22
   of the exact value, and rootn, within 1 ULP of a long double
   reference); on the none tier the scalar rootn against long double; and
   on every tier exact powers b^n, whose root must be exactly b. */
#ifndef RSIMD_NO_F64_SIMD

static double approx_worst[2] = {0, 0}; /* relative error: recip, rsqrt */
static double rootn_worst = 0;          /* ULP */
/* rootn's bound against the long double reference, which is exact enough
   (to about 0.15 ULP for 64-bit long double) except where long double is
   double. */
#define ROOTN_BOUND (LDBL_MANT_DIG >= 113 ? 1.0 : 1.15)

/* A finite double from random bits: every binade (and the subnormals)
   about equally likely, both signs. */
static double random_double(void) {
  double v;
  do {
    v = from_bits(next_rand());
  } while (!isfinite(v));
  return v;
}

/* The layer ops over normal numbers 2^-1022 <= |x| < 2^1022 (recip with
   both signs), through the predicated loop. */
static void test_approx_layer(void) {
  static double v[N + 1], r1[N + 1], r2[N + 1];
  static const double edge[] = {0x1p-1022,  0x1.fffffffffffffp1021, 1.0,  0x1.fffffffffffffp0,
                                2.0,        0x1.0000000000001p0,    3.0,  0x1.8p-1000,
                                0x1p1000,   1e-300,                 1e300, 0.1,
                                7.0,        0x1.fffffffffffffp-1,   4.0,  0x1.fffffffffffffp1,
                                0x1.6a09e667f3bcdp0, 0x1.6a09e667f3bcdp1};
  int rep;
  ptrdiff_t i, j;
  char buf[160];
  for (rep = 0; rep < 100; rep++) {
    const double sign = rep % 2 ? -1.0 : 1.0;
    for (j = 0; j < N; j++) {
      if (rep < 2 && j < (ptrdiff_t) (sizeof edge / sizeof edge[0])) {
        v[j] = edge[j];
      } else {
        v[j] = ldexp(1 + (double) (next_rand() >> 11) * 0x1p-53, (int) (next_rand() % 2044) - 1022);
      }
    }
    for (i = 0; i + RSIMD_LANES_64 <= N; i += RSIMD_LANES_64) {
      rsimd_vf64 a = rsimd_vf64_loadu(v + i);
      rsimd_vf64_storeu(r1 + i, rsimd_vf64_recip_approx(rsimd_vf64_mul(rsimd_vf64_set1(sign), a)));
      rsimd_vf64_storeu(r2 + i, rsimd_vf64_rsqrt_approx(a));
    }
    if (i < N) {
      rsimd_p64 pg = rsimd_p64_while(i, N);
      rsimd_vf64 a = rsimd_vf64_loadu_p(pg, v + i, v[i]);
      rsimd_vf64_storeu_p(pg, r1 + i,
                          rsimd_vf64_recip_approx(rsimd_vf64_mul(rsimd_vf64_set1(sign), a)));
      rsimd_vf64_storeu_p(pg, r2 + i, rsimd_vf64_rsqrt_approx(a));
    }
    for (j = 0; j < N; j++) {
      long double s = sqrtl((long double) v[j]);
      double e1 = (double) fabsl(((long double) r1[j] * sign) * v[j] - 1.0L);
      double e2 = (double) fabsl((long double) r2[j] * s - 1.0L);
      n_checks += 2;
      if (e1 > approx_worst[0]) approx_worst[0] = e1;
      if (e2 > approx_worst[1]) approx_worst[1] = e2;
      if (!(e1 <= 0x1p-22)) {
        snprintf(buf, sizeof buf, "input %a: got %a (error %g)", sign * v[j], r1[j], e1);
        fail("layer recip_approx", N, j, buf);
      }
      if (!(e2 <= 0x1p-22)) {
        snprintf(buf, sizeof buf, "input %a: got %a (error %g)", v[j], r2[j], e2);
        fail("layer rsqrt_approx", N, j, buf);
      }
    }
  }
}

/* ilogb_f64 with double and int32 input, broadcast or not. */
static double ilx[N + 1];
static int32_t ilxi[N + 1], ilout[N + 1];
static void test_ilogb(ptrdiff_t n) {
  static const double specials[] = {0.0,       -0.0, INFINITY, -INFINITY, NAN, 5e-324, -5e-324,
                                    0x1p-1022, DBL_MAX, 1.0, -1.0, 0x1.fffffffffffffp-1023, 0.75};
  const int ns = (int) (sizeof specials / sizeof specials[0]);
  ptrdiff_t j;
  int f;
  char what[64];
  for (j = 0; j < n; j++) {
    ilx[j] = j < ns ? specials[j] : random_double();
    ilxi[j] = j % 5 == 1 ? 0 : (int32_t) next_rand();
  }
  if (n > ns) ilx[ns] = rsimd_na_real();
  if (n > 3) ilxi[3] = RSIMD_NA_I32;
  for (f = 0; f < 4; f++) {
    const int flags = ((f & 1) ? RSIMD_EW_SCALAR(0) : 0) | ((f & 2) ? RSIMD_EW_I32(0) : 0);
    for (j = 0; j <= n; j++) ilout[j] = SENTINEL_I32;
    RSIMD_KERNEL(ilogb_f64)((f & 2) ? (const void *) ilxi : ilx, n, flags, ilout);
    snprintf(what, sizeof what, "ilogb_f64 flags %d", flags);
    for (j = 0; j < n; j++) {
      ptrdiff_t k = (f & 1) ? 0 : j;
      double a = (f & 2) ? (ilxi[k] == RSIMD_NA_I32 ? rsimd_na_real() : ilxi[k]) : ilx[k];
      int32_t want = (isnan(a) || isinf(a) || a == 0) ? RSIMD_NA_I32 : ilogb(a);
      check_int(what, n, j, ilout[j], want);
    }
    check_int(what, n, n, ilout[n], SENTINEL_I32);
  }
}

static long double ref_rootn(double x, double n) {
  long double r = powl(fabsl((long double) x), 1.0L / (long double) n);
  return x < 0 ? -r : r;
}

/* The ULP error of got (for x, n) against long double, checked against
   ROOTN_BOUND where the result is finite and nonzero. */
static void check_rootn_ulp(const char *what, ptrdiff_t n, ptrdiff_t i, double x, double nn,
                            double got) {
  long double ref;
  double w, err;
  char buf[200];
  if (LDBL_MANT_DIG < 64 || isnan(got) || isinf(got) || got == 0) return;
  ref = ref_rootn(x, nn);
  w = (double) ref;
  /* The ULP of a subnormal is the smallest subnormal. */
  err = (double) (fabsl((long double) got - ref) / fmax(ldexp(1.0, ilogb(w) - 52), DBL_TRUE_MIN));
  n_checks++;
  if (err > rootn_worst) rootn_worst = err;
  if (!(err <= ROOTN_BOUND)) {
    snprintf(buf, sizeof buf, "rootn(%a, %.0f): got %a, %.3f ULP from %.20Lg", x, nn, got, err, ref);
    fail(what, n, i, buf);
  }
}

/* The roots of exact powers: rootn(b^n, n) is b and rootn(-b^n, n) is -b
   for odd n, for n in 3..60 and b up to 3000 (and b^n exactly
   representable, below 2^1000), with the scalar form, and on the tiers with
   SLEEF with the math2 kernel too. */
static double ppx[N + 1], ppn[N + 1], ppb[N + 1];
#ifdef RSIMD_HAVE_SLEEF
static double ppout[N + 1];
#endif
static void rootn_exact_batch(ptrdiff_t m) {
  ptrdiff_t j;
  char buf[160];
  for (j = 0; j < m; j++) {
    double r = rsimd_rootn_f64(ppx[j], ppn[j]);
    n_checks++;
    if (bits(r) != bits(ppb[j])) {
      snprintf(buf, sizeof buf, "rootn(%a, %.0f): got %a want %a", ppx[j], ppn[j], r, ppb[j]);
      fail("rootn exact powers scalar", m, j, buf);
    }
  }
#ifdef RSIMD_HAVE_SLEEF
  ppout[m] = SENTINEL_F64;
  RSIMD_KERNEL(math2_f64)(RSIMD_MATH_ROOTN, ppx, ppn, m, 0, ppout);
  for (j = 0; j < m; j++) {
    n_checks++;
    if (bits(ppout[j]) != bits(ppb[j])) {
      snprintf(buf, sizeof buf, "rootn(%a, %.0f): got %a want %a", ppx[j], ppn[j], ppout[j],
               ppb[j]);
      fail("rootn exact powers kernel", m, j, buf);
    }
  }
#endif
}
static void test_rootn_exact(void) {
  ptrdiff_t m = 0;
  int n, b, k;
  for (n = 3; n <= 60; n++) {
    for (b = 2; b <= 3000; b++) {
      double p = b;
      for (k = 1; k < n; k++) {
        double q = p * b;
        if (fma(p, (double) b, -q) != 0 || q > 0x1p1000) break;
        p = q;
      }
      if (k < n) break;
      ppx[m] = p;
      ppn[m] = n;
      ppb[m] = b;
      if (++m == N) {
        rootn_exact_batch(m);
        m = 0;
      }
      if (n % 2) {
        ppx[m] = -p;
        ppn[m] = n;
        ppb[m] = -b;
        if (++m == N) {
          rootn_exact_batch(m);
          m = 0;
        }
      }
    }
  }
  if (m > 0) rootn_exact_batch(m);
}

#if RSIMD_TIER_IS(none)
/* The scalar rootn (the none tier's) against long double, over random
   finite x and a set of n. */
static void test_rootn_scalar(void) {
  static const double ns[] = {2,    3,    4,    5,    6,    7,     8,           9,
                              10,   11,   17,   64,   100,  511,   512,         513,
                              1000, 123457, 2147483647.0, -2, -3,  -4,          -5,
                              -7,   -10,  -100, -512, -513, -1000, -2147483647.0};
  size_t k;
  int j;
  for (k = 0; k < sizeof ns / sizeof ns[0]; k++) {
    const double nn = ns[k];
    const int odd = fmod(nn, 2.0) != 0;
    for (j = 0; j < 4000; j++) {
      double x = fabs(random_double());
      if (odd && (next_rand() & 1)) x = -x;
      check_rootn_ulp("scalar rootn", 4000, j, x, nn, rsimd_rootn_f64(x, nn));
    }
  }
}
#endif

#ifdef RSIMD_HAVE_SLEEF
/* got within relative error `bound` of want where want is a nonzero
   finite number, otherwise equal (NaN kinds compared); the worst error
   goes to *worst. */
static void check_rel(const char *what, ptrdiff_t n, const double *got, const double *want,
                      const double *a, double bound, double *worst) {
  ptrdiff_t i;
  char buf[200];
  for (i = 0; i < n; i++) {
    double g = got[i], w = want[i];
    int ok;
    n_checks++;
    if (isnan(w) || isnan(g) || isinf(w) || w == 0) {
      ok = bits(g) == bits(w) || (isnan(g) && isnan(w) && is_na_ref(g) == is_na_ref(w));
    } else {
      double e = (double) fabsl(((long double) g - w) / w);
      ok = e <= bound;
      if (e > *worst) *worst = e;
    }
    if (!ok) {
      snprintf(buf, sizeof buf, "input %a: got %a want %a", a[i], g, w);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (bits(got[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}

/* math1 extras over every length, double and int32 input, broadcast or
   not, with and without RSIMD_MATH_FAST (which they ignore). */
static void test_math1_extras(void) {
  static const ptrdiff_t lens[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 65, N};
  static const double specials[] = {0.0,    -0.0,    INFINITY,  -INFINITY, NAN,
                                    5e-324, -5e-324, 0x1p-1022, DBL_MAX,   -DBL_MAX,
                                    0x1p1022, -0x1p1022, 0x1.fffffffffffffp1021, 1.0, -1.0,
                                    0x1.fffffffffffffp-1023, 0x1p-1023};
  const int nsp = (int) (sizeof specials / sizeof specials[0]);
  size_t li;
  int op, f, st;
  ptrdiff_t j, n;
  char what[96];
  for (li = 0; li < sizeof lens / sizeof lens[0]; li++) {
    n = lens[li];
    for (j = 0; j < n; j++) {
      ma[j] = j < nsp ? specials[j] : random_double();
      mia[j] = (int32_t) (next_rand() % 2001) - 1000;
    }
    if (n > nsp) ma[nsp] = rsimd_na_real();
    if (n > 3) mia[3] = RSIMD_NA_I32;
    for (op = RSIMD_MATH_NEXT_UP; op <= RSIMD_MATH_RSQRT_APPROX; op++) {
      const int approx = op == RSIMD_MATH_RECIP_APPROX || op == RSIMD_MATH_RSQRT_APPROX;
      for (f = 0; f < 8; f++) {
        const int flags = ((f & 1) ? RSIMD_EW_SCALAR(0) : 0) | ((f & 2) ? RSIMD_EW_I32(0) : 0);
        reset_out();
        st = RSIMD_KERNEL(math1_f64)(op | ((f & 4) ? RSIMD_MATH_FAST : 0),
                                     (f & 2) ? (const void *) mia : ma, n, flags, 1.0, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          ptrdiff_t k = (f & 1) ? 0 : j;
          mref2[j] = (f & 2) ? (mia[k] == RSIMD_NA_I32 ? rsimd_na_real() : mia[k]) : ma[k];
          mref[j] = math_ref1(op, mref2[j], 1.0);
        }
        snprintf(what, sizeof what, "math1 op %d flags %d", op, f);
        if (approx) {
          check_rel(what, n, fout, mref, mref2, 0x1p-22,
                    &approx_worst[op == RSIMD_MATH_RSQRT_APPROX]);
        } else {
          check_math(what, n, fout, mref, mref2, NULL, 0, 0);
        }
        check_int(what, n, 0, st, math_ref_status);
      }
    }
  }
}

/* Inputs of the math2 extras: x in ma (and int32 in mia), y or n in mb
   (and mib). */
static void math2_extra_inputs(int op, ptrdiff_t n) {
  static const double sx[] = {0.0,  -0.0, INFINITY, -INFINITY, NAN,    5e-324,
                              -5e-324, 1.0, -1.0,  DBL_MAX,   -DBL_MAX, 0x1p-1022};
  static const double sy[] = {0.0, -0.0, INFINITY, -INFINITY, NAN, 1.0, 3.0, 5e-324, -0x1p-1022};
  static const double sn[] = {0, 1, -1, 2, -2, 3, -3, 4, 5, 1023, -1022, 1074, -1074, -1075,
                              2098, -2098, 2147483647.0, -2147483647.0, 512, 513, -513, 1000};
  const int nsx = (int) (sizeof sx / sizeof sx[0]), nsy = (int) (sizeof sy / sizeof sy[0]),
            nsn = (int) (sizeof sn / sizeof sn[0]);
  ptrdiff_t j;
  for (j = 0; j < n; j++) {
    double x = j < nsx ? sx[j] : random_double(), y;
    switch (op) {
    case RSIMD_MATH_NEXTAFTER:
      y = j % 4 == 0 ? x : j % 4 == 1 ? sy[(j / 4) % nsy] : random_double();
      break;
    case RSIMD_MATH_REMAINDER:
      y = j % 5 == 0 ? sy[(j / 5) % nsy]
          : j % 5 == 1 ? x / (double) (1 + next_rand() % 1000)
          : j % 5 == 2 ? ldexp(x, -(int) (next_rand() % 60))
                       : random_double();
      break;
    case RSIMD_MATH_SCALEB:
      y = j % 3 == 0 ? sn[(j / 3) % nsn] : (double) ((int) (next_rand() % 4401) - 2200);
      break;
    default: /* ROOTN: x log-uniform, some exact powers of the small n */
      if (j >= nsx) {
        x = pow(10.0, -320 + (double) (next_rand() >> 11) * 0x1p-53 * 628);
        if (next_rand() & 1) x = -x;
      }
      y = j % 3 == 0 ? sn[(j / 3) % nsn] : (double) ((int) (next_rand() % 81) - 40);
      if (j % 7 == 5) x = pow(3.0, fabs(y) > 30 ? 3 : fabs(y)) * (x < 0 ? -1 : 1);
      break;
    }
    ma[j] = x;
    mb[j] = y;
    mia[j] = (int32_t) (next_rand() % 2001) - 1000;
    mib[j] = (op == RSIMD_MATH_SCALEB || op == RSIMD_MATH_ROOTN)
               ? (int32_t) y
               : (int32_t) (next_rand() % 21) - 10;
  }
  if (n > 9) ma[9] = rsimd_na_real();
  if (n > 12) mb[12] = rsimd_na_real();
  if (n > 5) mia[5] = RSIMD_NA_I32;
  if (n > 7) mib[7] = RSIMD_NA_I32;
}

static void test_math2_extras(void) {
  static const ptrdiff_t lens[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 65, N};
  size_t li;
  int op, f, st;
  ptrdiff_t j, n;
  char what[96];
  for (li = 0; li < sizeof lens / sizeof lens[0]; li++) {
    n = lens[li];
    for (op = RSIMD_MATH_NEXTAFTER; op <= RSIMD_MATH_ROOTN; op++) {
      math2_extra_inputs(op, n);
      /* f: bits 0-1 broadcast x, y; bits 2-3 int32 x, y; bit 4 RSIMD_MATH_FAST. */
      for (f = 0; f < 32; f++) {
        const void *x = (f & 4) ? (const void *) mia : ma, *y = (f & 8) ? (const void *) mib : mb;
        int flags = ((f & 1) ? RSIMD_EW_SCALAR(0) : 0) | ((f & 2) ? RSIMD_EW_SCALAR(1) : 0) |
                    ((f & 4) ? RSIMD_EW_I32(0) : 0) | ((f & 8) ? RSIMD_EW_I32(1) : 0);
        double *xa = mref2, *yb = mout2;
        reset_out();
        st = RSIMD_KERNEL(math2_f64)(op | ((f & 16) ? RSIMD_MATH_FAST : 0), x, y, n, flags, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          ptrdiff_t jx = (f & 1) ? 0 : j, jy = (f & 2) ? 0 : j;
          xa[j] = (f & 4) ? (mia[jx] == RSIMD_NA_I32 ? rsimd_na_real() : mia[jx]) : ma[jx];
          yb[j] = (f & 8) ? (mib[jy] == RSIMD_NA_I32 ? rsimd_na_real() : mib[jy]) : mb[jy];
          mref[j] = math_ref2(op, xa[j], yb[j]);
        }
        snprintf(what, sizeof what, "math2 op %d flags %d", op, f);
        if (op == RSIMD_MATH_ROOTN) {
          /* Special values exactly as none, the rest within the bound. */
          for (j = 0; j < n; j++) {
            if (isnan(mref[j]) || isinf(mref[j]) || mref[j] == 0) {
              n_checks++;
              if (!(bits(fout[j]) == bits(mref[j]) ||
                    (isnan(fout[j]) && isnan(mref[j]) && is_na_ref(fout[j]) == is_na_ref(mref[j])))) {
                char buf[160];
                snprintf(buf, sizeof buf, "rootn(%a, %a): got %a want %a", xa[j], yb[j], fout[j],
                         mref[j]);
                fail(what, n, j, buf);
              }
            } else {
              check_rootn_ulp(what, n, j, xa[j], yb[j], fout[j]);
            }
          }
          n_checks++;
          if (bits(fout[n]) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
        } else {
          check_math(what, n, fout, mref, xa, yb, 0, 0);
        }
        check_int(what, n, 0, st, math_ref_status);
      }
    }
  }
}
#endif /* RSIMD_HAVE_SLEEF */

static void test_math_extras(void) {
  test_approx_layer();
  test_rootn_exact();
#if RSIMD_TIER_IS(none)
  test_rootn_scalar();
#endif
#ifdef RSIMD_HAVE_SLEEF
  test_math1_extras();
  test_math2_extras();
#endif
}
#endif /* RSIMD_NO_F64_SIMD */

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
    {
      int32_t got[2 * RSIMD_MAX_LANES_64];
      rsimd_vi32_storeu(got, rsimd_vi32_inc(x, m));
      for (j = 0; j < L; j++) {
        check_int("i32 inc", L, i + j, got[j],
                  (int32_t) ((uint32_t) ia[i + j] + (ia[i + j] > ib[i + j])));
      }
    }
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
    {
      int64_t got[RSIMD_MAX_LANES_64];
      rsimd_vi64_storeu(got, rsimd_vi64_inc(x, m));
      for (j = 0; j < L; j++) {
        n_checks++;
        if (got[j] != (int64_t) ((uint64_t) la[i + j] + (la[i + j] > lb[i + j]))) {
          fail("i64 inc", L, i + j, "mismatch");
        }
      }
    }
    check_int("mi64 any", L, i, rsimd_mi64_any(m), cnt > 0);
    check_int("mi64 all", L, i, rsimd_mi64_all(m), cnt == L);
    check_int("mi64 all (true)", L, i, rsimd_mi64_all(rsimd_vi64_cmp_eq(x, x)), 1);
    n_checks += 3;
    if (rsimd_vi64_reduce_add(x) != (int64_t) s) fail("i64 reduce_add", L, i, "mismatch");
    if (rsimd_vi64_reduce_min(x) != lo) fail("i64 reduce_min", L, i, "mismatch");
    if (rsimd_vi64_reduce_max(x) != hi) fail("i64 reduce_max", L, i, "mismatch");
  }
}

/* x %/% d and x %% d by multiplication (rsimd_ew_intdiv_const_i32) against
   the scalar rsimd_intdiv_i32(), with and without the NA check, for small
   divisors, powers of two and their neighbours, the int32 extremes and
   random ones, on dividends around the multiples of d, the extremes and
   random ones; lengths 1 .. 2 vectors + 1 cover the predicated tail. */
static void test_intdiv_const(void) {
  static int32_t xs[N + 1], got[N + 1], want[N + 1];
  int32_t ds[2600];
  int nd = 0, k, check, mod;
  ptrdiff_t j, n;
  for (k = -1100; k <= 1100; k++) {
    if (k != 0) ds[nd++] = k;
  }
  for (k = 11; k < 31; k++) {
    int32_t p = (int32_t) 1 << k;
    ds[nd++] = p;
    ds[nd++] = -p;
    ds[nd++] = p - 1;
    ds[nd++] = -(p - 1);
    ds[nd++] = p + 1;
    ds[nd++] = -(p + 1);
  }
  ds[nd++] = INT32_MAX;
  ds[nd++] = -INT32_MAX;
  ds[nd++] = INT32_MAX - 1;
  ds[nd++] = -INT32_MAX + 1;
  while (nd < 2600) {
    int32_t d = (int32_t) (uint32_t) next_rand();
    if (d != 0 && d != INT32_MIN) ds[nd++] = d;
  }
  for (k = 0; k < nd; k++) {
    const int32_t d = ds[k];
    for (j = 0; j < N; j++) {
      switch (j % 8) {
      case 0: xs[j] = (int32_t) (uint32_t) next_rand(); break;
      case 1: xs[j] = (int32_t) ((uint32_t) d * (uint32_t) (j - 100)); break;
      case 2: xs[j] = (int32_t) ((uint32_t) d * (uint32_t) (j - 100) + 1u); break;
      case 3: xs[j] = (int32_t) ((uint32_t) d * (uint32_t) (j - 100) - 1u); break;
      case 4: xs[j] = (int32_t) (next_rand() % 4001) - 2000; break;
      case 5: xs[j] = (j / 8) % 2 ? INT32_MAX - (int32_t) (j / 16) : INT32_MIN + (int32_t) (j / 16); break;
      case 6: xs[j] = j % 3 == 0 ? RSIMD_NA_I32 : (int32_t) (j / 8) - 12; break;
      default: xs[j] = (int32_t) (uint32_t) (next_rand() >> (next_rand() % 32)); break;
      }
    }
    for (check = 0; check < 2; check++) {
      for (mod = 0; mod < 2; mod++) {
        char what[64];
        snprintf(what, sizeof what, "intdiv_const d=%ld mod %d check %d", (long) d, mod, check);
        for (j = 0; j < N; j++) want[j] = rsimd_intdiv_i32(xs[j], d, mod, check);
        n = k % 8 == 0 ? N : 1 + (ptrdiff_t) (k % (2 * RSIMD_LANES_32 + 1));
        for (j = 0; j <= N; j++) got[j] = SENTINEL_I32;
#if RSIMD_TIER_IS(none)
        {
          const rsimd_divmagic_i32 g = rsimd_divmagic_i32_make(d);
          for (j = 0; j < n; j++) got[j] = rsimd_vi32_intdiv_const(xs[j], &g, mod, check);
        }
#else
        rsimd_ew_intdiv_const_i32(xs, d, n, mod, check, got);
#endif
        check_i32(what, n, got, want);
      }
    }
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
          rsimd_scan_state s = {inits[op], 0.0, 0, 0};
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
      rsimd_scan_state s = {0.0, 0.0, 0, 0};
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

/* ---- Elementwise kernels ------------------------------------------------ */

#ifndef RSIMD_NO_F64_SIMD
/* Bit-identical, except that two NaNs match when both or neither are NA. */
static void check_f64_kind(const char *what, ptrdiff_t n, const double *got, const double *want) {
  ptrdiff_t i;
  char buf[160];
  for (i = 0; i < n; i++) {
    int ok = bits(got[i]) == bits(want[i]) ||
             (isnan(got[i]) && isnan(want[i]) && is_na_ref(got[i]) == is_na_ref(want[i]));
    n_checks++;
    if (!ok) {
      snprintf(buf, sizeof buf, "got %a (%016llx) want %a (%016llx)", got[i],
               (unsigned long long) bits(got[i]), want[i], (unsigned long long) bits(want[i]));
      fail(what, n, i, buf);
    }
  }
}

/* Element i of an elementwise operand: double or int32 (NA becoming
   NA_real_), broadcast or not. */
static double ew_get(const void *p, int i32, int scalar, ptrdiff_t i) {
  ptrdiff_t j = scalar ? 0 : i;
  if (i32) {
    int32_t v = ((const int32_t *) p)[j];
    return v == RSIMD_NA_I32 ? rsimd_na_real() : (double) v;
  }
  return ((const double *) p)[j];
}

/* Doubles for %% with huge and overflowing quotients and exact ties. */
static double fmodx[N + 1], fmody[N + 1];

static void init_mod_inputs(void) {
  const double xs[] = {1e20, -1e20, 1e308, -1e308, 3 + 0x1p-51, -(3 + 0x1p-51), 7.5, -0.0,
                       5e-324, 1e300, HUGE_VAL, 0x1p60};
  const double ys[] = {3, 0.1, 5e-324, -5e-324, 1 + 0x1p-52, 1 + 0x1p-52, -2.5, 2, 1e-300,
                       -1e-10, 2, 7};
  int i, k = (int) (sizeof xs / sizeof xs[0]);
  for (i = 0; i < N; i++) {
    fmodx[i] = i < 2 * k ? xs[i % k] : fa[i] * 1e15;
    fmody[i] = i < 2 * k ? ys[(i + i / k) % k] : fb[i];
  }
}

static double ref_ew1(int op, double a) {
  switch (op) {
  case RSIMD_EW_NEG: return -a;
  case RSIMD_EW_ABS: return fabs(a);
  case RSIMD_EW_SIGN: return rsimd_sign_f64(a);
  case RSIMD_EW_RECIP: return 1.0 / a;
  case RSIMD_EW_SQRT: return sqrt(a);
  case RSIMD_EW_FLOOR: return floor(a);
  case RSIMD_EW_CEIL: return ceil(a);
  case RSIMD_EW_TRUNC: return trunc(a);
  default: return nearbyint(a);
  }
}

static double ref_ew2(int op, double a, double b) {
  double r;
  switch (op) {
  case RSIMD_EW_ADD: r = a + b; break;
  case RSIMD_EW_SUB: r = a - b; break;
  case RSIMD_EW_MUL: r = a * b; break;
  case RSIMD_EW_DIV: r = a / b; break;
  case RSIMD_EW_IDIV: r = rsimd_idiv_f64(a, b); break;
  case RSIMD_EW_MOD: r = rsimd_mod_f64(a, b); break;
  case RSIMD_EW_PMIN: return rsimd_pmin_f64(a, b);
  case RSIMD_EW_PMAX: return rsimd_pmax_f64(a, b);
  case RSIMD_EW_PMIN_NUM: return rsimd_pmin_num_f64(a, b);
  case RSIMD_EW_PMAX_NUM: return rsimd_pmax_num_f64(a, b);
  default: r = rsimd_copysign_f64(a, b); break;
  }
  return rsimd_na_merge_f64(r, a, b);
}

static double ref_ew3(int op, double a, double b, double c) {
  switch (op) {
  case RSIMD_EW_FMA: return rsimd_na_merge3_f64(fma(a, b, c), a, b, c);
  case RSIMD_EW_MUL_ADD: return rsimd_na_merge3_f64(a * b + c, a, b, c);
  case RSIMD_EW_ADD_MUL: return rsimd_na_merge3_f64((a + b) * c, a, b, c);
  case RSIMD_EW_LERP: return rsimd_na_merge3_f64(fma(c, b, (1.0 - c) * a), a, b, c);
  /* Fused where the tier has a native fused multiply-add (none has not). */
  case RSIMD_EW_MUL_ADD_APPROX:
    return rsimd_na_merge3_f64(RSIMD_NATIVE_FMA ? fma(a, b, c) : a * b + c, a, b, c);
  default: return rsimd_pmin_f64(rsimd_pmax_f64(a, b), c);
  }
}

#endif /* RSIMD_NO_F64_SIMD */

static int32_t ref_ew2_i32(int op, int32_t a, int32_t b, int check, int *ovf) {
  switch (op) {
  case RSIMD_EW_ADD: return rsimd_add_i32(a, b, check, ovf);
  case RSIMD_EW_SUB: return rsimd_sub_i32(a, b, check, ovf);
  case RSIMD_EW_MUL: return rsimd_mul_i32(a, b, check, ovf);
  case RSIMD_EW_ADD_WRAP: return rsimd_add_wrap_i32(a, b, check);
  case RSIMD_EW_SUB_WRAP: return rsimd_sub_wrap_i32(a, b, check);
  case RSIMD_EW_MUL_WRAP: return rsimd_mul_wrap_i32(a, b, check);
  case RSIMD_EW_IDIV: return rsimd_intdiv_i32(a, b, 0, check);
  case RSIMD_EW_MOD: return rsimd_intdiv_i32(a, b, 1, check);
  case RSIMD_EW_PMIN: return rsimd_pmin_i32(a, b);
  case RSIMD_EW_PMAX: return rsimd_pmax_i32(a, b);
  case RSIMD_EW_PMIN_NUM: return rsimd_pmin_num_i32(a, b);
  default: return rsimd_pmax_num_i32(a, b);
  }
}

static int32_t ref_ew3_i32(int op, int32_t a, int32_t b, int32_t c, int check, int *ovf,
                           int *lohi) {
  int32_t t;
  switch (op) {
  case RSIMD_EW_MUL_ADD:
    t = rsimd_mul_i32(a, b, check, ovf);
    return t == RSIMD_NA_I32 ? t : rsimd_add_i32(t, c, check, ovf);
  case RSIMD_EW_ADD_MUL:
    t = rsimd_add_i32(a, b, check, ovf);
    return t == RSIMD_NA_I32 ? t : rsimd_mul_i32(t, c, check, ovf);
  default:
    if (b != RSIMD_NA_I32 && c != RSIMD_NA_I32 && b > c) *lohi = 1;
    return rsimd_pmin_i32(rsimd_pmax_i32(a, b), c);
  }
}

static void test_arith(ptrdiff_t n) {
  const rsimd_opts o = {0, 1, RSIMD_PREC_FAST}, o_nocheck = {0, 0, RSIMD_PREC_FAST};
  int op, f, st, want_st, check;
  ptrdiff_t j;
  char what[96];
  if (n == 0) return; /* kernels are never called for no elements */
  /* int32 kernels: every op, each operand broadcast or not. */
  /* int32 kernels: every op, each operand broadcast or not, with and
     without NA checks. */
  for (check = 0; check < 2; check++) {
    const rsimd_opts *oc = check ? &o : &o_nocheck;
    for (op = RSIMD_EW_ADD; op <= RSIMD_EW_MUL_WRAP; op++) {
      if (op == RSIMD_EW_DIV || op == RSIMD_EW_COPYSIGN) continue;
      for (f = 0; f < 3; f++) {
        const int32_t *x = op == RSIMD_EW_MUL || op == RSIMD_EW_MUL_WRAP ? ismall : ia;
        const int32_t *y = op == RSIMD_EW_MUL || op == RSIMD_EW_MUL_WRAP ? ismall2 : ib;
        int ovf = 0;
        reset_out();
        st = RSIMD_KERNEL(ew2_i32)(op, x, y, n, f, iout, oc);
        for (j = 0; j < n; j++) {
          iref[j] = ref_ew2_i32(op, x[f & 1 ? 0 : j], y[f & 2 ? 0 : j], check, &ovf);
        }
        snprintf(what, sizeof what, "ew2_i32 op %d flags %d check %d", op, f, check);
        check_i32(what, n, iout, iref);
        /* Without checks, NA operands may or may not count as overflow. */
        if (check) check_int(what, n, 0, st, ovf ? RSIMD_EW_OVERFLOW : 0);
      }
    }
    for (op = RSIMD_EW_MUL_ADD; op <= RSIMD_EW_CLAMP; op++) {
      if (op == RSIMD_EW_LERP) continue;
      for (f = 0; f < 7; f++) {
        const int32_t *z = op == RSIMD_EW_CLAMP ? ib : ia;
        int ovf = 0, lohi = 0;
        reset_out();
        st = RSIMD_KERNEL(ew3_i32)(op, ismall, ismall2, z, n, f, iout, oc);
        for (j = 0; j < n; j++) {
          iref[j] = ref_ew3_i32(op, ismall[f & 1 ? 0 : j], ismall2[f & 2 ? 0 : j],
                                z[f & 4 ? 0 : j], check, &ovf, &lohi);
        }
        snprintf(what, sizeof what, "ew3_i32 op %d flags %d check %d", op, f, check);
        check_i32(what, n, iout, iref);
        if (check) {
          check_int(what, n, 0, st, (ovf ? RSIMD_EW_OVERFLOW : 0) | (lohi ? RSIMD_EW_LO_GT_HI : 0));
        }
      }
    }
  }
  /* The inactive lanes of the last vector must not raise a status bit:
     (x + 65536) * 32768 is 0 for x = -65536 but overflows for x = 0, and
     clamp's broadcast lo = 5 is above 0 but not above hi = 10. */
  {
    static int32_t xm[N + 1];
    const int32_t y1 = 65536, z1 = 32768;
    for (j = 0; j < n; j++) xm[j] = -65536;
    st = RSIMD_KERNEL(ew3_i32)(RSIMD_EW_ADD_MUL, xm, &y1, &z1, n,
                               RSIMD_EW_SCALAR(1) | RSIMD_EW_SCALAR(2), iout, &o);
    check_int("ew3_i32 tail status", n, 0, st, 0);
    st = RSIMD_KERNEL(ew2_i32)(RSIMD_EW_MUL, xm, &y1, n, RSIMD_EW_SCALAR(1), iout, &o);
    check_int("ew2_i32 tail status", n, 0, st, RSIMD_EW_OVERFLOW);
#ifndef RSIMD_NO_F64_SIMD
    {
      static double hi[N + 1];
      const double lo = 5.0;
      for (j = 0; j < n; j++) hi[j] = 10.0;
      st = RSIMD_KERNEL(ew3_f64)(RSIMD_EW_CLAMP, fa, &lo, hi, n, RSIMD_EW_SCALAR(1), fout, &o);
      check_int("ew3_f64 tail status", n, 0, st, 0);
    }
#endif
  }
  for (op = RSIMD_EW_NEG; op <= RSIMD_EW_ABS; op++) {
    reset_out();
    RSIMD_KERNEL(ew1_i32)(op, ia, n, iout, &o);
    for (j = 0; j < n; j++) iref[j] = op == RSIMD_EW_ABS ? rsimd_abs_i32(ia[j]) : rsimd_neg_i32(ia[j]);
    check_i32(op == RSIMD_EW_ABS ? "ew1_i32 abs" : "ew1_i32 neg", n, iout, iref);
  }
#ifndef RSIMD_NO_F64_SIMD
  /* double kernels: every op, each operand double or int32, broadcast or
     not. Flag bits 0-2 broadcast, 3-5 int32. */
  for (op = RSIMD_EW_NEG; op <= RSIMD_EW_ROUND; op++) {
    for (f = 0; f < 2; f++) {
      const void *x = f ? (const void *) ia : (const void *) fa;
      int nan_made = 0;
      reset_out();
      st = RSIMD_KERNEL(ew1_f64)(op, x, n, f ? RSIMD_EW_I32(0) : 0, fout, &o);
      for (j = 0; j < n; j++) {
        double a = ew_get(x, f, 0, j);
        if (op == RSIMD_EW_SQRT && a < 0) nan_made = 1;
        fref[j] = ref_ew1(op, a);
      }
      snprintf(what, sizeof what, "ew1_f64 op %d i32 %d", op, f);
      check_f64_kind(what, n, fout, fref);
      check_int(what, n, 0, st, nan_made ? RSIMD_EW_NAN_PRODUCED : 0);
    }
  }
  for (op = RSIMD_EW_ADD; op <= RSIMD_EW_COPYSIGN; op++) {
    for (f = 0; f < 64; f++) {
      int sc = f & 7, ti = f >> 3;
      const double *fx = op == RSIMD_EW_MOD || op == RSIMD_EW_IDIV ? fmodx : fa;
      const double *fy = op == RSIMD_EW_MOD || op == RSIMD_EW_IDIV ? fmody : fb;
      const void *x = ti & 1 ? (const void *) ia : (const void *) fx;
      const void *y = ti & 2 ? (const void *) ib : (const void *) fy;
      if ((sc & 3) == 3 || (sc & 4) || (ti & 4)) continue;
      reset_out();
      RSIMD_KERNEL(ew2_f64)(op, x, y, n, f, fout, &o);
      for (j = 0; j < n; j++) {
        fref[j] = ref_ew2(op, ew_get(x, ti & 1, sc & 1, j), ew_get(y, ti & 2, sc & 2, j));
      }
      snprintf(what, sizeof what, "ew2_f64 op %d flags %d", op, f);
      check_f64_kind(what, n, fout, fref);
    }
  }
  for (op = RSIMD_EW_FMA; op <= RSIMD_EW_MUL_ADD_APPROX; op++) {
    for (f = 0; f < 64; f++) {
      int sc = f & 7, ti = f >> 3, lohi = 0;
      const void *x = ti & 1 ? (const void *) ia : (const void *) fa;
      const void *y = ti & 2 ? (const void *) ib : (const void *) fb;
      const void *z = ti & 4 ? (const void *) ismall : (const void *) fc;
      if (sc == 7) continue;
      reset_out();
      st = RSIMD_KERNEL(ew3_f64)(op, x, y, z, n, f, fout, &o);
      for (j = 0; j < n; j++) {
        double a = ew_get(x, ti & 1, sc & 1, j), b = ew_get(y, ti & 2, sc & 2, j),
               c = ew_get(z, ti & 4, sc & 4, j);
        if (op == RSIMD_EW_CLAMP && b > c) lohi = 1;
        fref[j] = ref_ew3(op, a, b, c);
      }
      want_st = lohi ? RSIMD_EW_LO_GT_HI : 0;
      snprintf(what, sizeof what, "ew3_f64 op %d flags %d", op, f);
      check_f64_kind(what, n, fout, fref);
      check_int(what, n, 0, st, want_st);
    }
  }
#else
  (void) want_st;
#endif
}

/* ---- Predicates, comparisons, logic, bitwise ops and conversions ---------- */

#define SENTINEL_U8 0xA5
static double fcv[N + 1];                /* boundaries of the conversions */
static int32_t icv[N + 1];               /* around 0..255, with NA */
static uint8_t ua[N + 1], ub[N + 1], uout[N + 1], uref[N + 1];
static int32_t pa[N + 1];                /* any/all decided by the last element */
static int32_t pk[N + 1];                /* number classes: powers of two, odd, even */
#ifndef RSIMD_NO_F64_SIMD
static double pf[N + 1];
static double pcl[N + 1];                /* number classes: subnormal, whole, 2^k ... */
#endif

static void init_logical_inputs(void) {
  const double cvals[] = {
    2147483647.9, 2147483648.0, -2147483648.0, -2147483648.5, -2147483647.5, 2147483647.0,
    -2147483647.0, 255.9, 256.0, 255.0, -0.5, -1.0, -0.0, 0.0, 0.999, 4294967301.0,
    -4294967301.0, 1e19, -1e19, 9007199254740993.0, 1e300, HUGE_VAL, -HUGE_VAL,
    from_bits(UINT64_C(0x7FF8000000000000)), from_bits(UINT64_C(0x7FF00000000007A2)),
    4.9e-324, -4.9e-324, 65535.5, -255.5, 2.5, 3e9, 2147483904.0, -2147483904.0, 6442450944.0};
  const int nc = (int) (sizeof cvals / sizeof cvals[0]);
  const int32_t kvals[] = {0, 1, -1, 2, -2, RSIMD_NA_I32, INT32_MAX, -INT32_MAX, 1 << 30, 6, 7,
                           3, 4, 8, 12, 16, -16, 1 << 20, (1 << 20) + 1, 64};
  const int nk = (int) (sizeof kvals / sizeof kvals[0]);
#ifndef RSIMD_NO_F64_SIMD
  const double clvals[] = {1, 2, 0.5, 3, -4, 4.9e-324, DBL_MIN, DBL_MIN / 2, 0.0, -0.0,
                           HUGE_VAL, -HUGE_VAL, rsimd_na_real(), from_bits(UINT64_C(0x7FF8000000000000)),
                           0x1p1023, 0x1p53, 0x1p53 + 2, 1.5, -4.9e-324, DBL_MAX, -3, 7,
                           0x1p52 + 1, 0x1p52 + 0.5, -0x1p-1022, 0x1p-1060, 0x1.8p-1060, 1024,
                           -1024, 0.25, 6, 2.5, 0x1.fffffffffffffp52, 1e300, 0x1p-1074 * 3};
  const int ncl = (int) (sizeof clvals / sizeof clvals[0]);
#endif
  int i;
  for (i = 0; i < N; i++) {
    if (i < 2 * nc) {
      fcv[i] = cvals[(i * 7) % nc];
    } else {
      double m = (double) (int64_t) (next_rand() >> 11) / 9007199254740992.0 - 0.5;
      fcv[i] = ldexp(m, (int) (next_rand() % 72));
    }
    icv[i] = i % 11 == 5 ? RSIMD_NA_I32 : (int32_t) (next_rand() % 400) - 70;
    pk[i] = kvals[(i * 5) % nk];
#ifndef RSIMD_NO_F64_SIMD
    pcl[i] = clvals[(i * 7) % ncl];
#endif
    ua[i] = (uint8_t) (i * 37 + 11);
    ub[i] = i % 9 == 0 ? 0 : (uint8_t) next_rand();
  }
  fcv[N] = 0;
}

static void check_u8(const char *what, ptrdiff_t n, const uint8_t *got, const uint8_t *want) {
  ptrdiff_t i;
  char buf[96];
  for (i = 0; i < n; i++) {
    n_checks++;
    if (got[i] != want[i]) {
      snprintf(buf, sizeof buf, "got %d want %d", got[i], want[i]);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (got[n] != SENTINEL_U8) fail(what, n, n, "wrote past the end");
}

static void reset_out_u8(void) {
  int i;
  reset_out();
  for (i = 0; i <= N; i++) uout[i] = SENTINEL_U8;
}

/* The answer of a predicate in `mode` from the elementwise results. */
static int pred_fold(const int32_t *v, ptrdiff_t n, int mode) {
  ptrdiff_t j;
  for (j = 0; j < n; j++) {
    if (mode == RSIMD_PRED_ANY && v[j]) return 1;
    if (mode == RSIMD_PRED_ALL && !v[j]) return 0;
  }
  return mode == RSIMD_PRED_ALL;
}

static void test_logical(ptrdiff_t n) {
  const rsimd_opts o = {0, 1, RSIMD_PREC_FAST}, onc = {0, 0, RSIMD_PREC_FAST},
                   orm = {1, 1, RSIMD_PREC_FAST};
  const rsimd_opts *opts[3] = {&o, &onc, &orm};
  int op, f, mode, k, check, got, st, want_st;
  ptrdiff_t j;
  char what[96];
  if (n == 0) return; /* kernels are never called for no elements */

  /* Predicates, on inputs whose any/all is decided by the last element
     too (one NA or zero at the end of non-missing, non-zero values). */
  for (j = 0; j < n; j++) pa[j] = j + 1 == n ? (n % 2 ? RSIMD_NA_I32 : 0) : 3;
  for (op = RSIMD_PRED_NA; op <= RSIMD_PRED_POW2; op++) {
    for (f = 0; f < 3; f++) {
      const int32_t *x = f == 0 ? ia : f == 1 ? pa : pk;
      for (j = 0; j < n; j++) iref[j] = rsimd_pred_i32_1(op, x[j]);
      for (mode = RSIMD_PRED_ELT; mode <= RSIMD_PRED_ALL; mode++) {
        reset_out();
        got = RSIMD_KERNEL(pred_i32)(op, x, n, mode, iout);
        snprintf(what, sizeof what, "pred_i32 op %d mode %d input %d", op, mode, f);
        if (mode == RSIMD_PRED_ELT) check_i32(what, n, iout, iref);
        else check_int(what, n, 0, got, pred_fold(iref, n, mode));
      }
    }
  }
#ifndef RSIMD_SKIP_pred_f64
  for (j = 0; j < n; j++) {
    pf[j] = j + 1 == n ? (n % 3 == 0 ? rsimd_na_real() : n % 3 == 1 ? -0.0 : HUGE_VAL) : 2.5;
  }
  for (op = RSIMD_PRED_NA; op <= RSIMD_PRED_POW2; op++) {
    for (f = 0; f < 4; f++) {
      const double *x = f == 0 ? fa : f == 1 ? fcv : f == 2 ? pf : pcl;
      for (j = 0; j < n; j++) iref[j] = rsimd_pred_f64_1(op, x[j]);
      for (mode = RSIMD_PRED_ELT; mode <= RSIMD_PRED_ALL; mode++) {
        reset_out();
        got = RSIMD_KERNEL(pred_f64)(op, x, n, mode, iout);
        snprintf(what, sizeof what, "pred_f64 op %d mode %d input %d", op, mode, f);
        if (mode == RSIMD_PRED_ELT) check_i32(what, n, iout, iref);
        else check_int(what, n, 0, got, pred_fold(iref, n, mode));
      }
    }
  }
#endif

  /* Comparisons and three-valued logic, every operand broadcast or not;
     cmp_f64 and logic_f64 also with int32 operands (flags bits 3-4). */
  for (op = RSIMD_CMP_EQ; op <= RSIMD_CMP_GE; op++) {
    for (f = 0; f < 3; f++) {
      reset_out();
      RSIMD_KERNEL(cmp_i32)(op, ia, ib, n, f, iout);
      for (j = 0; j < n; j++) {
        int32_t a = ia[f & 1 ? 0 : j], b = ib[f & 2 ? 0 : j];
        iref[j] = a == RSIMD_NA_I32 || b == RSIMD_NA_I32 ? RSIMD_NA_I32 : rsimd_cmp_1(op, a, b);
      }
      snprintf(what, sizeof what, "cmp_i32 op %d flags %d", op, f);
      check_i32(what, n, iout, iref);
      reset_out();
      RSIMD_KERNEL(cmp_u8)(op, ua, ub, n, f, iout);
      for (j = 0; j < n; j++) iref[j] = rsimd_cmp_1(op, ua[f & 1 ? 0 : j], ub[f & 2 ? 0 : j]);
      snprintf(what, sizeof what, "cmp_u8 op %d flags %d", op, f);
      check_i32(what, n, iout, iref);
    }
#ifndef RSIMD_SKIP_cmp_f64
    for (f = 0; f < 32; f++) {
      int sc = f & 3, ti = f >> 3;
      const void *x = ti & 1 ? (const void *) ia : (const void *) fa;
      const void *y = ti & 2 ? (const void *) ib : (const void *) fb;
      if (sc == 3 || (f & 4)) continue;
      reset_out();
      RSIMD_KERNEL(cmp_f64)(op, x, y, n, f, iout);
      for (j = 0; j < n; j++) {
        double a = ew_get(x, ti & 1, sc & 1, j), b = ew_get(y, ti & 2, sc & 2, j);
        iref[j] = isnan(a) || isnan(b) ? RSIMD_NA_I32 : rsimd_cmp_1(op, a, b);
      }
      snprintf(what, sizeof what, "cmp_f64 op %d flags %d", op, f);
      check_i32(what, n, iout, iref);
    }
#endif
  }
  for (op = RSIMD_LOGIC_AND; op <= RSIMD_LOGIC_NOT; op++) {
    for (f = 0; f < 3; f++) {
      const int32_t *y = op == RSIMD_LOGIC_NOT ? NULL : ib;
      if (op == RSIMD_LOGIC_NOT && f & 2) continue;
      reset_out();
      RSIMD_KERNEL(logic_i32)(op, ia, y, n, f, iout);
      for (j = 0; j < n; j++) {
        int a = rsimd_lgl_of_i32(ia[f & 1 ? 0 : j]);
        int b = y == NULL ? 0 : rsimd_lgl_of_i32(ib[f & 2 ? 0 : j]);
        iref[j] = rsimd_logic_1(op, a, b);
      }
      snprintf(what, sizeof what, "logic_i32 op %d flags %d", op, f);
      check_i32(what, n, iout, iref);
    }
#ifndef RSIMD_SKIP_logic_f64
    for (f = 0; f < 32; f++) {
      int sc = f & 3, ti = f >> 3;
      const void *x = ti & 1 ? (const void *) ia : (const void *) fa;
      const void *y = op == RSIMD_LOGIC_NOT ? NULL : ti & 2 ? (const void *) ib : (const void *) fb;
      if (sc == 3 || (f & 4) || (op == RSIMD_LOGIC_NOT && (sc & 2 || ti & 2))) continue;
      reset_out();
      RSIMD_KERNEL(logic_f64)(op, x, y, n, f, iout);
      for (j = 0; j < n; j++) {
        int a = rsimd_lgl_of_f64(ew_get(x, ti & 1, sc & 1, j));
        int b = y == NULL ? 0 : rsimd_lgl_of_f64(ew_get(y, ti & 2, sc & 2, j));
        iref[j] = rsimd_logic_1(op, a, b);
      }
      snprintf(what, sizeof what, "logic_f64 op %d flags %d", op, f);
      check_i32(what, n, iout, iref);
    }
#endif
  }

  /* Bitwise ops: binary ops with every broadcast, the others with every
     count, with and without NA checks. */
  for (check = 0; check < 2; check++) {
    for (op = RSIMD_BIT_AND; op <= RSIMD_BIT_TZCNT; op++) {
      const int binary = op <= RSIMD_BIT_XOR;
      const int kmax = op >= RSIMD_BIT_SHL && op <= RSIMD_BIT_ROTR ? 31 : 0;
      for (f = 0; f < (binary ? 3 : 1); f++) {
        for (k = 0; k <= kmax; k++) {
          reset_out();
          RSIMD_KERNEL(bit_i32)(op, ia, binary ? ib : NULL, n, f, k, iout, check ? &o : &onc);
          for (j = 0; j < n; j++) {
            int32_t a = ia[f & 1 ? 0 : j], b = binary ? ib[f & 2 ? 0 : j] : 0;
            int na = check && (a == RSIMD_NA_I32 || b == RSIMD_NA_I32);
            iref[j] = na ? RSIMD_NA_I32 : rsimd_bit_i32_1(op, a, b, k);
          }
          snprintf(what, sizeof what, "bit_i32 op %d flags %d k %d check %d", op, f, k, check);
          check_i32(what, n, iout, iref);
        }
      }
    }
  }
  for (op = RSIMD_BIT_AND; op <= RSIMD_BIT_TZCNT; op++) {
    const int binary = op <= RSIMD_BIT_XOR, counts = rsimd_bit_counts(op);
    const int kmax = op == RSIMD_BIT_SHL || op == RSIMD_BIT_SHR ? 8
                   : op == RSIMD_BIT_ROTL || op == RSIMD_BIT_ROTR ? 7 : 0;
    if (op == RSIMD_BIT_SAR) continue;
    for (f = 0; f < (binary ? 3 : 1); f++) {
      for (k = 0; k <= kmax; k++) {
        reset_out_u8();
        RSIMD_KERNEL(bit_u8)(op, ua, binary ? ub : NULL, n, f, k,
                             counts ? (void *) iout : (void *) uout);
        for (j = 0; j < n; j++) {
          int32_t r = rsimd_bit_u8_1(op, ua[f & 1 ? 0 : j], binary ? ub[f & 2 ? 0 : j] : 0, k);
          if (counts) iref[j] = r;
          else uref[j] = (uint8_t) r;
        }
        snprintf(what, sizeof what, "bit_u8 op %d flags %d k %d", op, f, k);
        if (counts) check_i32(what, n, iout, iref);
        else check_u8(what, n, uout, uref);
      }
    }
  }
  for (f = 0; f < 3; f++) {
    rsimd_reduce_result r, want;
    memset(&r, 0, sizeof r);
    memset(&want, 0, sizeof want);
    RSIMD_KERNEL(popcnt_sum_i32)(ia, n, &r, opts[f]);
    for (j = 0; j < n; j++) {
      if (opts[f]->na_check && ia[j] == RSIMD_NA_I32) {
        want.saw_na = 1;
        if (!opts[f]->na_rm) break;
      } else {
        want.i64 += rsimd_popcnt_1((uint32_t) ia[j]);
      }
    }
    snprintf(what, sizeof what, "popcnt_sum_i32 opts %d", f);
    check_int(what, n, 0, r.saw_na, want.saw_na);
    /* Without na.rm the total is discarded once an NA is seen. */
    if (!want.saw_na || opts[f]->na_rm) check_int(what, n, 0, (long) r.i64, (long) want.i64);
  }
  {
    rsimd_reduce_result r;
    long want = 0;
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(popcnt_sum_u8)(ua, n, &r);
    for (j = 0; j < n; j++) want += rsimd_popcnt_1(ua[j]);
    check_int("popcnt_sum_u8", n, 0, (long) r.i64, want);
  }

  /* Conversions, every op and mode, with their status bits. */
#ifndef RSIMD_SKIP_convert
  for (op = RSIMD_CVT_F64_I32; op <= RSIMD_CVT_U8_LGL; op++) {
    for (mode = RSIMD_CVT_CHECKED; mode <= RSIMD_CVT_TRUNCATING; mode++) {
      for (f = 0; f < 2; f++) {
        const int from_f64 = op <= RSIMD_CVT_F64_LGL;
        const int from_i32 = op >= RSIMD_CVT_I32_F64 && op <= RSIMD_CVT_I32_LGL;
        const void *x = from_f64 ? (const void *) (f ? fa : fcv)
                      : from_i32 ? (const void *) (f ? ia : icv) : (const void *) (f ? ub : ua);
        void *out = op == RSIMD_CVT_I32_F64 || op == RSIMD_CVT_U8_F64  ? (void *) fout
                    : op == RSIMD_CVT_F64_U8 || op == RSIMD_CVT_I32_U8 ? (void *) uout
                                                                       : (void *) iout;
        static double fwant[N + 1];
        reset_out_u8();
        st = RSIMD_KERNEL(convert)(op, mode, x, n, out);
        want_st = 0;
        for (j = 0; j < n; j++) {
          void *w = out == (void *) fout   ? (void *) fwant
                    : out == (void *) uout ? (void *) uref
                                           : (void *) iref;
          want_st |= rsimd_cvt_1(op, mode, x, j, w);
        }
        snprintf(what, sizeof what, "convert op %d mode %d input %d", op, mode, f);
        if (out == (void *) fout) {
          fwant[n] = SENTINEL_F64;
          check_f64_kind(what, n, fout, fwant);
          check_int(what, n, n, bits(fout[n]) == bits(SENTINEL_F64), 1);
        } else if (out == (void *) uout) {
          check_u8(what, n, uout, uref);
        } else {
          check_i32(what, n, iout, iref);
        }
        check_int(what, n, 0, st, want_st);
      }
    }
  }
#else
  (void) st;
  (void) want_st;
#endif
}

/* Every byte value through every byte op and count. */
static void test_bytes_exhaustive(void) {
  static uint8_t x[257], y[257], out8[257];
  static int32_t out32[257];
  int op, k, f;
  ptrdiff_t j;
  char what[96];
  for (j = 0; j < 256; j++) {
    x[j] = (uint8_t) j;
    y[j] = (uint8_t) (255 - j);
  }
  for (op = RSIMD_BIT_AND; op <= RSIMD_BIT_TZCNT; op++) {
    const int binary = op <= RSIMD_BIT_XOR, counts = rsimd_bit_counts(op);
    const int kmax = op == RSIMD_BIT_SHL || op == RSIMD_BIT_SHR ? 8
                   : op == RSIMD_BIT_ROTL || op == RSIMD_BIT_ROTR ? 7 : 0;
    if (op == RSIMD_BIT_SAR) continue;
    for (f = 0; f < (binary ? 2 : 1); f++) {
      for (k = 0; k <= kmax; k++) {
        out8[256] = SENTINEL_U8;
        out32[256] = SENTINEL_I32;
        RSIMD_KERNEL(bit_u8)(op, x, binary ? y : NULL, 256, f ? RSIMD_EW_SCALAR(1) : 0, k,
                             counts ? (void *) out32 : (void *) out8);
        for (j = 0; j < 256; j++) {
          int32_t r = rsimd_bit_u8_1(op, x[j], binary ? y[f ? 0 : j] : 0, k);
          snprintf(what, sizeof what, "bit_u8 exhaustive op %d flags %d k %d", op, f, k);
          check_int(what, 256, j, counts ? out32[j] : out8[j], r);
        }
      }
    }
  }
}

/* ---- Complex kernels ------------------------------------------------------ */

#ifndef RSIMD_NO_F64_SIMD
static Rcomplex cz[N + 1], czout[N + 1];
static double cre[N + 1], cim[N + 1];

/* The layer's deinterleave on two full vectors. */
static void test_uzp(void) {
  ptrdiff_t W = RSIMD_LANES_64, j;
  double a[2 * RSIMD_MAX_LANES_64], e[RSIMD_MAX_LANES_64 + 1], o[RSIMD_MAX_LANES_64 + 1];
  for (j = 0; j < 2 * W; j++) a[j] = (double) j;
  e[W] = o[W] = SENTINEL_F64;
  rsimd_vf64_storeu(e, rsimd_vf64_uzp_even(rsimd_vf64_loadu(a), rsimd_vf64_loadu(a + W)));
  rsimd_vf64_storeu(o, rsimd_vf64_uzp_odd(rsimd_vf64_loadu(a), rsimd_vf64_loadu(a + W)));
  for (j = 0; j < W; j++) {
    check_int("uzp_even", W, j, (long) e[j], (long) (2 * j));
    check_int("uzp_odd", W, j, (long) o[j], (long) (2 * j + 1));
  }
}

/* Complex numbers from the reduction inputs: rd[i] and another element of
   rd, so that NA and NaN turn up in either part or both. */
static void fill_complex(void) {
  int i;
  for (i = 0; i < N; i++) {
    cz[i].r = rd[i];
    cz[i].i = rd[(i * 7 + 3) % N];
  }
}

static int same_or_nan(double a, double b) { return bits(a) == bits(b) || (isnan(a) && isnan(b)); }

static void test_complex(ptrdiff_t n) {
  static const int ops[] = {RSIMD_PRED_NA, RSIMD_PRED_NAN, RSIMD_PRED_FINITE, RSIMD_PRED_INFINITE};
  int variant, k, mode, narm, check, prec, way;
  char what[64];
  ptrdiff_t i;
  for (variant = 0; variant < 4; variant++) {
    fill_reduce(variant);
    fill_complex();

    /* conj: bit-exact, the imaginary sign flipped */
    for (i = 0; i <= n; i++) czout[i].r = czout[i].i = SENTINEL_F64;
    RSIMD_KERNEL(conj_c128)(cz, n, czout);
    for (i = 0; i < n; i++) {
      n_checks++;
      if (bits(czout[i].r) != bits(cz[i].r) ||
          bits(czout[i].i) != (bits(cz[i].i) ^ UINT64_C(0x8000000000000000))) {
        fail("conj_c128", n, i, "mismatch");
      }
    }
    n_checks++;
    if (bits(czout[n].r) != bits(SENTINEL_F64)) fail("conj_c128", n, n, "wrote past the end");

    /* real and imaginary parts */
    for (k = 0; k < 2; k++) {
      for (i = 0; i < n; i++) fref[i] = k ? cz[i].i : cz[i].r;
      fout[n] = SENTINEL_F64;
      RSIMD_KERNEL(part_c128)(k, cz, n, fout);
      check_f64(k ? "part_c128 im" : "part_c128 re", n, fout, fref, 1);
    }

    /* predicates */
    for (k = 0; k < 4; k++) {
      int want_any = 0, want_all = 1;
      for (i = 0; i < n; i++) {
        iref[i] = rsimd_pred_c128_1(ops[k], cz[i]);
        want_any |= iref[i];
        want_all &= iref[i];
      }
      iout[n] = SENTINEL_I32;
      snprintf(what, sizeof what, "pred_c128 op%d v%d", ops[k], variant);
      RSIMD_KERNEL(pred_c128)(ops[k], cz, n, RSIMD_PRED_ELT, iout);
      check_i32(what, n, iout, iref);
      check_int(what, n, -1, RSIMD_KERNEL(pred_c128)(ops[k], cz, n, RSIMD_PRED_ANY, NULL), want_any);
      check_int(what, n, -2, RSIMD_KERNEL(pred_c128)(ops[k], cz, n, RSIMD_PRED_ALL, NULL), want_all);
    }

    /* missing-value scans */
    for (mode = 0; mode < 3; mode++) {
      rsimd_reduce_result r;
      long cnt = 0;
      int any = 0;
      for (i = 0; i < n; i++) {
        if (isnan(cz[i].r) || isnan(cz[i].i)) iref[cnt++] = (int) (i + 1);
      }
      any = cnt > 0;
      reduce_init(&r, 0.0);
      r.i64 = 0;
      snprintf(what, sizeof what, "na_c128 mode%d v%d", mode, variant);
      RSIMD_KERNEL(na_c128)(cz, n, mode, 0, iout, &r);
      if (mode == RSIMD_NAMODE_ANY) {
        check_int(what, n, 0, r.saw_na, any);
      } else if (mode == RSIMD_NAMODE_COUNT) {
        check_int(what, n, 0, (long) r.i64, cnt);
      } else {
        check_int(what, n, 0, (long) r.i64, cnt);
        for (i = 0; i < cnt; i++) check_int(what, n, i + 1, iout[i], iref[i]);
      }
    }

    /* sum: each part as the scalar fold of that part (the values are exact
       in any order, so only Inf - Inf can differ, as NaN) */
    for (narm = 0; narm < 2; narm++) {
      for (check = 0; check < 2; check++) {
        for (prec = 0; prec < 3; prec++) {
          rsimd_opts o = {narm, check, prec}, po = {0, check, prec};
          rsimd_reduce_result r[2], want[2];
          ptrdiff_t removed = 0;
          if (!check && !narm && (variant == 1 || variant == 3)) continue;
          for (i = 0; i < n; i++) {
            cre[i] = cz[i].r;
            cim[i] = cz[i].i;
            if (narm && (isnan(cre[i]) || isnan(cim[i]))) {
              cre[i] = cim[i] = 0.0;
              removed++;
            }
          }
          if (narm) po.na_check = 0;
          reduce_init(&want[0], 0.0);
          reduce_init(&want[1], 0.0);
          rsimd_fold_f64(cre, NULL, n, RSIMD_TERM_X, &want[0], &po);
          rsimd_fold_f64(cim, NULL, n, RSIMD_TERM_X, &want[1], &po);
          want[0].count -= removed;
          want[1].count -= removed;
#define REF_CHECK                                                                         \
  for (k = 0; k < 2; k++) {                                                               \
    snprintf(what, sizeof what, "sum_c128 part%d v%d narm%d chk%d prec%d way%d", k, variant, \
             narm, check, prec, way);                                                     \
    check_flags(what, n, &r[k], &want[k]);                                                \
    n_checks++;                                                                           \
    if (!same_or_nan(rsimd_reduce_value(&r[k], prec), rsimd_reduce_value(&want[k], prec))) { \
      fail(what, n, 0, "sum differs");                                                    \
    }                                                                                     \
  }
          /* Pairwise chunks must start at a multiple of the leaf size. */
          for (way = 0; way < (prec == RSIMD_PREC_PAIRWISE ? 1 : 2); way++) {
            ptrdiff_t kk = way ? n / 3 : n, off, len;
            reduce_init(&r[0], 0.0);
            reduce_init(&r[1], 0.0);
            for (off = 0; off < n; off += len) {
              len = off < kk ? kk - off : n - off;
              RSIMD_KERNEL(sum_c128)(cz + off, len, r, &o);
            }
            REF_CHECK;
          }
#undef REF_CHECK
        }
      }
    }
  }
}
#endif

#ifndef RSIMD_NO_F64_SIMD
/* The layer's interleave on two full vectors. */
static void test_zip(void) {
  ptrdiff_t W = RSIMD_LANES_64, j;
  double a[RSIMD_MAX_LANES_64], b[RSIMD_MAX_LANES_64], lo[RSIMD_MAX_LANES_64 + 1],
    hi[RSIMD_MAX_LANES_64 + 1];
  for (j = 0; j < W; j++) {
    a[j] = (double) (2 * j);
    b[j] = (double) (2 * j + 1);
  }
  lo[W] = hi[W] = SENTINEL_F64;
  rsimd_vf64_storeu(lo, rsimd_vf64_zip_lo(rsimd_vf64_loadu(a), rsimd_vf64_loadu(b)));
  rsimd_vf64_storeu(hi, rsimd_vf64_zip_hi(rsimd_vf64_loadu(a), rsimd_vf64_loadu(b)));
  for (j = 0; j < W; j++) {
    check_int("zip_lo", W, j, (long) lo[j], (long) j);
    check_int("zip_hi", W, j, (long) hi[j], (long) (W + j));
  }
}

/* Base R's operators, as the runtime's __muldc3 and __divdc3 compute them
   (what C99's operators call where the inline product has two NaN parts,
   and for every quotient). */
extern double _Complex __muldc3(double, double, double, double);
extern double _Complex __divdc3(double, double, double, double);
static void h_mul1(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  double _Complex z = __muldc3(x->r, x->i, y->r, y->i);
  memcpy(out, &z, sizeof *out);
}
static void h_div1(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  double _Complex z = __divdc3(x->r, x->i, y->r, y->i);
  memcpy(out, &z, sizeof *out);
}
static void h_cp1(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  Rcomplex t = *y;
  out->r = x->r * t.r - x->i * t.i;
  out->i = x->r * t.i + x->i * t.r;
}

static Rcomplex cy[N + 1], cref[N + 1];

static int same_c128(Rcomplex a, Rcomplex b) {
  return bits(a.r) == bits(b.r) && bits(a.i) == bits(b.i);
}

/* nan_blind: any NaN matches any NaN (NA included). */
static void check_c128(const char *what, ptrdiff_t n, const Rcomplex *got, const Rcomplex *want,
                       int nan_blind) {
  ptrdiff_t i;
  char buf[200];
  for (i = 0; i < n; i++) {
    n_checks++;
    if (!same_c128(got[i], want[i]) &&
        !(nan_blind && same_or_nan(got[i].r, want[i].r) && same_or_nan(got[i].i, want[i].i))) {
      snprintf(buf, sizeof buf, "got %a%+ai want %a%+ai", got[i].r, got[i].i, want[i].r,
               want[i].i);
      fail(what, n, i, buf);
    }
  }
  n_checks++;
  if (bits(got[n].r) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}

/* Operands for the complex arithmetic: 0 the reduction inputs (ties,
   zeros, infinities, NA and NaN), 1 finite numbers of moderate magnitude,
   2 numbers over the whole exponent range with zero parts (every scaling
   case of the division). */
static double carith_double(int lo, int hi) {
  uint64_t r = next_rand();
  double m = 1.0 + (double) (r >> 12) / 4503599627370496.0;
  return (r & 1 ? -1.0 : 1.0) * ldexp(m, lo + (int) (next_rand() % (uint64_t) (hi - lo + 1)));
}

static void fill_carith(int kind) {
  const int lo = kind == 1 ? -8 : -1074, hi = kind == 1 ? 8 : 1023;
  ptrdiff_t i;
  for (i = 0; i < N; i++) {
    if (kind == 0) {
      cy[i].r = rd[(i * 5 + 1) % N];
      cy[i].i = rd[(i * 3 + 2) % N];
      continue;
    }
    cz[i].r = carith_double(lo, hi);
    cz[i].i = carith_double(lo, hi);
    cy[i].r = carith_double(lo, hi);
    cy[i].i = carith_double(lo, hi);
    if (kind == 2 && i % 5 == 0) cz[i].i = 0.0;
    if (kind == 2 && i % 7 == 0) cy[i].r = 0.0;
    if (kind == 2 && i % 11 == 0) cy[i].i = -0.0;
  }
}

static const int cmul_variants[10][2] = {{1, 1}, {0, 0}, {2, 1}, {1, 2}, {2, 2},
                                         {1, 0}, {0, 1}, {2, 0}, {0, 2}, {-1, -1}};

/* The fast-mode product as the vector tiers compute it: W lane products of
   the elements of full vectors (a skipped element as 1), combined
   pairwise, then the tail in order, then into s->p. */
static void cprod_fast_ref(const Rcomplex *x, ptrdiff_t n, rsimd_cprod_state *s,
                           const rsimd_c128_arith *a, const rsimd_opts *o) {
  const ptrdiff_t W = RSIMD_LANES_64;
  Rcomplex lane[RSIMD_MAX_LANES_64];
  rsimd_cprod_state t;
  ptrdiff_t i, j, k;
  for (j = 0; j < W; j++) {
    lane[j].r = 1.0;
    lane[j].i = 0.0;
  }
  for (i = 0; i + W <= n; i += W) {
    for (j = 0; j < W; j++) {
      Rcomplex e = x[i + j];
      if (rsimd_cprod_skip_(e, s, o)) {
        e.r = 1.0;
        e.i = 0.0;
      }
      lane[j] = rsimd_cmul_1_(a, lane[j], e);
    }
  }
  for (k = 1; k < W; k *= 2) {
    for (j = 0; j + k < W; j += 2 * k) lane[j] = rsimd_cmul_1_(a, lane[j], lane[j + k]);
  }
  t = *s;
  t.p = lane[0];
  rsimd_cprod_seq_(x, i, n, &t, a, o);
  s->saw_na = t.saw_na;
  s->saw_nan = t.saw_nan;
  s->p = rsimd_cmul_1_(a, s->p, t.p);
}

static void test_carith(ptrdiff_t n) {
  rsimd_c128_arith a = {0, 0, 0, 0, 0, h_mul1, h_div1, h_cp1};
  char what[96];
  int kind, v, f, d, prec, check, narm;
  ptrdiff_t i;
  for (kind = 0; kind < 3; kind++) {
    fill_reduce(kind == 0 ? 1 : 0);
    fill_complex();
    fill_carith(kind);
    /* x * y and x / y, every variant and broadcast */
    for (v = 0; v < 10; v++) {
      a.mul_re = cmul_variants[v][0];
      a.mul_im = cmul_variants[v][1];
      a.div = v < 3 ? v : RSIMD_CDIV_FMA;
      for (f = 0; f < 3; f++) {
        const int flags = f == 1 ? RSIMD_EW_SCALAR(0) : f == 2 ? RSIMD_EW_SCALAR(1) : 0;
        for (d = 0; d < 2; d++) {
          if (d && v >= 3) continue;
          for (i = 0; i < n; i++) {
            Rcomplex x = cz[f == 1 ? 0 : i], y = cy[f == 2 ? 0 : i];
            cref[i] = d ? rsimd_cdiv_1_(&a, x, y) : rsimd_cmul_1_(&a, x, y);
          }
          for (i = 0; i <= n; i++) czout[i].r = czout[i].i = SENTINEL_F64;
          RSIMD_KERNEL(ew2_c128)(d ? RSIMD_EW_DIV : RSIMD_EW_MUL, cz, cy, n, flags, czout, &a);
          snprintf(what, sizeof what, "ew2_c128 %s kind%d variant%d flags%d", d ? "div" : "mul",
                   kind, v, flags);
          check_c128(what, n, czout, cref, 0);
        }
      }
#ifndef RSIMD_SKIP_formula_c128
      /* The formulas alone, against the scalar helpers where those use
         them (no NaN part). */
      if (v < 9) {
        RSIMD_KERNEL(formula_c128)(RSIMD_EW_MUL, a.mul_re, a.mul_im, cz, cy, n, czout);
        for (i = 0; i < n; i++) {
          Rcomplex r = rsimd_cmul_1_(&a, cz[i], cy[i]);
          n_checks++;
          if (!isnan(czout[i].r) && !isnan(czout[i].i) && !same_c128(czout[i], r)) {
            fail("formula_c128 mul", n, i, "differs");
          }
        }
      }
#endif
    }
    /* products */
    a.mul_re = a.mul_im = RSIMD_CMUL_FMA1;
    for (v = 0; v < 3; v++) {
      if (v == 1) a.mul_re = a.mul_im = RSIMD_CMUL_UNFUSED;
      if (v == 2) a.mul_re = a.mul_im = RSIMD_CMUL_SCALAR;
      for (prec = 0; prec < 3; prec++) {
        for (check = 0; check < 2; check++) {
          for (narm = 0; narm < 2; narm++) {
            rsimd_opts o = {narm, check, prec};
            rsimd_cprod_state got, want;
            ptrdiff_t k = n / 3;
            memset(&got, 0, sizeof got);
            got.p.r = 1.0;
            want = got;
            /* Two chunks, as the entry point would pass them. */
            RSIMD_KERNEL(prod_c128)(cz, k, &got, &a, &o);
            RSIMD_KERNEL(prod_c128)(cz + k, n - k, &got, &a, &o);
            if (prec == RSIMD_PREC_FAST && a.mul_re != RSIMD_CMUL_SCALAR &&
                RSIMD_LANES_64 > 1) {
              cprod_fast_ref(cz, k, &want, &a, &o);
              cprod_fast_ref(cz + k, n - k, &want, &a, &o);
            } else {
              rsimd_cprod_seq_(cz, 0, n, &want, &a, &o);
            }
            snprintf(what, sizeof what, "prod_c128 kind%d variant%d prec%d chk%d narm%d", kind,
                     v, prec, check, narm);
            n_checks++;
            if (!(o.na_check || o.na_rm) && (want.saw_nan || isnan(want.p.r))) continue;
            if (!same_c128(got.p, want.p) || !same_c128(got.lo, want.lo) ||
                got.saw_na != want.saw_na || got.saw_nan != want.saw_nan ||
                got.leaves != want.leaves) {
              fail(what, n, 0, "product or flags differ");
            }
          }
        }
      }
    }
    /* eq and ne */
    for (v = 0; v < 2; v++) {
      for (f = 0; f < 3; f++) {
        const int flags = f == 1 ? RSIMD_EW_SCALAR(0) : f == 2 ? RSIMD_EW_SCALAR(1) : 0;
        /* Finite operands: half of y equal to x. */
        if (kind != 0) {
          for (i = 0; i < n; i++) cy[i] = cz[(i / 2) * 2];
        }
        for (i = 0; i < n; i++) {
          Rcomplex x = cz[f == 1 ? 0 : i], y = cy[f == 2 ? 0 : i];
          if (isnan(x.r) || isnan(x.i) || isnan(y.r) || isnan(y.i)) iref[i] = RSIMD_NA_I32;
          else iref[i] = (x.r == y.r && x.i == y.i) == (v == 0);
        }
        iout[n] = SENTINEL_I32;
        RSIMD_KERNEL(cmp_c128)(v == 0 ? RSIMD_CMP_EQ : RSIMD_CMP_NE, cz, cy, n, flags, iout);
        snprintf(what, sizeof what, "cmp_c128 kind%d op%d flags%d", kind, v, flags);
        check_i32(what, n, iout, iref);
      }
    }
#ifndef RSIMD_SKIP_scan_c128
    /* cumsum and cumprod in two chunks */
    for (v = 0; v < 2; v++) {
      Rcomplex sg, sw;
      a.cp_re = v ? RSIMD_CMUL_FMA2 : RSIMD_CMUL_SCALAR;
      a.cp_im = v ? RSIMD_CMUL_FMA1 : RSIMD_CMUL_SCALAR;
      for (d = 0; d < 2; d++) {
        sg.r = sw.r = d ? 1.0 : 0.0;
        sg.i = sw.i = 0.0;
        for (i = 0; i < n; i++) {
          if (d == 0) {
            sw.r += cz[i].r;
            sw.i += cz[i].i;
          } else if (v == 0) {
            Rcomplex u = sw;
            h_cp1(cz + i, &u, &sw);
          } else {
            double re = rsimd_cmul_re_(a.cp_re, cz[i].r, cz[i].i, sw.r, sw.i);
            sw.i = rsimd_cmul_im_(a.cp_im, cz[i].r, cz[i].i, sw.r, sw.i);
            sw.r = re;
          }
          cref[i] = sw;
        }
        for (i = 0; i <= n; i++) czout[i].r = czout[i].i = SENTINEL_F64;
        RSIMD_KERNEL(scan_c128)(d, cz, n / 2, czout, &sg, &a);
        RSIMD_KERNEL(scan_c128)(d, cz + n / 2, n - n / 2, czout + n / 2, &sg, &a);
        /* NA versus NaN is the entry point's fix-up, and two inlined
           copies of the same step may propagate different payloads. */
        snprintf(what, sizeof what, "scan_c128 kind%d op%d variant%d", kind, d, v);
        check_c128(what, n, czout, cref, 1);
      }
    }
#endif
  }
}

#if !defined(RSIMD_SKIP_math1_c128) && defined(RSIMD_HAVE_SLEEF)
/* Mod and Arg: libm exactly where a part is not finite, else within
   SLEEF's bounds of the long double reference. */
static void test_cmath(ptrdiff_t n) {
  int kind, op, fast;
  static char what[64]; /* check_ulp keeps the name of the worst case */
  ptrdiff_t i;
  for (kind = 0; kind < 3; kind++) {
    fill_reduce(kind == 0 ? 1 : 0);
    fill_complex();
    fill_carith(kind);
    for (fast = 0; fast < 2; fast++) {
      for (op = 0; op < 2; op++) {
        for (i = 0; i < n; i++) {
          fa[i] = cz[i].r;
          fb[i] = cz[i].i;
          ldref[i] = op == RSIMD_CMATH_MOD ? hypotl(fa[i], fb[i]) : atan2l(fb[i], fa[i]);
        }
        fout[n] = SENTINEL_F64;
        RSIMD_KERNEL(math1_c128)(op | (fast ? RSIMD_MATH_FAST : 0), cz, n, fout);
        snprintf(what, sizeof what, "math1_c128 kind%d op%d fast%d", kind, op, fast);
        for (i = 0; i < n; i++) {
          if (!isfinite(cz[i].r) || !isfinite(cz[i].i)) {
            fref[i] = op == RSIMD_CMATH_MOD ? hypot(cz[i].r, cz[i].i) : atan2(cz[i].i, cz[i].r);
            n_checks++;
            if (bits(fout[i]) != bits(fref[i]) && !(isnan(fout[i]) && isnan(fref[i]))) {
              fail(what, n, i, "non-finite element differs from libm");
            }
            ldref[i] = (long double) fout[i];
          }
        }
        /* The tail is libm's, which is within 1 ULP. */
        check_ulp(what, n, fout, fast ? 3.5 : 1.0);
      }
    }
  }
}
#endif
#endif

/* ---- Complex elementary functions (cmath.inc.c) ------------------------- */

#ifndef RSIMD_NO_F64_SIMD
#include <complex.h>
#include <stdlib.h>

/* The base functions the kernels call for the elements they leave out: C99
   libm, and base R's algorithms for whole-number powers, log with a base
   and atan2 (with the multiply and divide of the current variant). */
/* The division variant matches this platform's __divdc3 (h_div1), which
   the base functions call: libgcc is built with fused multiply-adds on
   aarch64 and without on x86-64. */
#if defined(__aarch64__)
#define CM_DIV RSIMD_CDIV_FMA
#elif defined(__i386__)
/* i686 libgcc divides in x87 extended precision: no variant matches. */
#define CM_DIV RSIMD_CDIV_SCALAR
#else
#define CM_DIV RSIMD_CDIV_UNFUSED
#endif
static rsimd_c128_arith cm_arith = {RSIMD_CMUL_FMA1, RSIMD_CMUL_FMA1, CM_DIV, 0, 0,
                                    h_mul1, h_div1, h_cp1};

static double complex cm_c99(Rcomplex z) { return CMPLX(z.r, z.i); }
static void cm_out(double complex z, Rcomplex *out) {
  out->r = creal(z);
  out->i = cimag(z);
}

#define CM_H1(name, f)                                                                       \
  static void name(const Rcomplex *x, Rcomplex *out) { cm_out(f(cm_c99(*x)), out); }
/* Base R's own code (src/main/complex.c) for tan, asin, acos, atan and the
   inverse hyperbolic sine and tangent. */
static double complex z_tan(double complex z) {
  double y = cimag(z);
  double complex r = ctan(z);
  if (isfinite(y) && fabs(y) > 25.0) r = CMPLX(creal(r), y < 0 ? -1.0 : 1.0);
  return r;
}
static double complex z_asin(double complex z) {
  if (cimag(z) == 0 && fabs(creal(z)) > 1) {
    double alpha, t1, t2, x = creal(z), ri;
    t1 = 0.5 * fabs(x + 1);
    t2 = 0.5 * fabs(x - 1);
    alpha = t1 + t2;
    ri = log(alpha + sqrt(alpha * alpha - 1));
    if (x > 1) ri *= -1;
    return asin(t1 - t2) + ri * I;
  }
  return casin(z);
}
static double complex z_acos(double complex z) {
  if (cimag(z) == 0 && fabs(creal(z)) > 1) return M_PI_2 - z_asin(z);
  return cacos(z);
}
static double complex z_atan(double complex z) {
  if (creal(z) == 0 && fabs(cimag(z)) > 1) {
    double y = cimag(z), rr, ri;
    rr = (y > 0) ? M_PI_2 : -M_PI_2;
    ri = 0.25 * log(((y + 1) * (y + 1)) / ((y - 1) * (y - 1)));
    return rr + ri * I;
  }
  return catan(z);
}
static double complex z_asinh(double complex z) { return -I * z_asin(z * I); }
static double complex z_atanh(double complex z) { return -I * z_atan(z * I); }

CM_H1(h_csqrt, csqrt)
CM_H1(h_cexp, cexp)
CM_H1(h_clog, clog)
CM_H1(h_csin, csin)
CM_H1(h_ccos, ccos)
CM_H1(h_ctan, z_tan)
CM_H1(h_csinh, csinh)
CM_H1(h_ccosh, ccosh)
CM_H1(h_ctanh, ctanh)
CM_H1(h_casin, z_asin)
CM_H1(h_cacos, z_acos)
CM_H1(h_catan, z_atan)
CM_H1(h_casinh, z_asinh)
CM_H1(h_cacosh, cacosh)
CM_H1(h_catanh, z_atanh)

/* x^k by base R's binary powering, with the variant's multiply/divide. */
static Rcomplex h_ipow(Rcomplex x, int k) {
  Rcomplex z, one;
  int m = k < 0 ? -k : k;
  z.r = one.r = 1.0;
  z.i = one.i = 0.0;
  if (m == 0) return z;
  if (m == 1) z = x;
  else {
    while (m > 0) {
      if (m & 1) z = rsimd_cmul_1_(&cm_arith, z, x);
      if (m == 1) break;
      m >>= 1;
      x = rsimd_cmul_1_(&cm_arith, x, x);
    }
  }
  return k < 0 ? rsimd_cdiv_1_(&cm_arith, one, z) : z;
}

static int cm_whole(Rcomplex y) {
  return y.i == 0 && y.r == trunc(y.r) && fabs(y.r) <= 65536;
}

static void h_cpow(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  if (x->r == 0 && x->i == 0) {
    if (y->i == 0) {
      out->r = pow(0.0, y->r);
      out->i = 0.0;
    } else {
      out->r = out->i = NAN;
    }
  } else if (cm_whole(*y)) {
    *out = h_ipow(*x, (int) y->r);
  } else {
    cm_out(cpow(cm_c99(*x), cm_c99(*y)), out);
  }
}
static void h_clogb(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  Rcomplex a, b;
  cm_out(clog(cm_c99(*x)), &a);
  cm_out(clog(cm_c99(*y)), &b);
  h_div1(&a, &b, out);
}
static void h_catan2(const Rcomplex *x, const Rcomplex *y, Rcomplex *out) {
  Rcomplex q;
  double complex d;
  if (y->r == 0 && y->i == 0) {
    if (x->r == 0 && x->i == 0) {
      out->r = out->i = from_bits(UINT64_C(0x7FF00000000007A2));
      return;
    }
    out->r = isnan(x->r) ? x->r : x->r >= 0 ? M_PI_2 : -M_PI_2;
    out->i = 0.0;
    return;
  }
  h_div1(x, y, &q);
  d = catan(cm_c99(q));
  if (y->r < 0) d += M_PI;
  if (creal(d) > M_PI) d -= 2 * M_PI;
  cm_out(d, out);
}

static const rsimd_cmath_base cm_base = {
  {h_csqrt, h_cexp, h_clog, h_csin, h_ccos, h_ctan, h_csinh, h_ccosh, h_ctanh, h_casin, h_cacos,
   h_catan, h_casinh, h_cacosh, h_catanh},
  {h_cpow, h_clogb, h_catan2}};

static const char *const cm_names[RSIMD_CM_COUNT + RSIMD_CM2_COUNT] = {
  "sqrt", "exp", "log", "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin",
  "acos", "atan", "asinh", "acosh", "atanh", "pow", "logb", "atan2"};

/* The exact value (long double) of op at x (and y). */
static long double complex cm_ref(int op, Rcomplex x, Rcomplex y) {
  const long double complex z = CMPLXL(x.r, x.i), w = CMPLXL(y.r, y.i);
  long double complex d;
  switch (op) {
  case RSIMD_CM_SQRT: return csqrtl(z);
  case RSIMD_CM_EXP: return cexpl(z);
  case RSIMD_CM_LOG: return clogl(z);
  case RSIMD_CM_SIN: return csinl(z);
  case RSIMD_CM_COS: return ccosl(z);
  case RSIMD_CM_TAN: return ctanl(z);
  case RSIMD_CM_SINH: return csinhl(z);
  case RSIMD_CM_COSH: return ccoshl(z);
  case RSIMD_CM_TANH: return ctanhl(z);
  case RSIMD_CM_ASIN: return casinl(z);
  case RSIMD_CM_ACOS: return cacosl(z);
  case RSIMD_CM_ATAN: return catanl(z);
  case RSIMD_CM_ASINH: return casinhl(z);
  case RSIMD_CM_ACOSH: return cacoshl(z);
  case RSIMD_CM_ATANH: return catanhl(z);
  case RSIMD_CM_COUNT + RSIMD_CM2_POW: return cpowl(z, w);
  case RSIMD_CM_COUNT + RSIMD_CM2_LOGB: return clogl(z) / clogl(w);
  default: {
    /* Base R rounds the quotient before catan, so that is the input. */
    Rcomplex q = rsimd_cdiv_1_(&cm_arith, x, y);
    /* And adjusts by M_PI, deciding on the rounded double sum. */
    double t;
    d = catanl(CMPLXL(q.r, q.i));
    t = (double) creall(d);
    if (y.r < 0) {
      d += (long double) M_PI;
      t += M_PI;
    }
    if (t > M_PI) d -= 2 * (long double) M_PI;
    return d;
  }
  }
}

static double cm_ulp(double x) {
  int e;
  if (x == 0 || !isfinite(x)) return DBL_TRUE_MIN;
  frexp(x, &e);
  return fmax(ldexp(1.0, e - 53), DBL_TRUE_MIN);
}

/* The error of part g against the exact r, in ULPs of r or, for a part
   much smaller than the whole result, of |result| (mod), whichever is
   smaller. */
static double cm_err(double g, long double r, long double mod) {
  long double d;
  double ec, em;
  if (!isfinite(g) || !isfinite((double) r)) return (double) g == (double) r ? 0 : HUGE_VAL;
  d = fabsl((long double) g - r);
  if (d == 0) return 0;
  ec = (double) (d / cm_ulp((double) r));
  em = (double) (d / cm_ulp((double) mod));
  return ec < em ? ec : em;
}

#define CM_M 2048
static Rcomplex cmx[CM_M + 1], cmy[CM_M + 1], cmo[CM_M + 1];
/* Worst errors: [op][fast][0 rsimd's formulas, 1 the base functions]. */
static double cm_worst[RSIMD_CM_COUNT + RSIMD_CM2_COUNT][2][2];
static Rcomplex cm_worst_at[RSIMD_CM_COUNT + RSIMD_CM2_COUNT][2][2];
/* The worst pow error divided by 1 + |y log x|. */
static double cm_pow_scaled[2];

static double cm_rand_part(int kind) {
  uint64_t r = next_rand();
  double u = (double) (r >> 11) * 0x1.0p-53, s = r & 1 ? -1.0 : 1.0;
  switch (kind) {
  case 0: return s * u * 4; /* moderate */
  case 1: return carith_double(-1074, 1023); /* the whole range */
  case 2: return s * ldexp(1.0 + u, (int) (next_rand() % 60) - 30);
  default: return s * ldexp(1.0 + u, (int) (next_rand() % 20) - 10);
  }
}

/* Inputs of kind k: 0 moderate, 1 the whole exponent range, 2 modulus near
   1, 3 near the branch points +-1 and +-i, 4 on and near the axes (signed
   zeros and tiny parts), 5 a mix of specials (NA, NaN, infinities, zeros),
   6 magnitudes 2^-30..2^30, 7 near the overflow and underflow edges of the
   formulas (|parts| around 20, 25, 708, 710, 2^500). */
static void cm_fill(Rcomplex *v, ptrdiff_t n, int kind) {
  static const double edges[] = {20.0, 25.0, 708.0, 709.5, 710.0, 0x1p500, 0x1p-500, 0x1p-1000,
                                 0x1p1020, 1.0, 0.5};
  static const double sp[] = {0.0, -0.0, 1.0, -1.0, 2.0, -0.5, HUGE_VAL, -HUGE_VAL, NAN};
  const double na = from_bits(UINT64_C(0x7FF00000000007A2));
  ptrdiff_t i;
  for (i = 0; i < n; i++) {
    uint64_t r = next_rand();
    double u = (double) (r >> 11) * 0x1.0p-53, t = 2 * M_PI * u - M_PI;
    switch (kind) {
    case 0: v[i].r = cm_rand_part(0); v[i].i = cm_rand_part(0); break;
    case 1: v[i].r = cm_rand_part(1); v[i].i = cm_rand_part(1); break;
    case 2: {
      double m = 1.0 + ldexp(cm_rand_part(0), -(int) (next_rand() % 50));
      v[i].r = m * cos(t);
      v[i].i = m * sin(t);
      break;
    }
    case 3: {
      double e1 = ldexp(cm_rand_part(0), -(int) (next_rand() % 60)),
             e2 = ldexp(cm_rand_part(0), -(int) (next_rand() % 60));
      if (r & 2) {
        v[i].r = (r & 4 ? 1.0 : -1.0) + e1;
        v[i].i = r & 8 ? e2 : 0.0;
      } else {
        v[i].r = r & 8 ? e1 : 0.0;
        v[i].i = (r & 4 ? 1.0 : -1.0) + e2;
      }
      break;
    }
    case 4:
      v[i].r = r & 2 ? (r & 4 ? -0.0 : 0.0) : cm_rand_part(r & 8 ? 0 : 1);
      v[i].i = r & 2 ? cm_rand_part(r & 8 ? 0 : 1) : (r & 4 ? -0.0 : 0.0);
      if (r & 16) v[i].i = ldexp(cm_rand_part(0), -1000 - (int) (next_rand() % 70));
      break;
    case 5:
      v[i].r = r % 11 == 0 ? na : sp[(r >> 8) % 9];
      v[i].i = (r >> 4) % 11 == 0 ? na : sp[(r >> 16) % 9];
      if ((r >> 24) % 3 == 0) v[i].r = cm_rand_part(0);
      if ((r >> 28) % 3 == 0) v[i].i = cm_rand_part(0);
      break;
    case 6: v[i].r = cm_rand_part(2); v[i].i = cm_rand_part(2); break;
    default:
      v[i].r = (r & 2 ? -1 : 1) * edges[(r >> 4) % 11] * (1 + ldexp(u - 0.5, -(int) (next_rand() % 40)));
      v[i].i = (r & 64 ? -1 : 1) * edges[(r >> 12) % 11] * (1 + ldexp(u - 0.5, -(int) (next_rand() % 40)));
      if (r & 128) v[i].i = cm_rand_part(0);
      break;
    }
  }
}

/* Exponents for pow: kind 0 whole numbers in -70..70 (some +-65536,
   65537), 1 moderate complex, 2 moderate real, 3 the base kinds. */
static void cm_fill_pow(Rcomplex *v, ptrdiff_t n, int kind) {
  ptrdiff_t i;
  for (i = 0; i < n; i++) {
    uint64_t r = next_rand();
    if (kind == 0) {
      static const double big[] = {65536, -65536, 65537, 2, -1, 1, 0, -2};
      v[i].r = r % 9 == 0 ? big[(r >> 8) % 8] : (double) ((int) ((r >> 8) % 141) - 70);
      v[i].i = 0.0;
    } else if (kind == 1) {
      v[i].r = cm_rand_part(0);
      v[i].i = cm_rand_part(0);
    } else if (kind == 2) {
      v[i].r = cm_rand_part(0);
      v[i].i = r & 2 ? -0.0 : 0.0;
    } else {
      cm_fill(v + i, 1, (int) (r % 8));
    }
  }
}

static int cm_same(Rcomplex a, Rcomplex b) {
  return same_or_nan(a.r, b.r) && same_or_nan(a.i, b.i) &&
         rsimd_is_na_f64(a.r) == rsimd_is_na_f64(b.r) && rsimd_is_na_f64(a.i) == rsimd_is_na_f64(b.i);
}

/* Runs op (cm_names index) over n elements of cmx (and cmy, broadcast as
   flags say) and checks every element: what the base functions give, under
   base R's NA rules, wherever a part is not finite or the result is not
   finite or zero; within `bound` of the exact value elsewhere (bound < 0:
   bit-identical, for whole-number powers), and records the worst errors. */
static void cm_check(int op, int fast, ptrdiff_t n, int flags, double bound) {
  const int two = op >= RSIMD_CM_COUNT, o2 = two ? op - RSIMD_CM_COUNT : -1;
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  static char what[64];
  char buf[260];
  ptrdiff_t i;
  cmo[n].r = cmo[n].i = SENTINEL_F64;
  if (two) {
    RSIMD_KERNEL(cmath2_c128)(o2 | (fast ? RSIMD_MATH_FAST : 0), cmx, cmy, n, flags, cmo, &cm_base, &cm_arith);
  } else {
    RSIMD_KERNEL(cmath1_c128)(op | (fast ? RSIMD_MATH_FAST : 0), cmx, n, cmo, &cm_base);
  }
  snprintf(what, sizeof what, "cmath %s fast%d flags%d", cm_names[op], fast, flags);
  for (i = 0; i < n; i++) {
    const Rcomplex x = cmx[sx ? 0 : i], y = cmy[sy ? 0 : i];
    Rcomplex b;
    int is_base, special;
    n_checks++;
    if (two) rsimd_cmath2_1(&cm_base, o2, &x, &y, &b);
    else rsimd_cmath1_1(&cm_base, op, &x, &b);
    is_base = cm_same(cmo[i], b);
    special = !isfinite(x.r) || !isfinite(x.i) || (two && (!isfinite(y.r) || !isfinite(y.i))) ||
              !isfinite(b.r) || !isfinite(b.i);
    if (special || bound < 0 || (o2 == RSIMD_CM2_POW && cm_whole(y))) {
      if (!is_base) {
        snprintf(buf, sizeof buf, "input %a%+ai, %a%+ai: got %a%+ai want %a%+ai", x.r, x.i, y.r,
                 y.i, cmo[i].r, cmo[i].i, b.r, b.i);
        fail(what, n, i, buf);
      }
      continue;
    }
    {
      const long double complex r = cm_ref(op, x, y);
      const long double mod = cabsl(r);
      const double e = fmax(cm_err(cmo[i].r, creall(r), mod), cm_err(cmo[i].i, cimagl(r), mod));
      const double eb = fmax(cm_err(b.r, creall(r), mod), cm_err(b.i, cimagl(r), mod));
      double lim = bound;
      if (op == RSIMD_CM_COUNT + RSIMD_CM2_POW) {
        /* The error of exp(y log x) grows with |y log x|. */
        const long double complex w = CMPLXL(y.r, y.i) * clogl(CMPLXL(x.r, x.i));
        const double s = (double) (1 + cabsl(w));
        lim = bound * s;
        if (!is_base && e / s > cm_pow_scaled[fast]) cm_pow_scaled[fast] = e / s;
      }
      if (!isfinite((double) creall(r)) || !isfinite((double) cimagl(r))) continue;
      if (!is_base && e > cm_worst[op][fast][0]) {
        cm_worst[op][fast][0] = e;
        cm_worst_at[op][fast][0] = x;
      }
      /* Base R's values on the axes follow its own branch-cut code, not
         C99's. */
      if (eb > cm_worst[op][fast][1] && x.r != 0 && x.i != 0) {
        cm_worst[op][fast][1] = eb;
        cm_worst_at[op][fast][1] = x;
      }
      /* Zero parts must have the base functions' sign (a part that is a
         subnormal instead of zero is an error like any other), except in
         general powers, where it depends on how the C library rounds y log
         x. */
      if (!is_base && op != RSIMD_CM_COUNT + RSIMD_CM2_POW &&
          ((b.r == 0 && cmo[i].r == 0 && bits(cmo[i].r) != bits(b.r)) ||
           (b.i == 0 && cmo[i].i == 0 && bits(cmo[i].i) != bits(b.i)))) {
        snprintf(buf, sizeof buf, "input %a%+ai, %a%+ai: got %a%+ai want zero signs of %a%+ai",
                 x.r, x.i, y.r, y.i, cmo[i].r, cmo[i].i, b.r, b.i);
        fail(what, n, i, buf);
      } else if (!is_base && e > lim) {
        snprintf(buf, sizeof buf, "input %a%+ai, %a%+ai: got %a%+ai want %La%+Lai (%.2f ULP)",
                 x.r, x.i, y.r, y.i, cmo[i].r, cmo[i].i, creall(r), cimagl(r), e);
        fail(what, n, i, buf);
      }
    }
  }
  n_checks++;
  if (bits(cmo[n].r) != bits(SENTINEL_F64)) fail(what, n, n, "wrote past the end");
}

/* rsimd's bound per function in ULPs (accurate, fast), for the parts the
   formulas compute (the measured worst cases over every tier, rounded up);
   pow's is per 1 + |y log x|. The ?simd_exp complex section states them. */
static const double cm_bound[RSIMD_CM_COUNT + RSIMD_CM2_COUNT][2] = {
  {2, 3},     /* sqrt */
  {2.5, 3},   /* exp */
  {3, 5},     /* log */
  {3, 4},     /* sin */
  {3, 4},     /* cos */
  {4.5, 6},   /* tan */
  {3, 4},     /* sinh */
  {3, 4},     /* cosh */
  {5, 6},     /* tanh */
  {3.5, 4.5}, /* asin */
  {3.5, 4.5}, /* acos */
  {3, 3.5},   /* atan */
  {3.5, 4.5}, /* asinh */
  {3.5, 4.5}, /* acosh */
  {3, 3.5},   /* atanh */
  {3, 5},     /* pow, per 1 + |y log x| */
  {4.5, 6.5}, /* log with a base */
  {3.5, 4}};  /* atan2 */

/* Every function on every input kind at length n (tails), in both modes. */
static void test_cmath_fns(ptrdiff_t n) {
  int op, fast, kind, k2;
  for (kind = 0; kind < 8; kind++) {
    for (op = 0; op < RSIMD_CM_COUNT + RSIMD_CM2_COUNT; op++) {
      for (fast = 0; fast < 2; fast++) {
        if (op < RSIMD_CM_COUNT) {
          cm_fill(cmx, n, kind);
          cm_check(op, fast, n, 0, cm_bound[op][fast]);
        } else if (op == RSIMD_CM_COUNT + RSIMD_CM2_POW) {
          for (k2 = 0; k2 < 4; k2++) {
            cm_fill(cmx, n, kind);
            cm_fill_pow(cmy, n, k2);
            cm_check(op, fast, n, 0, cm_bound[op][fast]);
            /* A broadcast exponent, and a broadcast base. */
            cm_check(op, fast, n, RSIMD_EW_SCALAR(1), cm_bound[op][fast]);
            cm_check(op, fast, n, RSIMD_EW_SCALAR(0), cm_bound[op][fast]);
          }
        } else {
          cm_fill(cmx, n, kind);
          cm_fill(cmy, n, (kind + 3) % 8);
          cm_check(op, fast, n, 0, cm_bound[op][fast]);
          cm_check(op, fast, n, RSIMD_EW_SCALAR(1), cm_bound[op][fast]);
        }
      }
    }
  }
}

/* Many more inputs at full length, for the error statistics, then whole
   powers in the unfused multiply variant (x86-64 builds of R). */
static void test_cmath_sweep(void) {
  const char *e = getenv("RSIMD_CMATH_REPS");
  int rep, reps = e ? atoi(e) : 2;
  for (rep = 0; rep < reps; rep++) test_cmath_fns(CM_M);
  cm_arith.mul_re = cm_arith.mul_im = RSIMD_CMUL_UNFUSED;
  cm_arith.div = RSIMD_CDIV_UNFUSED;
  cm_fill(cmx, CM_M, 0);
  cm_fill_pow(cmy, CM_M, 0);
  cm_check(RSIMD_CM_COUNT + RSIMD_CM2_POW, 0, CM_M, 0, -1);
  cm_arith.mul_re = cm_arith.mul_im = RSIMD_CMUL_FMA1;
  cm_arith.div = CM_DIV;
}

static void cm_report(void) {
  int op, fast;
  if (!getenv("RSIMD_CMATH_REPORT")) return;
  printf("cmath worst errors in ULP (rsimd formulas | base functions), accurate / fast:\n");
  for (op = 0; op < RSIMD_CM_COUNT + RSIMD_CM2_COUNT; op++) {
    for (fast = 0; fast < 2; fast++) {
      printf("  %-6s %s rsimd %8.3f at %a%+ai | base %8.3f at %a%+ai\n", cm_names[op],
             fast ? "fast" : "acc ", cm_worst[op][fast][0], cm_worst_at[op][fast][0].r,
             cm_worst_at[op][fast][0].i, cm_worst[op][fast][1], cm_worst_at[op][fast][1].r,
             cm_worst_at[op][fast][1].i);
    }
  }
  printf("  pow per 1 + |y log x|: %.3f / %.3f\n", cm_pow_scaled[0], cm_pow_scaled[1]);
}
#endif

/* ---- integer64 kernels (int64.inc.c), against the scalar forms ---------- */

static int64_t lsmall[N + 1], lsmall2[N + 1]; /* within +-2^31, and edges */
static double fbig[N + 1];                     /* doubles around the int64 range */

static void init_i64_inputs(void) {
  const int64_t edges[] = {INT64_C(3037000500), INT64_C(3037000499), -INT64_C(3037000499),
                           INT64_C(4294967296), INT64_C(2147483648), -INT64_C(2147483648),
                           INT64_MAX, -INT64_MAX, INT64_MIN, 0, 1, -1};
  const double dedges[] = {9223372036854775808.0, -9223372036854775808.0, 9.5, -9.5, 0.0, -0.0,
                           HUGE_VAL, -HUGE_VAL, from_bits(UINT64_C(0x7FF8000000000000)),
                           18446744073709555712.0, -9223372036854779904.0,
                           9223372036854777856.0, 1e300, -1e300, 4503599627370497.5};
  const int ne = (int) (sizeof edges / sizeof edges[0]);
  const int nd = (int) (sizeof dedges / sizeof dedges[0]);
  int i;
  for (i = 0; i < N; i++) {
    if (i < 2 * ne) {
      lsmall[i] = edges[i % ne];
      lsmall2[i] = edges[(i * 5 + 3) % ne];
    } else {
      lsmall[i] = (int64_t) (next_rand() % UINT64_C(4294967295)) - INT64_C(2147483647);
      lsmall2[i] = i % 9 == 0 ? INT64_MIN
                              : (int64_t) (next_rand() % UINT64_C(4294967295)) - INT64_C(2147483647);
    }
    fbig[i] = i < 2 * nd ? dedges[i % nd]
                         : ((double) (int64_t) (next_rand() >> 11) / 4503599627370496.0 - 1.0) *
                             18446744073709551616.0 * 1.25;
  }
  /* Results exactly INT64_MIN: add, sub, mul. */
  lsmall[0] = -INT64_MAX;
  lsmall2[0] = -1;
  lsmall[1] = -INT64_MAX;
  lsmall2[1] = 1;
  lsmall[2] = -INT64_C(4294967296);
  lsmall2[2] = INT64_C(2147483648);
}

static void test_int64(ptrdiff_t n) {
  static const rsimd_opts oc = {0, 1, 0}, onc = {0, 0, 0}, orm = {1, 1, 0};
  const rsimd_opts *opts[3] = {&onc, &oc, &orm};
  char what[96];
  ptrdiff_t j;
  int op, f, k, v, mode, st, want;

  /* Reductions, with each of na_check off, on, and na.rm. */
  for (v = 0; v < 3; v++) {
    const rsimd_opts *o = opts[v];
    const int check = o->na_check || o->na_rm;
    int stop;
    rsimd_reduce_result r, ref;
    const int64_t *x = v == 2 ? lsmall2 : la;
    memset(&r, 0, sizeof r);
    memset(&ref, 0, sizeof ref);
    RSIMD_KERNEL(sum_i64)(x, n, &r, o);
    for (j = 0; j < n; j++) {
      if (check && x[j] == INT64_MIN) {
        ref.saw_na = 1;
        if (o->na_rm) continue;
      }
      rsimd_add_carry_i64(&ref.i64, &ref.carry, x[j]);
      ref.count++;
    }
    snprintf(what, sizeof what, "sum_i64 variant %d", v);
    check_int(what, n, 0, (long) (r.i64 >> 32), (long) (ref.i64 >> 32));
    check_int(what, n, 1, (long) (r.i64 & 0xFFFFFFFF), (long) (ref.i64 & 0xFFFFFFFF));
    check_int(what, n, 2, (long) r.carry, (long) ref.carry);
    check_flags(what, n, &r, &ref);

    reduce_init(&r, 0.0);
    reduce_init(&ref, 0.0);
    RSIMD_KERNEL(minmax_i64)(x, n, &r, o);
    for (j = 0; j < n; j++) {
      if (check && x[j] == INT64_MIN) {
        ref.saw_na = 1;
        continue;
      }
      if (x[j] < ref.i64) ref.i64 = x[j];
      if (x[j] > ref.i64_hi) ref.i64_hi = x[j];
      ref.count++;
    }
    snprintf(what, sizeof what, "minmax_i64 variant %d", v);
    check_int(what, n, 0, r.i64 == ref.i64, 1);
    check_int(what, n, 1, r.i64_hi == ref.i64_hi, 1);
    check_flags(what, n, &r, &ref);

    for (stop = RSIMD_STOP_NONE; stop <= RSIMD_STOP_FALSE; stop++) {
      memset(&r, 0, sizeof r);
      memset(&ref, 0, sizeof ref);
      RSIMD_KERNEL(anyall_i64)(x, n, stop, &r, o);
      for (j = 0; j < n; j++) {
        if (check && x[j] == INT64_MIN) ref.saw_na = 1;
        else if (x[j] == 0) {
          ref.any_false = 1;
          if (stop == RSIMD_STOP_FALSE) break;
        } else {
          ref.any_true = 1;
          if (stop == RSIMD_STOP_TRUE) break;
        }
      }
      /* Only the answer that the stop decides is defined after an early
         stop (any for TRUE, all for FALSE). */
      snprintf(what, sizeof what, "anyall_i64 variant %d stop %d", v, stop);
      if (stop != RSIMD_STOP_FALSE) {
        check_int(what, n, 0, lgl3(&r, 0, o->na_rm), lgl3(&ref, 0, o->na_rm));
      }
      if (stop != RSIMD_STOP_TRUE) {
        check_int(what, n, 1, lgl3(&r, 1, o->na_rm), lgl3(&ref, 1, o->na_rm));
      }
    }

    memset(&r, 0, sizeof r);
    memset(&ref, 0, sizeof ref);
    RSIMD_KERNEL(popcnt_sum_i64)(x, n, &r, o);
    for (j = 0; j < n; j++) {
      if (check && x[j] == INT64_MIN) {
        ref.saw_na = 1;
        if (!o->na_rm) break;
        continue;
      }
      ref.i64 += rsimd_popcnt64_1((uint64_t) x[j]);
    }
    snprintf(what, sizeof what, "popcnt_sum_i64 variant %d", v);
    check_int(what, n, 0, r.saw_na, ref.saw_na);
    if (!ref.saw_na || o->na_rm) check_int(what, n, 1, (long) r.i64, (long) ref.i64);
  }
  for (j = 0; j < n; j += 7) {
    R_xlen_t got = RSIMD_KERNEL(find_i64)(la, n, la[j]), ref = 0;
    while (la[ref] != la[j]) ref++;
    check_int("find_i64", n, j, (long) got, (long) ref);
  }
  check_int("find_i64 absent", n, 0, (long) RSIMD_KERNEL(find_i64)(la, n, INT64_C(77)), -1L);
  for (mode = RSIMD_NAMODE_ANY; mode <= RSIMD_NAMODE_WHICH_I32; mode++) {
    rsimd_reduce_result r;
    int32_t idx[N + 1];
    ptrdiff_t count = 0;
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(na_i64)(lsmall2, n, mode, 10, idx, &r);
    for (j = 0; j < n; j++) {
      if (lsmall2[j] != INT64_MIN) continue;
      if (mode == RSIMD_NAMODE_WHICH_I32) check_int("na_i64 which", n, j, idx[count], (long) j + 11);
      count++;
    }
    snprintf(what, sizeof what, "na_i64 mode %d", mode);
    if (mode == RSIMD_NAMODE_ANY) check_int(what, n, 0, r.saw_na, count > 0);
    else check_int(what, n, 0, (long) r.i64, (long) count);
  }
#ifndef RSIMD_SKIP_scan_i64
  for (op = 0; op <= 3; op += op == 0 ? 2 : 1) {
    rsimd_scan_state s, sref;
    R_xlen_t stopped;
    memset(&s, 0, sizeof s);
    s.i64 = op == 2 ? INT64_MAX : op == 3 ? -INT64_MAX : 0;
    sref = s;
    reset_out();
    stopped = RSIMD_KERNEL(scan_i64)(op, lsmall2, n, lout, &s);
    for (j = 0; j < n; j++) {
      if (rsimd_scan_i64_1(op, lsmall2[j], lref + j, &sref)) break;
    }
    snprintf(what, sizeof what, "scan_i64 op %d", op);
    check_int(what, n, 0, (long) stopped, j < n ? (long) j : -1L);
    for (; j <= n; j++) lref[j] = lout[j];
    check_i64(what, n, lout, lref);
  }
#endif

  /* Elementwise: every op with each broadcast and int32-operand flag
     combination, on full-range and small inputs, with and without NA
     checks. */
  for (op = RSIMD_EW_NEG; op <= RSIMD_EW_SIGN; op++) {
    reset_out();
    RSIMD_KERNEL(ew1_i64)(op, la, n, lout, &oc);
    for (j = 0; j < n; j++) {
      lref[j] = op == RSIMD_EW_NEG ? rsimd_neg_i64(la[j])
              : op == RSIMD_EW_ABS ? rsimd_abs_i64(la[j]) : rsimd_sign_i64(la[j]);
    }
    snprintf(what, sizeof what, "ew1_i64 op %d", op);
    check_i64(what, n, lout, lref);
  }
  for (v = 0; v < 2; v++) {
    const int64_t *X = v ? lsmall : la, *Y = v ? lsmall2 : lb;
    for (op = RSIMD_EW_ADD; op <= RSIMD_EW_MUL_WRAP; op++) {
      if (op == RSIMD_EW_DIV || op == RSIMD_EW_COPYSIGN) continue;
      for (f = 0; f < 32; f++) {
        const int sc = f & 3, ti = f >> 3;
        const void *x = ti & 1 ? (const void *) ia : (const void *) X;
        const void *y = ti & 2 ? (const void *) ib : (const void *) Y;
        int c;
        if (sc == 3 || (f & 4)) continue;
        for (c = 0; c < 2; c++) {
          const int check = op >= RSIMD_EW_PMIN && op <= RSIMD_EW_PMAX_NUM ? 1 : c;
          reset_out();
          st = RSIMD_KERNEL(ew2_i64)(op, x, y, n, f, lout, c ? &oc : &onc);
          want = 0;
          for (j = 0; j < n; j++) {
            lref[j] = rsimd_ew2_i64_1(op, rsimd_i64_get(x, f, 0, j, check),
                                      rsimd_i64_get(y, f, 1, j, check), check, &want);
          }
          snprintf(what, sizeof what, "ew2_i64 op %d flags %d check %d input %d", op, f, c, v);
          check_i64(what, n, lout, lref);
          check_int(what, n, 0, st, want);
        }
      }
    }
    for (op = RSIMD_EW_MUL_ADD; op <= RSIMD_EW_CLAMP; op++) {
      if (op == RSIMD_EW_LERP) continue;
      for (f = 0; f < 7; f++) {
        const int64_t *Z = op == RSIMD_EW_CLAMP ? Y : la;
        const int check = 1;
        reset_out();
        st = RSIMD_KERNEL(ew3_i64)(op, X, op == RSIMD_EW_CLAMP ? X : Y, Z, n, f, lout, &oc);
        want = 0;
        for (j = 0; j < n; j++) {
          lref[j] = rsimd_ew3_i64_1(op, rsimd_i64_get(X, f, 0, j, check),
                                    rsimd_i64_get(op == RSIMD_EW_CLAMP ? X : Y, f, 1, j, check),
                                    rsimd_i64_get(Z, f, 2, j, check), check, &want);
        }
        snprintf(what, sizeof what, "ew3_i64 op %d flags %d input %d", op, f, v);
        check_i64(what, n, lout, lref);
        check_int(what, n, 0, st, want);
      }
    }
    for (op = RSIMD_CMP_EQ; op <= RSIMD_CMP_GE; op++) {
      for (f = 0; f < 32; f++) {
        const int sc = f & 3, ti = f >> 3;
        const void *x = ti & 1 ? (const void *) ia : (const void *) X;
        const void *y = ti & 2 ? (const void *) ib : (const void *) Y;
        if (sc == 3 || (f & 4)) continue;
        reset_out();
        RSIMD_KERNEL(cmp_i64)(op, x, y, n, f, iout);
        for (j = 0; j < n; j++) {
          iref[j] = rsimd_cmp_i64_1(op, rsimd_i64_get(x, f, 0, j, 1), rsimd_i64_get(y, f, 1, j, 1));
        }
        snprintf(what, sizeof what, "cmp_i64 op %d flags %d input %d", op, f, v);
        check_i32(what, n, iout, iref);
      }
    }
    for (op = RSIMD_PRED_NA; op <= RSIMD_PRED_POW2; op++) {
      int pm;
      if (op == RSIMD_PRED_NAN || op == RSIMD_PRED_INFINITE) continue;
      for (pm = RSIMD_PRED_ELT; pm <= RSIMD_PRED_ALL; pm++) {
        int got, ref = pm == RSIMD_PRED_ALL;
        reset_out();
        got = RSIMD_KERNEL(pred_i64)(op, Y, n, pm, iout);
        for (j = 0; j < n; j++) {
          iref[j] = rsimd_pred_i64_1(op, Y[j]);
          if (pm == RSIMD_PRED_ANY && iref[j]) ref = 1;
          if (pm == RSIMD_PRED_ALL && !iref[j]) ref = 0;
        }
        snprintf(what, sizeof what, "pred_i64 op %d mode %d input %d", op, pm, v);
        if (pm == RSIMD_PRED_ELT) check_i32(what, n, iout, iref);
        else check_int(what, n, 0, got, ref);
      }
    }
  }
  for (op = RSIMD_BIT_AND; op <= RSIMD_BIT_TZCNT; op++) {
    const int binary = op <= RSIMD_BIT_XOR, counts = rsimd_bit_counts_i64(op);
    const int kmax = op >= RSIMD_BIT_SHL && op <= RSIMD_BIT_ROTR ? 63 : 0;
    for (f = 0; f < (binary ? 32 : 1); f++) {
      const int sc = f & 3, ti = f >> 3;
      const void *x = ti & 1 ? (const void *) ia : (const void *) la;
      const void *y = binary ? (ti & 2 ? (const void *) ib : (const void *) lb) : NULL;
      int c;
      if (sc == 3 || (f & 4)) continue;
      for (k = 0; k <= kmax; k++) {
        for (c = 0; c < 2; c++) {
          reset_out();
          RSIMD_KERNEL(bit_i64)(op, x, y, n, f, k, counts ? (void *) iout : (void *) lout,
                                c ? &oc : &onc);
          for (j = 0; j < n; j++) {
            int64_t a = rsimd_i64_get(x, f, 0, j, 1), b = binary ? rsimd_i64_get(y, f, 1, j, 1) : 0;
            int na = c && rsimd_na2_i64(a, b);
            if (counts) iref[j] = na ? RSIMD_NA_I32 : (int32_t) rsimd_bit_i64_1(op, a, b, k);
            else lref[j] = na ? INT64_MIN : rsimd_bit_i64_1(op, a, b, k);
          }
          snprintf(what, sizeof what, "bit_i64 op %d flags %d k %d check %d", op, f, k, c);
          if (counts) check_i32(what, n, iout, iref);
          else check_i64(what, n, lout, lref);
        }
      }
    }
  }

  /* Conversions, every mode, through the slot's integer64 part. */
  for (mode = RSIMD_CVT_CHECKED; mode <= RSIMD_CVT_TRUNCATING; mode++) {
    for (op = RSIMD_CVT_I64_F64; op <= RSIMD_CVT_U8_I64; op++) {
      for (v = 0; v < 3; v++) {
        const void *x = op == RSIMD_CVT_F64_I64 ? (const void *) (v == 0 ? fbig : v == 1 ? fa : fconv)
                      : op == RSIMD_CVT_I32_I64 ? (const void *) ia
                      : op == RSIMD_CVT_U8_I64 ? (const void *) la
                      : (const void *) (v == 0 ? la : v == 1 ? lsmall : lsmall2);
        void *out, *ref;
        size_t size;
        int i;
        if (op == RSIMD_CVT_I64_F64 || op == RSIMD_CVT_F64_I64) size = 8;
        else size = op == RSIMD_CVT_I64_U8 ? 1 : op == RSIMD_CVT_I32_I64 || op == RSIMD_CVT_U8_I64 ? 8 : 4;
        out = size == 8 ? (void *) lout : (void *) iout;
        ref = size == 8 ? (void *) lref : (void *) iref;
        reset_out();
        for (i = 0; i <= N; i++) {
          lref[i] = SENTINEL_I64;
          iref[i] = SENTINEL_I32;
        }
        if (size == 1) {
          memset(iout, 0x5A, sizeof iout);
          memset(iref, 0x5A, sizeof iref);
        }
        st = RSIMD_KERNEL(convert_i64_)(op, mode, x, n, out);
        want = 0;
        for (j = 0; j < n; j++) want |= rsimd_cvt_i64_1(op, mode, x, j, ref);
        snprintf(what, sizeof what, "convert_i64 op %d mode %d input %d", op, mode, v);
        check_int(what, n, 0, st, want);
        if (size == 8) check_i64(what, n, lout, lref);
        else if (size == 4) check_i32(what, n, iout, iref);
        else check_int(what, n, 1, memcmp(iout, iref, (size_t) n + 1) == 0, 1);
      }
    }
  }
}

/* ---- Hamming distances ----------------------------------------------------- */

/* Operands from few values, so that pairs are often equal: clean (no
   missing value) and dirty (about one in 13 missing) variants. */
static double hf[2][2][N + 1];
static int32_t hi[2][2][N + 1];
static int64_t hl[2][2][N + 1];
static Rcomplex hz[2][2][N + 1];

static void init_hamming_inputs(void) {
  const double fv[] = {1.0, 2.0, -0.0, 0.0, 3.5};
  const int32_t iv[] = {1, 2, 0, -1, 7};
  const int64_t lv[] = {1, 2, 0, INT64_MAX, -INT64_MAX};
  int d, k, i;
  for (d = 0; d < 2; d++) {
    for (k = 0; k < 2; k++) {
      for (i = 0; i <= N; i++) {
        unsigned r = (unsigned) (next_rand() % 65), a = r % 5;
        int miss = d && r % 13 == 0;
        double nan = r % 2 ? rsimd_na_real() : from_bits(UINT64_C(0x7FF8000000000000));
        hf[d][k][i] = miss ? nan : fv[a];
        hi[d][k][i] = miss ? RSIMD_NA_I32 : iv[a];
        hl[d][k][i] = miss ? RSIMD_NA_I64 : lv[a];
        hz[d][k][i].r = miss && r % 3 == 0 ? nan : fv[a];
        hz[d][k][i].i = miss && r % 3 != 0 ? nan : fv[(r / 5) % 2];
      }
    }
  }
}

/* r against the reference count (want, missing): the missing flag
   always, the count unless the kernel may have stopped at a missing
   pair. */
static void check_hamming(const char *what, ptrdiff_t n, const rsimd_reduce_result *r, long want,
                          int missing, int na_rm) {
  check_int(what, n, 0, r->saw_na, missing);
  if (!missing || na_rm) check_int(what, n, 1, (long) r->i64, want);
}

static void test_hamming(ptrdiff_t n) {
  static const rsimd_opts keep = {0, 1, 0}, rm = {1, 1, 0};
  rsimd_reduce_result r;
  char what[96];
  ptrdiff_t j;
  int d, f, narm;
  if (n == 0) return;
  for (d = 0; d < 2; d++) {
    for (narm = 0; narm < 2; narm++) {
      const rsimd_opts *o = narm ? &rm : &keep;
      /* int32 operands, each broadcast or not */
      for (f = 0; f < 3; f++) {
        long want = 0;
        int miss = 0;
        for (j = 0; j < n; j++) {
          int32_t a = hi[d][0][f & 1 ? 0 : j], b = hi[d][1][f & 2 ? 0 : j];
          if (a == RSIMD_NA_I32 || b == RSIMD_NA_I32) miss = 1;
          else want += a != b;
        }
        memset(&r, 0, sizeof r);
        RSIMD_KERNEL(hamming_i32)(hi[d][0], hi[d][1], n, f, &r, o);
        snprintf(what, sizeof what, "hamming_i32 flags %d input %d na_rm %d", f, d, narm);
        check_hamming(what, n, &r, want, miss, narm);
      }
      /* doubles, with int32 operands (flags bits 3-4) */
#ifndef RSIMD_SKIP_hamming_f64
      for (f = 0; f < 32; f++) {
        const int sc = f & 3, ti = f >> 3;
        const void *x = ti & 1 ? (const void *) hi[d][0] : (const void *) hf[d][0];
        const void *y = ti & 2 ? (const void *) hi[d][1] : (const void *) hf[d][1];
        long want = 0;
        int miss = 0;
        if (sc == 3 || (f & 4)) continue;
        for (j = 0; j < n; j++) {
          double a = ew_get(x, ti & 1, sc & 1, j), b = ew_get(y, ti & 2, sc & 2, j);
          if (isnan(a) || isnan(b)) miss = 1;
          else want += a != b;
        }
        memset(&r, 0, sizeof r);
        RSIMD_KERNEL(hamming_f64)(x, y, n, f, &r, o);
        snprintf(what, sizeof what, "hamming_f64 flags %d input %d na_rm %d", f, d, narm);
        check_hamming(what, n, &r, want, miss, narm);
      }
#endif
      /* integer64, with int32 operands */
      for (f = 0; f < 32; f++) {
        const int sc = f & 3, ti = f >> 3;
        const void *x = ti & 1 ? (const void *) hi[d][0] : (const void *) hl[d][0];
        const void *y = ti & 2 ? (const void *) hi[d][1] : (const void *) hl[d][1];
        long want = 0;
        int miss = 0;
        if (sc == 3 || (f & 4)) continue;
        for (j = 0; j < n; j++) {
          int64_t a = rsimd_i64_get(x, f, 0, j, 1), b = rsimd_i64_get(y, f, 1, j, 1);
          if (a == RSIMD_NA_I64 || b == RSIMD_NA_I64) miss = 1;
          else want += a != b;
        }
        memset(&r, 0, sizeof r);
        RSIMD_KERNEL(hamming_i64)(x, y, n, f, &r, o);
        snprintf(what, sizeof what, "hamming_i64 flags %d input %d na_rm %d", f, d, narm);
        check_hamming(what, n, &r, want, miss, narm);
      }
      /* complex */
#ifndef RSIMD_SKIP_hamming_c128
      for (f = 0; f < 3; f++) {
        long want = 0;
        int miss = 0;
        for (j = 0; j < n; j++) {
          Rcomplex a = hz[d][0][f & 1 ? 0 : j], b = hz[d][1][f & 2 ? 0 : j];
          if (isnan(a.r) || isnan(a.i) || isnan(b.r) || isnan(b.i)) miss = 1;
          else want += a.r != b.r || a.i != b.i;
        }
        memset(&r, 0, sizeof r);
        RSIMD_KERNEL(hamming_c128)(hz[d][0], hz[d][1], n, f, &r, o);
        snprintf(what, sizeof what, "hamming_c128 flags %d input %d na_rm %d", f, d, narm);
        check_hamming(what, n, &r, want, miss, narm);
      }
#endif
    }
  }
  /* bytes and bit distances, each operand broadcast or not */
  for (f = 0; f < 3; f++) {
    long want = 0, wb = 0, wi = 0, wl = 0;
    for (j = 0; j < n; j++) {
      uint8_t a = ua[f & 1 ? 0 : j], b = ub[f & 2 ? 0 : j];
      want += a != b;
      wb += rsimd_popcount32((uint32_t) (a ^ b));
      wi += rsimd_popcount32((uint32_t) ia[f & 1 ? 0 : j] ^ (uint32_t) ib[f & 2 ? 0 : j]);
      wl += rsimd_popcnt64_1((uint64_t) la[f & 1 ? 0 : j] ^ (uint64_t) lb[f & 2 ? 0 : j]);
    }
    snprintf(what, sizeof what, "hamming_u8 flags %d", f);
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(hamming_u8)(ua, ub, n, f, &r);
    check_hamming(what, n, &r, want, 0, 0);
    snprintf(what, sizeof what, "hamming_bits_u8 flags %d", f);
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(hamming_bits_u8)(ua, ub, n, f, &r);
    check_hamming(what, n, &r, wb, 0, 0);
    snprintf(what, sizeof what, "hamming_bits_i32 flags %d", f);
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(hamming_bits_i32)(ia, ib, n, f, &r);
    check_hamming(what, n, &r, wi, 0, 0);
    snprintf(what, sizeof what, "hamming_bits_i64 flags %d", f);
    memset(&r, 0, sizeof r);
    RSIMD_KERNEL(hamming_bits_i64)(la, lb, n, f, &r);
    check_hamming(what, n, &r, wl, 0, 0);
  }
}

int main(void) {
  ptrdiff_t n;
  init_inputs();
  init_logical_inputs();
  init_i64_inputs();
  init_hamming_inputs();
#ifndef RSIMD_NO_F64_SIMD
  init_mod_inputs();
#endif
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
    test_arith(n);
    test_logical(n);
    test_int64(n);
#ifndef RSIMD_NO_F64_SIMD
    test_complex(n);
    test_carith(n);
#if !defined(RSIMD_SKIP_math1_c128) && defined(RSIMD_HAVE_SLEEF)
    test_cmath(n);
#endif
    if (n <= 2 * RSIMD_LANES_64 + 1 || n == N) test_cmath_fns(n);
#endif
    test_hamming(n);
#ifndef RSIMD_NO_F64_SIMD
    test_ilogb(n);
#endif
  }
  test_bytes_exhaustive();
#ifndef RSIMD_NO_F64_SIMD
  test_na_cancel();
#endif
#ifndef RSIMD_NO_F64_SIMD
  test_f64_horizontal();
  test_uzp();
  test_zip();
#endif
  test_int_horizontal();
  test_intdiv_const();
#ifndef RSIMD_NO_F64_SIMD
  test_math_extras();
  test_cmath_sweep();
  cm_report();
#endif
#ifdef RSIMD_HAVE_SLEEF
  test_sleef();
  test_pi_oracle();
  test_math();
  test_math_exact();
  test_softmax();
#endif
  printf("tier %s: lanes64=%ld lanes32=%ld, %ld checks, %ld failures", RSIMD_TIER_STRING,
         (long) RSIMD_LANES_64, (long) RSIMD_LANES_32, n_checks, n_fail);
#ifdef RSIMD_HAVE_SLEEF
  printf("; SLEEF worst %.3f ULP (%s), %.3f ULP (%s), %.3f ULP (%s)", sleef_worst[0],
         sleef_worst_name[0], sleef_worst[1], sleef_worst_name[1], sleef_worst[2],
         sleef_worst_name[2]);
  printf("; none sinpi/cospi/tanpi worst %.3f/%.3f/%.3f ULP", pi_oracle_worst[0],
         pi_oracle_worst[1], pi_oracle_worst[2]);
#endif
#ifndef RSIMD_NO_F64_SIMD
  printf("; approx worst recip 2^%.2f rsqrt 2^%.2f; rootn worst %.3f ULP", log2(approx_worst[0]),
         log2(approx_worst[1]), rootn_worst);
#endif
  printf("\n");
  return n_fail != 0;
}
