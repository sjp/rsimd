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
 * those files define for the none tier. On tiers built with SLEEF the
 * elementary-function wrappers are compared with long double libm, and the
 * elementary-function kernels (src/kernels/math.inc.c) with the scalar
 * forms that file defines for the none tier. Prints a
 * summary and exits non-zero on any mismatch.
 */

#include <float.h>
#include <stdio.h>
#include "kernels/common.inc.h"
#include "kernels/reduce.inc.c"
#include "kernels/scan.inc.c"
#include "kernels/arith.inc.c"
#include "kernels/predicates.inc.c"
#include "kernels/compare.inc.c"
#include "kernels/bitwise.inc.c"
#include "kernels/convert.inc.c"
#include "kernels/math.inc.c"

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

#ifdef RSIMD_HAVE_SLEEF
/* The SLEEF wrappers of common.inc.h, through the predicated loop, against
   SLEEF's documented error bounds: 1.0 ULP for the _u10 functions and 0.5
   ULP for the _u05 ones, measured against the long double libm function
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
static double sleef_worst[2] = {0, 0};
static const char *sleef_worst_name[2] = {"", ""};

static void check_ulp(const char *what, ptrdiff_t n, const double *got, double bound) {
  ptrdiff_t i;
  char buf[200];
  const int k = bound < 1;
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
   (libm, as the none tier): within 2 ULP (4 for cbrt), with NaN kinds (NA or NaN),
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
  default: r = rsimd_math2_na_f64(hypot(a, b), a, b); break;
  }
  if (isnan(r) && !isnan(a) && !isnan(b)) math_ref_status = RSIMD_EW_NAN_PRODUCED;
  return r;
}
static int math_exact_pow(double a, double b) {
  return b == 0 || b == 2 || b == 3 || b == 4 || a == 0 || a == 1 || !isfinite(a) ||
         !isfinite(b);
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

/* log2 and exp2 of the powers of 2 and integers, log10 and exp10 of
   10^(0:22) and 0:22: exact, as in base R. */
static void test_math_exact(void) {
  static double xs[2100], got[2100], want[2100];
  static const struct { int op, lo, hi, pow; } cases[] = {
    {RSIMD_MATH_LOG2, -1074, 1023, 2}, {RSIMD_MATH_EXP2, -1074, 1023, 0},
    {RSIMD_MATH_LOG10, 0, 22, 10}, {RSIMD_MATH_EXP10, 0, 22, 0}};
  size_t c;
  int k, m;
  char buf[96];
  for (c = 0; c < sizeof cases / sizeof cases[0]; c++) {
    for (m = 0, k = cases[c].lo; k <= cases[c].hi; k++, m++) {
      double p = cases[c].op == RSIMD_MATH_LOG10 || cases[c].op == RSIMD_MATH_EXP10
                   ? pow(10.0, k) : ldexp(1.0, k);
      xs[m] = cases[c].pow ? p : (double) k;
      want[m] = cases[c].pow ? (double) k : p;
    }
    RSIMD_KERNEL(math1_f64)(cases[c].op, xs, m, 0, 1.0, got);
    for (k = 0; k < m; k++) {
      n_checks++;
      if (bits(got[k]) != bits(want[k])) {
        snprintf(buf, sizeof buf, "input %a: got %a want %a", xs[k], got[k], want[k]);
        fail("math exact powers", m, k, buf);
      }
    }
  }
}

static void test_math(void) {
  static const ptrdiff_t lens[] = {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 65, N};
  int op, f, st;
  size_t li;
  ptrdiff_t j, n;
  char what[96];
  for (li = 0; li < sizeof lens / sizeof lens[0]; li++) {
    n = lens[li];
    for (op = RSIMD_MATH_EXP; op <= RSIMD_MATH_ATANH; op++) {
      const double p = op == RSIMD_MATH_LOGB ? log(3.0) : 1.0;
      math_inputs(op, ma, mia, n);
      for (f = 0; f < 2; f++) {
        reset_out();
        st = RSIMD_KERNEL(math1_f64)(op, f ? (const void *) mia : ma, n,
                                     f ? RSIMD_EW_I32(0) : 0, p, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          mref[j] = math_ref1(op, f ? (mia[j] == RSIMD_NA_I32 ? rsimd_na_real() : mia[j]) : ma[j], p);
        }
        snprintf(what, sizeof what, "math1 op %d int32 %d", op, f);
        /* glibc's cbrt is up to 3 ULP from the exact value on aarch64. */
        check_math(what, n, fout, mref, ma, NULL, 0, op == RSIMD_MATH_CBRT ? 4 : 2);
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
        st = RSIMD_KERNEL(math2_f64)(op, x, y, n, flags, fout);
        math_ref_status = 0;
        for (j = 0; j < n; j++) {
          ptrdiff_t jx = (f & 1) ? 0 : j, jy = (f & 2) ? 0 : j;
          xa[j] = (f & 4) ? (mia[jx] == RSIMD_NA_I32 ? rsimd_na_real() : mia[jx]) : ma[jx];
          yb[j] = (f & 8) ? (mib[jy] == RSIMD_NA_I32 ? rsimd_na_real() : mib[jy]) : mb[jy];
          mref[j] = math_ref2(op, xa[j], yb[j]);
        }
        snprintf(what, sizeof what, "math2 op %d flags %d", op, f);
        check_math(what, n, fout, mref, xa, yb, op == RSIMD_MATH_POW, 2);
        check_int(what, n, 0, st, math_ref_status);
      }
    }
    /* sincos: both outputs. */
    sleef_inputs(ma, n, -1e3, 1e3, 0, 1);
    if (n > 9) ma[9] = rsimd_na_real();
    reset_out();
    for (j = 0; j <= n; j++) mout2[j] = SENTINEL_F64;
    st = RSIMD_KERNEL(sincos_f64)(ma, n, 0, fout, mout2);
    math_ref_status = 0;
    for (j = 0; j < n; j++) {
      mref[j] = math_ref1(RSIMD_MATH_SIN, ma[j], 1);
      mref2[j] = math_ref1(RSIMD_MATH_COS, ma[j], 1);
    }
    check_math("sincos sin", n, fout, mref, ma, NULL, 0, 2);
    check_math("sincos cos", n, mout2, mref2, ma, NULL, 0, 2);
    check_int("sincos status", n, 0, st, math_ref_status);
  }
}
#endif /* RSIMD_HAVE_SLEEF */

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
  for (op = RSIMD_EW_FMA; op <= RSIMD_EW_CLAMP; op++) {
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
#ifndef RSIMD_NO_F64_SIMD
static double pf[N + 1];
#endif

static void init_logical_inputs(void) {
  const double cvals[] = {
    2147483647.9, 2147483648.0, -2147483648.0, -2147483648.5, -2147483647.5, 2147483647.0,
    -2147483647.0, 255.9, 256.0, 255.0, -0.5, -1.0, -0.0, 0.0, 0.999, 4294967301.0,
    -4294967301.0, 1e19, -1e19, 9007199254740993.0, 1e300, HUGE_VAL, -HUGE_VAL,
    from_bits(UINT64_C(0x7FF8000000000000)), from_bits(UINT64_C(0x7FF00000000007A2)),
    4.9e-324, -4.9e-324, 65535.5, -255.5, 2.5, 3e9, 2147483904.0, -2147483904.0, 6442450944.0};
  const int nc = (int) (sizeof cvals / sizeof cvals[0]);
  int i;
  for (i = 0; i < N; i++) {
    if (i < 2 * nc) {
      fcv[i] = cvals[(i * 7) % nc];
    } else {
      double m = (double) (int64_t) (next_rand() >> 11) / 9007199254740992.0 - 0.5;
      fcv[i] = ldexp(m, (int) (next_rand() % 72));
    }
    icv[i] = i % 11 == 5 ? RSIMD_NA_I32 : (int32_t) (next_rand() % 400) - 70;
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
  for (op = RSIMD_PRED_NA; op <= RSIMD_PRED_ZERO; op++) {
    if (op == RSIMD_PRED_NAN || op == RSIMD_PRED_INFINITE) continue;
    for (f = 0; f < 2; f++) {
      const int32_t *x = f ? pa : ia;
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
  for (op = RSIMD_PRED_NA; op <= RSIMD_PRED_ZERO; op++) {
    for (f = 0; f < 3; f++) {
      const double *x = f == 0 ? fa : f == 1 ? fcv : pf;
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

int main(void) {
  ptrdiff_t n;
  init_inputs();
  init_logical_inputs();
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
  }
  test_bytes_exhaustive();
#ifndef RSIMD_NO_F64_SIMD
  test_na_cancel();
#endif
#ifndef RSIMD_NO_F64_SIMD
  test_f64_horizontal();
#endif
  test_int_horizontal();
#ifdef RSIMD_HAVE_SLEEF
  test_sleef();
  test_pi_oracle();
  test_math();
  test_math_exact();
#endif
  printf("tier %s: lanes64=%ld lanes32=%ld, %ld checks, %ld failures", RSIMD_TIER_STRING,
         (long) RSIMD_LANES_64, (long) RSIMD_LANES_32, n_checks, n_fail);
#ifdef RSIMD_HAVE_SLEEF
  printf("; SLEEF worst %.3f ULP (%s), %.3f ULP (%s)", sleef_worst[0], sleef_worst_name[0],
         sleef_worst[1], sleef_worst_name[1]);
  printf("; none sinpi/cospi/tanpi worst %.3f/%.3f/%.3f ULP", pi_oracle_worst[0],
         pi_oracle_worst[1], pi_oracle_worst[2]);
#endif
  printf("\n");
  return n_fail != 0;
}
