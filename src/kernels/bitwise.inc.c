/* Logical and bitwise kernels (op codes in kernel_types.h).
 *
 * Three-valued logic (logic_f64, logic_i32): each operand element is read
 * as TRUE (non-zero), FALSE (zero) or NA (NA_integer_, or any NaN for
 * doubles), and
 *   x & y   is FALSE if either is FALSE, else NA if either is NA, else TRUE;
 *   x | y   is TRUE if either is TRUE, else NA if either is NA, else FALSE;
 *   xor     is NA if either is NA, else whether exactly one is TRUE;
 *   !x      is NA for NA, else whether x is FALSE.
 *
 * Integer bitwise ops (bit_i32) work on the uint32 bit pattern of the
 * elements, except SAR, which propagates the sign bit. With na_check an NA
 * operand gives NA; any result whose bit pattern is 0x80000000 is NA as
 * well, because that is NA_integer_'s pattern (bitwShiftL(1L, 31L) is NA in
 * base R too). Counts: POPCNT is the number of set bits, LZCNT and TZCNT
 * the number of leading and trailing zero bits (32 for 0).
 *
 * Byte ops (bit_u8) act on 8-bit values and never see NA: shifts lose the
 * bits shifted out (a shift by 8 gives 0), rotates are 8-bit rotates, and
 * the counts are in 0..8. The vector tiers compute bytes in 32-bit lanes.
 */

#include "logical.inc.h"

/* ---- Scalar forms: the none tier, and references ------------------------- */

static inline int32_t rsimd_popcnt_1(uint32_t x) { return rsimd_popcount32(x); }
static inline int32_t rsimd_lzcnt_1(uint32_t x) {
  int32_t k = 0;
  if (x == 0) return 32;
  while (!(x & 0x80000000u)) {
    x <<= 1;
    k++;
  }
  return k;
}
static inline int32_t rsimd_tzcnt_1(uint32_t x) {
  int32_t k = 0;
  if (x == 0) return 32;
  while (!(x & 1u)) {
    x >>= 1;
    k++;
  }
  return k;
}

/* Logical value of an element: 1, 0 or NA. */
static inline int rsimd_lgl_of_i32(int32_t x) { return x == RSIMD_NA_I32 ? x : x != 0; }
static inline int rsimd_lgl_of_f64(double x) { return isnan(x) ? RSIMD_NA_I32 : x != 0; }

/* Three-valued op of the logical values a and b. */
static inline int rsimd_logic_1(int op, int a, int b) {
  const int na = RSIMD_NA_I32;
  switch (op) {
  case RSIMD_LOGIC_AND:
    if (a == 0 || b == 0) return 0;
    return a == na || b == na ? na : 1;
  case RSIMD_LOGIC_OR:
    if (a == 1 || b == 1) return 1;
    return a == na || b == na ? na : 0;
  case RSIMD_LOGIC_XOR: return a == na || b == na ? na : a != b;
  default: return a == na ? na : !a;
  }
}

/* Bit op of int32 elements (uint32 patterns), without NA handling. */
static inline int32_t rsimd_bit_i32_1(int op, int32_t a, int32_t b, int k) {
  uint32_t u = (uint32_t) a, v = (uint32_t) b, r;
  switch (op) {
  case RSIMD_BIT_AND: r = u & v; break;
  case RSIMD_BIT_OR: r = u | v; break;
  case RSIMD_BIT_XOR: r = u ^ v; break;
  case RSIMD_BIT_NOT: r = ~u; break;
  case RSIMD_BIT_SHL: r = u << k; break;
  case RSIMD_BIT_SHR: r = u >> k; break;
  case RSIMD_BIT_SAR: r = a < 0 ? ~(~u >> k) : u >> k; break;
  case RSIMD_BIT_ROTL: r = k == 0 ? u : (u << k) | (u >> (32 - k)); break;
  case RSIMD_BIT_ROTR: r = k == 0 ? u : (u >> k) | (u << (32 - k)); break;
  case RSIMD_BIT_POPCNT: return rsimd_popcnt_1(u);
  case RSIMD_BIT_LZCNT: return rsimd_lzcnt_1(u);
  default: return rsimd_tzcnt_1(u);
  }
  /* Two's complement reinterpretation without implementation-defined
     conversion of values above INT32_MAX. */
  return r <= 0x7FFFFFFFu ? (int32_t) r : (int32_t) (r - 0x80000000u) - 0x7FFFFFFF - 1;
}

/* Bit op of bytes: the result byte, or the count for the counting ops. */
static inline int32_t rsimd_bit_u8_1(int op, uint8_t a, uint8_t b, int k) {
  unsigned u = a, v = b;
  switch (op) {
  case RSIMD_BIT_AND: return (int32_t) (u & v);
  case RSIMD_BIT_OR: return (int32_t) (u | v);
  case RSIMD_BIT_XOR: return (int32_t) (u ^ v);
  case RSIMD_BIT_NOT: return (int32_t) (~u & 0xFFu);
  case RSIMD_BIT_SHL: return (int32_t) ((u << k) & 0xFFu);
  case RSIMD_BIT_SHR: return (int32_t) (u >> k);
  case RSIMD_BIT_ROTL: return (int32_t) (((u << k) | (u >> (8 - k))) & 0xFFu);
  case RSIMD_BIT_ROTR: return (int32_t) (((u >> k) | (u << (8 - k))) & 0xFFu);
  case RSIMD_BIT_POPCNT: return rsimd_popcnt_1(u);
  case RSIMD_BIT_LZCNT: return rsimd_lzcnt_1(u) - 24;
  default: return u == 0 ? 8 : rsimd_tzcnt_1(u);
  }
}

static inline int rsimd_bit_counts(int op) {
  return op == RSIMD_BIT_POPCNT || op == RSIMD_BIT_LZCNT || op == RSIMD_BIT_TZCNT;
}

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(logic_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             int *out);
void RSIMD_KERNEL(logic_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int a = rsimd_lgl_of_f64(rsimd_ew_get(x, flags, 0, i, 1));
    int b = op == RSIMD_LOGIC_NOT ? 0 : rsimd_lgl_of_f64(rsimd_ew_get(y, flags, 1, i, 1));
    out[i] = rsimd_logic_1(op, a, b);
  }
}

void RSIMD_KERNEL(logic_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                             int *out);
void RSIMD_KERNEL(logic_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                             int *out) {
  const R_xlen_t sx = RSIMD_LGL_SCALAR(0) ? 0 : 1, sy = RSIMD_LGL_SCALAR(1) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int a = rsimd_lgl_of_i32(x[i * sx]);
    int b = op == RSIMD_LOGIC_NOT ? 0 : rsimd_lgl_of_i32(y[i * sy]);
    out[i] = rsimd_logic_1(op, a, b);
  }
}

void RSIMD_KERNEL(bit_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int k,
                           int *out, const rsimd_opts *o);
void RSIMD_KERNEL(bit_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int k,
                           int *out, const rsimd_opts *o) {
  const int binary = op <= RSIMD_BIT_XOR, check = o->na_check;
  const R_xlen_t sx = RSIMD_LGL_SCALAR(0) ? 0 : 1, sy = RSIMD_LGL_SCALAR(1) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int32_t a = x[i * sx], b = binary ? y[i * sy] : 0;
    if (check && (a == RSIMD_NA_I32 || b == RSIMD_NA_I32)) out[i] = RSIMD_NA_I32;
    else out[i] = rsimd_bit_i32_1(op, a, b, k);
  }
}

void RSIMD_KERNEL(bit_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags, int k,
                          void *out);
void RSIMD_KERNEL(bit_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags, int k,
                          void *out) {
  const int binary = op <= RSIMD_BIT_XOR;
  const R_xlen_t sx = RSIMD_LGL_SCALAR(0) ? 0 : 1, sy = RSIMD_LGL_SCALAR(1) ? 0 : 1;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int32_t r = rsimd_bit_u8_1(op, x[i * sx], binary ? y[i * sy] : 0, k);
    if (rsimd_bit_counts(op)) ((int *) out)[i] = r;
    else ((Rbyte *) out)[i] = (Rbyte) r;
  }
}

void RSIMD_KERNEL(popcnt_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o);
void RSIMD_KERNEL(popcnt_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I32) {
      r->saw_na = 1;
      if (!o->na_rm) return;
      continue;
    }
    r->i64 += rsimd_popcnt_1((uint32_t) x[i]);
  }
}

void RSIMD_KERNEL(popcnt_sum_u8)(const Rbyte *x, R_xlen_t n, rsimd_reduce_result *r);
void RSIMD_KERNEL(popcnt_sum_u8)(const Rbyte *x, R_xlen_t n, rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) r->i64 += rsimd_popcnt_1(x[i]);
}

#else /* vector tiers */

/* ---- Bit counts of int32 lanes -------------------------------------------- */

#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
RSIMD_INLINE rsimd_vi32 rsimd_vi32_popcnt(rsimd_vi32 x) {
  return svreinterpret_s32_u32(svcnt_s32_x(RSIMD_PT32, x));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_clz(rsimd_vi32 x) {
  return svreinterpret_s32_u32(svclz_s32_x(RSIMD_PT32, x));
}
#elif RSIMD_TIER_IS(neon) && defined(SIMDE_ARM_NEON_A32V7_NATIVE)
/* Byte counts summed pairwise into 32-bit lanes; NEON's own clz. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_popcnt(rsimd_vi32 x) {
  uint8x16_t c = vcntq_u8(vreinterpretq_u8_u32(simde__m128i_to_neon_u32(x)));
  return simde__m128i_from_neon_u32(vpaddlq_u16(vpaddlq_u8(c)));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_clz(rsimd_vi32 x) {
  return simde__m128i_from_neon_u32(vclzq_u32(simde__m128i_to_neon_u32(x)));
}
#else
/* SWAR population count (Hacker's Delight, figure 5-2), on x86 where no
   tier has a 32-bit lane popcount (AVX512-VPOPCNTDQ and AVX512-CD's lzcnt
   are not part of the avx512 tier). */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_popcnt(rsimd_vi32 x) {
  const rsimd_vi32 m1 = rsimd_vi32_set1(0x55555555), m2 = rsimd_vi32_set1(0x33333333),
                   m4 = rsimd_vi32_set1(0x0F0F0F0F);
  x = rsimd_vi32_sub(x, rsimd_vi32_and(rsimd_vi32_srl(x, 1), m1));
  x = rsimd_vi32_add(rsimd_vi32_and(x, m2), rsimd_vi32_and(rsimd_vi32_srl(x, 2), m2));
  x = rsimd_vi32_and(rsimd_vi32_add(x, rsimd_vi32_srl(x, 4)), m4);
  x = rsimd_vi32_add(x, rsimd_vi32_srl(x, 8));
  x = rsimd_vi32_add(x, rsimd_vi32_srl(x, 16));
  return rsimd_vi32_and(x, rsimd_vi32_set1(0x3F));
}
/* Leading zeros: smear the highest set bit down, then 32 - popcount. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_clz(rsimd_vi32 x) {
  x = rsimd_vi32_or(x, rsimd_vi32_srl(x, 1));
  x = rsimd_vi32_or(x, rsimd_vi32_srl(x, 2));
  x = rsimd_vi32_or(x, rsimd_vi32_srl(x, 4));
  x = rsimd_vi32_or(x, rsimd_vi32_srl(x, 8));
  x = rsimd_vi32_or(x, rsimd_vi32_srl(x, 16));
  return rsimd_vi32_sub(rsimd_vi32_set1(32), rsimd_vi32_popcnt(x));
}
#endif
/* Trailing zeros: the bits below the lowest set bit, ~x & (x - 1). */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_ctz(rsimd_vi32 x) {
  return rsimd_vi32_popcnt(rsimd_vi32_andnot(x, rsimd_vi32_sub(x, rsimd_vi32_set1(1))));
}

/* ---- Three-valued logic --------------------------------------------------- */

/* The result of op from the TRUE masks ta, tb and the NA masks na, nb
   (FALSE where neither): sets the masks t (TRUE) and m (NA). */
#define RSIMD_LOGIC_MASKS(M, t, m, op, ta, tb, na, nb)                                     \
  do {                                                                                     \
    if ((op) == RSIMD_LOGIC_AND) {                                                         \
      /* FALSE unless both are TRUE or NA */                                               \
      M f_ = M##_not(M##_and(M##_or(ta, na), M##_or(tb, nb)));                             \
      t = M##_and(ta, tb);                                                                 \
      m = M##_andnot(f_, M##_or(na, nb));                                                  \
    } else if ((op) == RSIMD_LOGIC_OR) {                                                   \
      t = M##_or(ta, tb);                                                                  \
      m = M##_andnot(t, M##_or(na, nb));                                                   \
    } else if ((op) == RSIMD_LOGIC_XOR) {                                                  \
      m = M##_or(na, nb);                                                                  \
      t = M##_or(M##_andnot(tb, ta), M##_andnot(ta, tb));                                  \
    } else {                                                                               \
      m = na;                                                                              \
      t = M##_not(M##_or(ta, na));                                                         \
    }                                                                                      \
  } while (0)

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(logic_i32_)(const int op, const int *x, const int *y,
                                                  R_xlen_t n, int flags, int *out) {
  RSIMD_LANE32_LOOP(int32_t, rsimd_vi32_loadu, rsimd_vi32_loadu_p, int32_t, rsimd_vi32_storeu,
                    rsimd_vi32_storeu_p, {
                      const rsimd_vi32 z_ = rsimd_vi32_zero();
                      rsimd_mi32 na = rsimd_vi32_is_na(a), nb = rsimd_vi32_is_na(b), t, m;
                      rsimd_mi32 ta = rsimd_mi32_not(rsimd_mi32_or(na, rsimd_vi32_cmp_eq(a, z_)));
                      rsimd_mi32 tb = rsimd_mi32_not(rsimd_mi32_or(nb, rsimd_vi32_cmp_eq(b, z_)));
                      RSIMD_LOGIC_MASKS(rsimd_mi32, t, m, op, ta, tb, na, nb);
                      r = rsimd_lgl_vi32(t, m);
                    });
}

void RSIMD_KERNEL(logic_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                             int *out);
void RSIMD_KERNEL(logic_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags,
                             int *out) {
  switch (op) {
  case RSIMD_LOGIC_AND: RSIMD_KERNEL(logic_i32_)(RSIMD_LOGIC_AND, x, y, n, flags, out); break;
  case RSIMD_LOGIC_OR: RSIMD_KERNEL(logic_i32_)(RSIMD_LOGIC_OR, x, y, n, flags, out); break;
  case RSIMD_LOGIC_XOR: RSIMD_KERNEL(logic_i32_)(RSIMD_LOGIC_XOR, x, y, n, flags, out); break;
  default: RSIMD_KERNEL(logic_i32_)(RSIMD_LOGIC_NOT, x, NULL, n, flags, out); break;
  }
}

#ifndef RSIMD_NO_F64_SIMD
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(logic_f64_)(const int op, const void *x, const void *y,
                                                  R_xlen_t n, int flags, int *out) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_ew_bcast(y, flags, 1, 1);
  RSIMD_LANE64_LOOP_FLAGS({
    const rsimd_vf64 z_ = rsimd_vf64_zero();
    rsimd_mf64 na = rsimd_vf64_is_nan(a), nb = rsimd_vf64_is_nan(b), t, m;
    rsimd_mf64 ta = rsimd_mf64_not(rsimd_mf64_or(na, rsimd_vf64_cmp_eq(a, z_)));
    rsimd_mf64 tb = rsimd_mf64_not(rsimd_mf64_or(nb, rsimd_vf64_cmp_eq(b, z_)));
    RSIMD_LOGIC_MASKS(rsimd_mf64, t, m, op, ta, tb, na, nb);
    r = rsimd_lgl_vi64(t, m);
  });
}

void RSIMD_KERNEL(logic_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             int *out);
void RSIMD_KERNEL(logic_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             int *out) {
  switch (op) {
  case RSIMD_LOGIC_AND: RSIMD_KERNEL(logic_f64_)(RSIMD_LOGIC_AND, x, y, n, flags, out); break;
  case RSIMD_LOGIC_OR: RSIMD_KERNEL(logic_f64_)(RSIMD_LOGIC_OR, x, y, n, flags, out); break;
  case RSIMD_LOGIC_XOR: RSIMD_KERNEL(logic_f64_)(RSIMD_LOGIC_XOR, x, y, n, flags, out); break;
  default: RSIMD_KERNEL(logic_f64_)(RSIMD_LOGIC_NOT, x, NULL, n, flags, out); break;
  }
}
#else
#define RSIMD_SKIP_logic_f64 1
#endif

#undef RSIMD_LOGIC_MASKS

/* ---- Bitwise ops ---------------------------------------------------------- */

/* op of the int32 lanes a and b (count k), without NA handling. */
RSIMD_ALWAYS_INLINE rsimd_vi32 rsimd_bit_vi32(const int op, rsimd_vi32 a, rsimd_vi32 b, int k) {
  switch (op) {
  case RSIMD_BIT_AND: return rsimd_vi32_and(a, b);
  case RSIMD_BIT_OR: return rsimd_vi32_or(a, b);
  case RSIMD_BIT_XOR: return rsimd_vi32_xor(a, b);
  case RSIMD_BIT_NOT: return rsimd_vi32_xor(a, rsimd_vi32_set1(-1));
  case RSIMD_BIT_SHL: return rsimd_vi32_sll(a, k);
  case RSIMD_BIT_SHR: return rsimd_vi32_srl(a, k);
  case RSIMD_BIT_SAR: return rsimd_vi32_sra(a, k);
  case RSIMD_BIT_ROTL:
    return k == 0 ? a : rsimd_vi32_or(rsimd_vi32_sll(a, k), rsimd_vi32_srl(a, 32 - k));
  case RSIMD_BIT_ROTR:
    return k == 0 ? a : rsimd_vi32_or(rsimd_vi32_srl(a, k), rsimd_vi32_sll(a, 32 - k));
  case RSIMD_BIT_POPCNT: return rsimd_vi32_popcnt(a);
  case RSIMD_BIT_LZCNT: return rsimd_vi32_clz(a);
  default: return rsimd_vi32_ctz(a);
  }
}

/* op of bytes held in 32-bit lanes: the result's low 8 bits are the byte
   (the store drops the rest), or the lane is the count. */
RSIMD_ALWAYS_INLINE rsimd_vi32 rsimd_bit_vu8(const int op, rsimd_vi32 a, rsimd_vi32 b, int k) {
  switch (op) {
  case RSIMD_BIT_AND: return rsimd_vi32_and(a, b);
  case RSIMD_BIT_OR: return rsimd_vi32_or(a, b);
  case RSIMD_BIT_XOR: return rsimd_vi32_xor(a, b);
  case RSIMD_BIT_NOT: return rsimd_vi32_xor(a, rsimd_vi32_set1(0xFF));
  case RSIMD_BIT_SHL: return rsimd_vi32_sll(a, k);
  case RSIMD_BIT_SHR: return rsimd_vi32_srl(a, k);
  case RSIMD_BIT_ROTL: return rsimd_vi32_or(rsimd_vi32_sll(a, k), rsimd_vi32_srl(a, 8 - k));
  case RSIMD_BIT_ROTR: return rsimd_vi32_or(rsimd_vi32_srl(a, k), rsimd_vi32_sll(a, 8 - k));
  case RSIMD_BIT_POPCNT: return rsimd_vi32_popcnt(a);
  case RSIMD_BIT_LZCNT: return rsimd_vi32_sub(rsimd_vi32_clz(a), rsimd_vi32_set1(24));
  default: return rsimd_vi32_min(rsimd_vi32_ctz(a), rsimd_vi32_set1(8));
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(bit_i32_)(const int op, const int *x, const int *y,
                                                R_xlen_t n, int flags, int k, int *out,
                                                const int check) {
  RSIMD_LANE32_LOOP(int32_t, rsimd_vi32_loadu, rsimd_vi32_loadu_p, int32_t, rsimd_vi32_storeu,
                    rsimd_vi32_storeu_p, {
                      r = rsimd_bit_vi32(op, a, b, k);
                      if (check) {
                        r = rsimd_vi32_set_na(r, op <= RSIMD_BIT_XOR ? rsimd_vi32_na2(a, b)
                                                                     : rsimd_vi32_is_na(a));
                      }
                    });
}

/* Calls KERNEL(op, ...) with op as a constant. */
#define RSIMD_BIT_DISPATCH(KERNEL, ...)                                                    \
  switch (op) {                                                                            \
  case RSIMD_BIT_AND: KERNEL(RSIMD_BIT_AND, __VA_ARGS__); break;                           \
  case RSIMD_BIT_OR: KERNEL(RSIMD_BIT_OR, __VA_ARGS__); break;                             \
  case RSIMD_BIT_XOR: KERNEL(RSIMD_BIT_XOR, __VA_ARGS__); break;                           \
  case RSIMD_BIT_NOT: KERNEL(RSIMD_BIT_NOT, __VA_ARGS__); break;                           \
  case RSIMD_BIT_SHL: KERNEL(RSIMD_BIT_SHL, __VA_ARGS__); break;                           \
  case RSIMD_BIT_SHR: KERNEL(RSIMD_BIT_SHR, __VA_ARGS__); break;                           \
  case RSIMD_BIT_SAR: KERNEL(RSIMD_BIT_SAR, __VA_ARGS__); break;                           \
  case RSIMD_BIT_ROTL: KERNEL(RSIMD_BIT_ROTL, __VA_ARGS__); break;                         \
  case RSIMD_BIT_ROTR: KERNEL(RSIMD_BIT_ROTR, __VA_ARGS__); break;                         \
  case RSIMD_BIT_POPCNT: KERNEL(RSIMD_BIT_POPCNT, __VA_ARGS__); break;                     \
  case RSIMD_BIT_LZCNT: KERNEL(RSIMD_BIT_LZCNT, __VA_ARGS__); break;                       \
  default: KERNEL(RSIMD_BIT_TZCNT, __VA_ARGS__); break;                                    \
  }

void RSIMD_KERNEL(bit_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int k,
                           int *out, const rsimd_opts *o);
void RSIMD_KERNEL(bit_i32)(int op, const int *x, const int *y, R_xlen_t n, int flags, int k,
                           int *out, const rsimd_opts *o) {
  if (op > RSIMD_BIT_XOR) y = NULL;
  if (o->na_check) {
    RSIMD_BIT_DISPATCH(RSIMD_KERNEL(bit_i32_), x, y, n, flags, k, out, 1)
  } else {
    RSIMD_BIT_DISPATCH(RSIMD_KERNEL(bit_i32_), x, y, n, flags, k, out, 0)
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(bit_u8_)(const int op, const Rbyte *x, const Rbyte *y,
                                               R_xlen_t n, int flags, int k, void *out) {
  if (rsimd_bit_counts(op)) {
    RSIMD_LANE32_LOOP(uint8_t, rsimd_vi32_loadu_u8, rsimd_vi32_loadu_u8_p, int32_t,
                      rsimd_vi32_storeu, rsimd_vi32_storeu_p, r = rsimd_bit_vu8(op, a, b, k));
  } else {
    RSIMD_LANE32_LOOP(uint8_t, rsimd_vi32_loadu_u8, rsimd_vi32_loadu_u8_p, uint8_t,
                      rsimd_vi32_storeu_u8, rsimd_vi32_storeu_u8_p, r = rsimd_bit_vu8(op, a, b, k));
  }
}

void RSIMD_KERNEL(bit_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags, int k,
                          void *out);
void RSIMD_KERNEL(bit_u8)(int op, const Rbyte *x, const Rbyte *y, R_xlen_t n, int flags, int k,
                          void *out) {
  if (op > RSIMD_BIT_XOR) y = NULL;
  RSIMD_BIT_DISPATCH(RSIMD_KERNEL(bit_u8_), x, y, n, flags, k, out)
}

#undef RSIMD_BIT_DISPATCH

/* ---- Population count totals ---------------------------------------------- */

/* Lane counts are added in 32-bit lanes and moved to r->i64 at least every
   2^24 vectors, long before a lane could overflow. */
#define RSIMD_POPCNT_FLUSH ((ptrdiff_t) 1 << 24)

void RSIMD_KERNEL(popcnt_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o);
void RSIMD_KERNEL(popcnt_sum_i32)(const int *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  rsimd_vi32 acc = rsimd_vi32_zero();
  ptrdiff_t i = 0, since = 0;
  for (; i < n; i += RSIMD_LANES_32) {
    /* The fill 0 counts nothing and is not NA. */
    rsimd_vi32 v = i + RSIMD_LANES_32 <= n ? rsimd_vi32_loadu(x + i)
                                           : rsimd_vi32_loadu_p(rsimd_p32_while(i, n), x + i, 0);
    rsimd_vi32 c = rsimd_vi32_popcnt(v);
    if (check) {
      rsimd_mi32 m = rsimd_vi32_is_na(v);
      if (rsimd_mi32_any(m)) {
        r->saw_na = 1;
        if (!o->na_rm) break;
        c = rsimd_vi32_blend(c, rsimd_vi32_zero(), m);
      }
    }
    acc = rsimd_vi32_add(acc, c);
    if (++since == RSIMD_POPCNT_FLUSH) {
      r->i64 += rsimd_vi32_reduce_add(acc);
      acc = rsimd_vi32_zero();
      since = 0;
    }
  }
  r->i64 += rsimd_vi32_reduce_add(acc);
}

void RSIMD_KERNEL(popcnt_sum_u8)(const Rbyte *x, R_xlen_t n, rsimd_reduce_result *r);
void RSIMD_KERNEL(popcnt_sum_u8)(const Rbyte *x, R_xlen_t n, rsimd_reduce_result *r) {
  rsimd_vi32 acc = rsimd_vi32_zero();
  ptrdiff_t i = 0, since = 0;
  for (; i < n; i += RSIMD_LANES_32) {
    rsimd_vi32 v = i + RSIMD_LANES_32 <= n ? rsimd_vi32_loadu_u8(x + i)
                                           : rsimd_vi32_loadu_u8_p(rsimd_p32_while(i, n), x + i, 0);
    acc = rsimd_vi32_add(acc, rsimd_vi32_popcnt(v));
    if (++since == RSIMD_POPCNT_FLUSH) {
      r->i64 += rsimd_vi32_reduce_add(acc);
      acc = rsimd_vi32_zero();
      since = 0;
    }
  }
  r->i64 += rsimd_vi32_reduce_add(acc);
}

#undef RSIMD_POPCNT_FLUSH

#endif /* vector tiers */
