/* bit64::integer64 kernels: int64_t elements whose NA is INT64_MIN (slots
 * in kernel_list.h, op codes and status bits in kernel_types.h).
 *
 * The scalar forms below (rsimd_*_i64) are shared by every tier: the none
 * tier's kernels are plain loops over them, and the vector tiers use them
 * for the ops and lanes they leave to scalar code, so both give identical
 * results. Operands of the elementwise, comparison and bit kernels may be
 * int32 elements (flag RSIMD_EW_I32(k), integer or logical), which are
 * sign-extended, NA_integer_ becoming NA. Arithmetic honours na_check
 * there (as the int32 kernels do); pmin, pmax, comparisons and the bit
 * ops always treat INT64_MIN as NA. In the predicated last vector the
 * inactive lanes repeat the first active element of each operand, as in
 * arith.inc.c, so they cannot raise a status bit that no element raised.
 *
 * How each op runs on each vector tier. "native": instructions of the
 * tier; "SIMDe": SIMDe's emulation of an instruction the tier lacks;
 * "partial": built from 32-bit partial products or shifts; "scalar": a
 * loop over the scalar forms. Choices on sse2, avx2, avx512 and sve were
 * made from the instruction sequences and are not measured (only neon is
 * timed on the development machine).
 *
 *   op                       sse2     avx2     avx512   neon     sve/sve2
 *   sum (exact, carries)     native   native   native   native   native
 *     (the NA test uses 64-bit compares: SIMDe on sse2, only when NA
 *     checking is on)
 *   min/max (also of |x|),   SIMDe    native   native   native   native
 *     pmin/pmax, clamp,      (64-bit compares are SSE4.2; min/max are
 *     comparisons            compare + blend except on avx512 and sve)
 *   any/all, is_na ...,      native   native   native   native   native
 *     any_na/count/which
 *   add, sub (checked/wrap)  native   native   native   native   native
 *     (overflow tests use 64-bit compares: SIMDe on sse2)
 *   mul_wrap                 partial  partial  native   partial  native
 *   mul (checked)            partial  native   native   native   native
 *     (lanes with operands beyond +-2^31 use the scalar form, except on
 *     sve, which has the high half of the product)
 *   neg, abs, sign           native   native   native   native   native
 *   idiv, mod                partial  partial  native   native   native
 *     (in double lanes, for blocks of operands within +-2^51; other
 *     blocks are scalar, by a multiply for a scalar divisor except on
 *     AArch64; the int64 <-> double conversions are partial on sse2 and
 *     avx2)
 *   mul_add, add_mul,        scalar   scalar   scalar   scalar   scalar
 *     cumsum, cummin, cummax
 *   and, or, xor, not,       native   native   native   native   native
 *     shl, shr, rotl, rotr
 *   sar                      partial  partial  partial  partial  native
 *   popcount, hamming_bits   partial  partial  partial  native   native
 *   is_normal, is_even,      SIMDe    native   native   native   native
 *     is_odd, is_pow2,
 *     hamming (64-bit
 *     compares)
 *   lzcnt, tzcnt             partial  partial  partial  partial  native
 *   integer64 -> double      partial  partial  native   native   native
 *   double -> integer64      scalar   scalar   scalar   native*  native*
 *     (* checked and saturating; truncating is scalar)
 *   integer64 -> integer,    native   native   native   native   native
 *     logical; integer ->
 *     integer64
 *   integer64 <-> raw        scalar   scalar   scalar   scalar   scalar
 *
 * "partial" sar, popcount and lzcnt: the arithmetic shift from logical
 * shifts and the sign, the population count by the SWAR method (bytes
 * summed with psadbw on sse2 and avx2), the leading zeros as 64 minus the
 * population count of the highest set bit smeared down. 32-bit ARM builds
 * the same kernels from the 128-bit layer (NEON's 64-bit lanes); there the
 * double conversions are scalar.
 */

#include "logical.inc.h"

/* ---- Scalar forms (every tier) ------------------------------------------- */

static inline int64_t rsimd_wrap_add_i64(int64_t x, int64_t y) {
  return (int64_t) ((uint64_t) x + (uint64_t) y);
}
static inline int64_t rsimd_wrap_sub_i64(int64_t x, int64_t y) {
  return (int64_t) ((uint64_t) x - (uint64_t) y);
}
static inline int64_t rsimd_wrap_mul_i64(int64_t x, int64_t y) {
  return (int64_t) ((uint64_t) x * (uint64_t) y);
}
static inline int rsimd_na2_i64(int64_t x, int64_t y) {
  return x == RSIMD_NA_I64 || y == RSIMD_NA_I64;
}

/* Element i of operand k (int64 or, with RSIMD_EW_I32(k), int32), for the
   elementwise kernels; check maps NA_integer_ to NA. */
static inline int64_t rsimd_i64_get(const void *p, int flags, int k, R_xlen_t i, int check) {
  R_xlen_t j = (flags & RSIMD_EW_SCALAR(k)) ? 0 : i;
  if (flags & RSIMD_EW_I32(k)) {
    int32_t v = ((const int32_t *) p)[j];
    return check && v == RSIMD_NA_I32 ? RSIMD_NA_I64 : (int64_t) v;
  }
  return ((const int64_t *) p)[j];
}

/* Checked arithmetic: NA for an NA operand (with check), NA and *ovf = 1
   for a result outside the int64 range or equal to INT64_MIN. */
static inline int64_t rsimd_add_i64(int64_t x, int64_t y, int check, int *ovf) {
  int64_t r = rsimd_wrap_add_i64(x, y);
  if (check && rsimd_na2_i64(x, y)) return RSIMD_NA_I64;
  if (((x ^ r) & (y ^ r)) < 0 || r == RSIMD_NA_I64) {
    *ovf = 1;
    return RSIMD_NA_I64;
  }
  return r;
}
static inline int64_t rsimd_sub_i64(int64_t x, int64_t y, int check, int *ovf) {
  int64_t r = rsimd_wrap_sub_i64(x, y);
  if (check && rsimd_na2_i64(x, y)) return RSIMD_NA_I64;
  if (((x ^ y) & (x ^ r)) < 0 || r == RSIMD_NA_I64) {
    *ovf = 1;
    return RSIMD_NA_I64;
  }
  return r;
}
/* 1 if x * y does not fit in int64 (r is the wrapped product). */
static inline int rsimd_mul_ovf_i64(int64_t x, int64_t y, int64_t *r) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_mul_overflow(x, y, r);
#else
  *r = rsimd_wrap_mul_i64(x, y);
  if (x == 0 || y == 0) return 0;
  if ((x == -1 && y == INT64_MIN) || (y == -1 && x == INT64_MIN)) return 1;
  return *r / y != x;
#endif
}
static inline int64_t rsimd_mul_i64(int64_t x, int64_t y, int check, int *ovf) {
  int64_t r;
  int o = rsimd_mul_ovf_i64(x, y, &r);
  if (check && rsimd_na2_i64(x, y)) return RSIMD_NA_I64;
  if (o || r == RSIMD_NA_I64) {
    *ovf = 1;
    return RSIMD_NA_I64;
  }
  return r;
}
/* Wrapping: NA for an NA operand (with check). Negation and absolute value
   map NA (INT64_MIN) to itself and need no check. */
static inline int64_t rsimd_add_wrap_i64(int64_t x, int64_t y, int check) {
  return check && rsimd_na2_i64(x, y) ? RSIMD_NA_I64 : rsimd_wrap_add_i64(x, y);
}
static inline int64_t rsimd_sub_wrap_i64(int64_t x, int64_t y, int check) {
  return check && rsimd_na2_i64(x, y) ? RSIMD_NA_I64 : rsimd_wrap_sub_i64(x, y);
}
static inline int64_t rsimd_mul_wrap_i64(int64_t x, int64_t y, int check) {
  return check && rsimd_na2_i64(x, y) ? RSIMD_NA_I64 : rsimd_wrap_mul_i64(x, y);
}
static inline int64_t rsimd_neg_i64(int64_t x) { return (int64_t) (0u - (uint64_t) x); }
static inline int64_t rsimd_abs_i64(int64_t x) { return x < 0 ? rsimd_neg_i64(x) : x; }
static inline int64_t rsimd_sign_i64(int64_t x) {
  return x == RSIMD_NA_I64 ? x : (int64_t) ((x > 0) - (x < 0));
}

/* x %/% y (mod = 0) or x %% y (mod = 1) as bit64: floor division, the
   remainder taking the sign of y. NA for an NA operand (with check), and
   NA with *dz = 1 for y == 0. INT64_MIN %/% -1 (possible only without
   check) is NA. */
static inline int64_t rsimd_intdiv_i64(int64_t x, int64_t y, int mod, int check, int *dz) {
  int64_t q, r;
  if (check && rsimd_na2_i64(x, y)) return RSIMD_NA_I64;
  if (y == 0) {
    *dz = 1;
    return RSIMD_NA_I64;
  }
  if (y == -1) return mod ? 0 : rsimd_neg_i64(x);
  q = x / y;
  r = x % y;
  if (r != 0 && ((r < 0) != (y < 0))) {
    q -= 1;
    r += y;
  }
  return mod ? r : q;
}

/* The high 64 bits of the 128-bit product a * b. */
static inline int64_t rsimd_mulhi_i64(int64_t a, int64_t b) {
#if defined(__SIZEOF_INT128__)
  __extension__ typedef __int128 rsimd_i128;
  return (int64_t) (((rsimd_i128) a * b) >> 64);
#else
  /* The unsigned high half from 32-bit partial products, less b for a < 0
     and a for b < 0 (mod 2^64). */
  const uint64_t ua = (uint64_t) a, ub = (uint64_t) b, lo = UINT64_C(0xFFFFFFFF);
  uint64_t t = (ua >> 32) * (ub & lo) + (((ua & lo) * (ub & lo)) >> 32);
  uint64_t w = (t & lo) + (ua & lo) * (ub >> 32);
  uint64_t hi = (ua >> 32) * (ub >> 32) + (t >> 32) + (w >> 32);
  if (a < 0) hi -= ub;
  if (b < 0) hi -= ua;
  return (int64_t) hi;
#endif
}

/* Division of int64 values by a constant d other than 0, 1, -1 and
   INT64_MIN without dividing: rsimd_divmagic_i32 (na.h) for 64 bits. The
   truncated quotient of x is ((mulhi(m, x) + add * x) >> s) plus one when
   that is negative. */
typedef struct {
  int64_t d, m;
  int s, add;
} rsimd_divmagic_i64;
static inline int rsimd_divmagic_i64_ok(int64_t d) {
  return d != 0 && d != 1 && d != -1 && d != RSIMD_NA_I64;
}
static inline rsimd_divmagic_i64 rsimd_divmagic_i64_make(int64_t d) {
  const uint64_t two63 = UINT64_C(0x8000000000000000);
  rsimd_divmagic_i64 g;
  uint64_t ad = d < 0 ? 0u - (uint64_t) d : (uint64_t) d, t, anc, q1, r1, q2, r2, delta, m;
  int p = 63;
  t = two63 + ((uint64_t) d >> 63);
  anc = t - 1 - t % ad;
  q1 = two63 / anc;
  r1 = two63 - q1 * anc;
  q2 = two63 / ad;
  r2 = two63 - q2 * ad;
  do {
    p++;
    q1 *= 2;
    r1 *= 2;
    if (r1 >= anc) {
      q1++;
      r1 -= anc;
    }
    q2 *= 2;
    r2 *= 2;
    if (r2 >= ad) {
      q2++;
      r2 -= ad;
    }
    delta = ad - r2;
  } while (q1 < delta || (q1 == delta && r1 == 0));
  m = q2 + 1;
  if (d < 0) m = 0u - m;
  g.d = d;
  memcpy(&g.m, &m, sizeof m);
  g.s = p - 64;
  g.add = d > 0 && g.m < 0 ? 1 : (d < 0 && g.m > 0 ? -1 : 0);
  return g;
}
/* x %/% g->d (mod = 0) or x %% g->d (mod = 1), as rsimd_intdiv_i64 for an
   x that is not NA. */
static inline int64_t rsimd_divmagic_i64_div(int64_t x, const rsimd_divmagic_i64 *g, int mod) {
  uint64_t u = (uint64_t) rsimd_mulhi_i64(g->m, x);
  int64_t q, r;
  if (g->add > 0) u += (uint64_t) x;
  else if (g->add < 0) u -= (uint64_t) x;
  q = (int64_t) u;
  q = q < 0 ? ~(~q >> g->s) : q >> g->s;
  q += (int64_t) ((uint64_t) q >> 63);
  r = x - q * g->d;
  if (r != 0 && (r ^ g->d) < 0) {
    q -= 1;
    r += g->d;
  }
  return mod ? r : q;
}

/* %/% (mod = 0) or %% (mod = 1) of elements [i, e) of the operands of an
   elementwise kernel into out, one division per element, or with g (a
   scalar divisor of g->d) none; returns the status bits. Inlined with
   constant arguments, so that the tests of flags, mod and check leave the
   loop. */
RSIMD_ALWAYS_INLINE int rsimd_intdiv_i64_loop(const void *x, const void *y, R_xlen_t i, R_xlen_t e,
                                              const int flags, const int mod, const int check,
                                              const rsimd_divmagic_i64 *g, int64_t *out) {
  int dz = 0;
  if (g) {
    /* A copy, which the stores to out cannot alias. */
    const rsimd_divmagic_i64 gc = *g;
    for (; i < e; i++) {
      int64_t a = rsimd_i64_get(x, flags, 0, i, check);
      out[i] = check && a == RSIMD_NA_I64 ? RSIMD_NA_I64 : rsimd_divmagic_i64_div(a, &gc, mod);
    }
    return 0;
  }
  for (; i < e; i++) {
    out[i] = rsimd_intdiv_i64(rsimd_i64_get(x, flags, 0, i, check),
                              rsimd_i64_get(y, flags, 1, i, check), mod, check, &dz);
  }
  return dz ? RSIMD_EW_DIV_ZERO : 0;
}
/* rsimd_intdiv_i64_loop for any flags, with each form of the operands
   (int64 or int32 vector, or int64 scalar: an int32 scalar is read as
   int64 first) inlined with constant flags; a test of the flags in every
   element costs about as much as the division. */
static int rsimd_intdiv_i64_loop_any(const void *x, const void *y, R_xlen_t i, R_xlen_t e, int flags,
                                     int mod, int check, const rsimd_divmagic_i64 *g, int64_t *out) {
#define RSIMD_INTDIV_I64_CASE(f)                                                           \
  case f: return rsimd_intdiv_i64_loop(x, y, i, e, f, mod, check, g, out)
  switch (flags) {
    RSIMD_INTDIV_I64_CASE(0);
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_SCALAR(1));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_I32(1));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_I32(0));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_I32(0) | RSIMD_EW_I32(1));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_I32(0) | RSIMD_EW_SCALAR(1));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_SCALAR(0));
    RSIMD_INTDIV_I64_CASE(RSIMD_EW_SCALAR(0) | RSIMD_EW_I32(1));
  default: return rsimd_intdiv_i64_loop(x, y, i, e, flags, mod, check, g, out);
  }
#undef RSIMD_INTDIV_I64_CASE
}
/* rsimd_intdiv_i64_loop, through rsimd_intdiv_i64_loop_any unless flags
   is one of the forms inlined with constants. */
RSIMD_ALWAYS_INLINE int rsimd_intdiv_i64_run(const void *x, const void *y, R_xlen_t i, R_xlen_t e,
                                             const int flags, const int mod, const int check,
                                             const rsimd_divmagic_i64 *g, int64_t *out) {
  if (flags == 0 || flags == RSIMD_EW_SCALAR(1) || flags == RSIMD_EW_I32(1)) {
    return rsimd_intdiv_i64_loop(x, y, i, e, flags, mod, check, g, out);
  }
  return rsimd_intdiv_i64_loop_any(x, y, i, e, flags, mod, check, g, out);
}

/* The divider for a scalar divisor y that has one (else NULL), in *g.
   AArch64 divides faster than the divider runs (1e6 elements on the
   Apple M-series: 660 us against 780 us), so its vector tiers divide;
   the none tier uses the divider on every machine, and so tests it. */
#if (defined(__aarch64__) || defined(_M_ARM64)) && !RSIMD_TIER_IS(none)
#define RSIMD_I64_DIVMAGIC 0
#else
#define RSIMD_I64_DIVMAGIC 1
#endif
static inline const rsimd_divmagic_i64 *rsimd_intdiv_i64_magic(const void *y, int flags, int check,
                                                               rsimd_divmagic_i64 *g) {
  int64_t d;
  if (!RSIMD_I64_DIVMAGIC || !(flags & RSIMD_EW_SCALAR(1))) return NULL;
  d = rsimd_i64_get(y, flags, 1, 0, check);
  if (!rsimd_divmagic_i64_ok(d)) return NULL;
  *g = rsimd_divmagic_i64_make(d);
  return g;
}
/* A scalar operand k read into *v (as int64, NA_integer_ becoming NA with
   check); an int32 one is replaced by *v, so that it takes the loops for
   int64 ones (the others stay, as the loops are faster reading them). */
static inline const void *rsimd_intdiv_i64_scalar_arg(const void *p, int *flags, int k, int check,
                                                      int64_t *v) {
  if (!(*flags & RSIMD_EW_SCALAR(k))) return p;
  *v = rsimd_i64_get(p, *flags, k, 0, check);
  if (!(*flags & RSIMD_EW_I32(k))) return p;
  *flags &= ~RSIMD_EW_I32(k);
  return v;
}
/* The scalar kernel of %/% and %%: rsimd_intdiv_i64_loop over [0, n) with
   constant arguments for the common cases. */
static int rsimd_ew_intdiv_i64_scalar(const void *x, const void *y, R_xlen_t n, int flags, int mod,
                                      int check, int64_t *out) {
  int64_t x0 = 0, y0 = 0;
  rsimd_divmagic_i64 gs;
  const rsimd_divmagic_i64 *g;
  x = rsimd_intdiv_i64_scalar_arg(x, &flags, 0, check, &x0);
  y = rsimd_intdiv_i64_scalar_arg(y, &flags, 1, check, &y0);
  g = rsimd_intdiv_i64_magic(y, flags, check, &gs);
  if (g) {
    if (flags == RSIMD_EW_SCALAR(1)) {
      return rsimd_intdiv_i64_loop(x, y, 0, n, RSIMD_EW_SCALAR(1), mod, check, g, out);
    }
    return rsimd_intdiv_i64_loop_any(x, y, 0, n, flags, mod, check, g, out);
  }
  if (flags == RSIMD_EW_I32(1)) {
    return rsimd_intdiv_i64_loop(x, y, 0, n, RSIMD_EW_I32(1), mod, check, NULL, out);
  }
  if (flags != 0) return rsimd_intdiv_i64_loop_any(x, y, 0, n, flags, mod, check, NULL, out);
  if (mod) {
    return check ? rsimd_intdiv_i64_loop(x, y, 0, n, 0, 1, 1, NULL, out)
                 : rsimd_intdiv_i64_loop(x, y, 0, n, 0, 1, 0, NULL, out);
  }
  return check ? rsimd_intdiv_i64_loop(x, y, 0, n, 0, 0, 1, NULL, out)
               : rsimd_intdiv_i64_loop(x, y, 0, n, 0, 0, 0, NULL, out);
}

/* pmin and pmax: NA if either operand is; the _num forms return the other
   operand when one is NA. */
static inline int64_t rsimd_pmin_i64(int64_t x, int64_t y) {
  return rsimd_na2_i64(x, y) ? RSIMD_NA_I64 : (y < x ? y : x);
}
static inline int64_t rsimd_pmax_i64(int64_t x, int64_t y) {
  return rsimd_na2_i64(x, y) ? RSIMD_NA_I64 : (y > x ? y : x);
}
static inline int64_t rsimd_pmin_num_i64(int64_t x, int64_t y) {
  if (x == RSIMD_NA_I64) return y;
  if (y == RSIMD_NA_I64) return x;
  return y < x ? y : x;
}
static inline int64_t rsimd_pmax_num_i64(int64_t x, int64_t y) {
  if (x == RSIMD_NA_I64) return y;
  if (y == RSIMD_NA_I64) return x;
  return y > x ? y : x;
}

/* Binary op `op` (RSIMD_EW_*) of one pair; *st collects the status bits. */
static inline int64_t rsimd_ew2_i64_1(int op, int64_t a, int64_t b, int check, int *st) {
  int ovf = 0, dz = 0;
  int64_t r;
  switch (op) {
  case RSIMD_EW_ADD: r = rsimd_add_i64(a, b, check, &ovf); break;
  case RSIMD_EW_SUB: r = rsimd_sub_i64(a, b, check, &ovf); break;
  case RSIMD_EW_MUL: r = rsimd_mul_i64(a, b, check, &ovf); break;
  case RSIMD_EW_ADD_WRAP: r = rsimd_add_wrap_i64(a, b, check); break;
  case RSIMD_EW_SUB_WRAP: r = rsimd_sub_wrap_i64(a, b, check); break;
  case RSIMD_EW_MUL_WRAP: r = rsimd_mul_wrap_i64(a, b, check); break;
  case RSIMD_EW_IDIV: r = rsimd_intdiv_i64(a, b, 0, check, &dz); break;
  case RSIMD_EW_MOD: r = rsimd_intdiv_i64(a, b, 1, check, &dz); break;
  case RSIMD_EW_PMIN: r = rsimd_pmin_i64(a, b); break;
  case RSIMD_EW_PMAX: r = rsimd_pmax_i64(a, b); break;
  case RSIMD_EW_PMIN_NUM: r = rsimd_pmin_num_i64(a, b); break;
  case RSIMD_EW_PMAX_NUM: r = rsimd_pmax_num_i64(a, b); break;
  default: r = RSIMD_NA_I64; break;
  }
  if (ovf) *st |= RSIMD_EW_OVERFLOW;
  if (dz) *st |= RSIMD_EW_DIV_ZERO;
  return r;
}
/* Ternary op of one triple: an NA intermediate (an overflow, or an NA
   operand) stays NA without counting as a second overflow. */
static inline int64_t rsimd_ew3_i64_1(int op, int64_t a, int64_t b, int64_t c, int check,
                                      int *st) {
  int ovf = 0;
  int64_t r;
  switch (op) {
  case RSIMD_EW_MUL_ADD:
    r = rsimd_mul_i64(a, b, check, &ovf);
    if (r != RSIMD_NA_I64) r = rsimd_add_i64(r, c, check, &ovf);
    break;
  case RSIMD_EW_ADD_MUL:
    r = rsimd_add_i64(a, b, check, &ovf);
    if (r != RSIMD_NA_I64) r = rsimd_mul_i64(r, c, check, &ovf);
    break;
  case RSIMD_EW_CLAMP:
    if (!rsimd_na2_i64(b, c) && b > c) *st |= RSIMD_EW_LO_GT_HI;
    r = rsimd_pmin_i64(rsimd_pmax_i64(a, b), c);
    break;
  default: r = RSIMD_NA_I64; break;
  }
  if (ovf) *st |= RSIMD_EW_OVERFLOW;
  return r;
}

/* Comparison of one pair as a logical value. */
static inline int rsimd_cmp_i64_1(int op, int64_t a, int64_t b) {
  if (rsimd_na2_i64(a, b)) return RSIMD_NA_I32;
  switch (op) {
  case RSIMD_CMP_EQ: return a == b;
  case RSIMD_CMP_NE: return a != b;
  case RSIMD_CMP_LT: return a < b;
  case RSIMD_CMP_LE: return a <= b;
  case RSIMD_CMP_GT: return a > b;
  default: return a >= b;
  }
}

static inline int rsimd_pred_i64_1(int op, int64_t x) {
  switch (op) {
  case RSIMD_PRED_NA: return x == RSIMD_NA_I64;
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return x != RSIMD_NA_I64;
  case RSIMD_PRED_NEGATIVE: return x < 0 && x != RSIMD_NA_I64;
  case RSIMD_PRED_NORMAL: return x != 0 && x != RSIMD_NA_I64;
  case RSIMD_PRED_EVEN: return x % 2 == 0 && x != RSIMD_NA_I64;
  case RSIMD_PRED_ODD: return x % 2 != 0;
  case RSIMD_PRED_POW2: return x > 0 && (x & (x - 1)) == 0;
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return 0;
  default: return x == 0;
  }
}

/* Bit counts of a uint64 pattern (64 for no bits). */
static inline int32_t rsimd_popcnt64_1(uint64_t x) {
  return rsimd_popcount32((uint32_t) x) + rsimd_popcount32((uint32_t) (x >> 32));
}
/* With the compiler's builtins where there are some, else a loop. */
static inline int32_t rsimd_lzcnt64_1(uint64_t x) {
  int32_t k = 0;
  if (x == 0) return 64;
#if defined(__GNUC__) || defined(__clang__)
  (void) k;
  return (int32_t) __builtin_clzll(x);
#else
  while (!(x & UINT64_C(0x8000000000000000))) {
    x <<= 1;
    k++;
  }
  return k;
#endif
}
static inline int32_t rsimd_tzcnt64_1(uint64_t x) {
  int32_t k = 0;
  if (x == 0) return 64;
#if defined(__GNUC__) || defined(__clang__)
  (void) k;
  return (int32_t) __builtin_ctzll(x);
#else
  while (!(x & 1u)) {
    x >>= 1;
    k++;
  }
  return k;
#endif
}

/* The scalar Hamming loops over elements [i, n): the none tier's kernels
   and the vector tiers' tails. rsimd_ham_i64_from returns 1 when it stops
   at a missing pair (without na_rm). */
static inline int rsimd_ham_i64_from(const void *x, const void *y, R_xlen_t i, R_xlen_t n,
                                     int flags, rsimd_reduce_result *r, int na_rm) {
  for (; i < n; i++) {
    int64_t a = rsimd_i64_get(x, flags, 0, i, 1), b = rsimd_i64_get(y, flags, 1, i, 1);
    if (rsimd_na2_i64(a, b)) {
      r->saw_na = 1;
      if (!na_rm) return 1;
      continue;
    }
    r->i64 += a != b;
  }
  return 0;
}
static inline void rsimd_ham_bits_i64_from(const int64_t *x, const int64_t *y, R_xlen_t i,
                                           R_xlen_t n, int flags, rsimd_reduce_result *r) {
  const R_xlen_t sx = (flags & RSIMD_EW_SCALAR(0)) ? 0 : 1, sy = (flags & RSIMD_EW_SCALAR(1)) ? 0 : 1;
  for (; i < n; i++) r->i64 += rsimd_popcnt64_1((uint64_t) x[i * sx] ^ (uint64_t) y[i * sy]);
}
/* Bit op of int64 elements (uint64 patterns, SAR sign-propagating), without
   NA handling; a result pattern of 0x8000000000000000 is NA by itself. The
   counts are returned as int64 values in 0..64. */
static inline int64_t rsimd_bit_i64_1(int op, int64_t a, int64_t b, int k) {
  uint64_t u = (uint64_t) a, v = (uint64_t) b, r;
  switch (op) {
  case RSIMD_BIT_AND: r = u & v; break;
  case RSIMD_BIT_OR: r = u | v; break;
  case RSIMD_BIT_XOR: r = u ^ v; break;
  case RSIMD_BIT_NOT: r = ~u; break;
  case RSIMD_BIT_SHL: r = u << k; break;
  case RSIMD_BIT_SHR: r = u >> k; break;
  case RSIMD_BIT_SAR: r = a < 0 ? ~(~u >> k) : u >> k; break;
  case RSIMD_BIT_ROTL: r = k == 0 ? u : (u << k) | (u >> (64 - k)); break;
  case RSIMD_BIT_ROTR: r = k == 0 ? u : (u >> k) | (u << (64 - k)); break;
  case RSIMD_BIT_POPCNT: return rsimd_popcnt64_1(u);
  case RSIMD_BIT_LZCNT: return rsimd_lzcnt64_1(u);
  default: return rsimd_tzcnt64_1(u);
  }
  return (int64_t) r;
}
static inline int rsimd_bit_counts_i64(int op) { return op >= RSIMD_BIT_POPCNT; }

/* Adds b to the wrapped sum *s, counting a wrap past INT64_MAX (+1) or
   INT64_MIN (-1) in *carry. */
static inline void rsimd_add_carry_i64(int64_t *s, int64_t *carry, int64_t b) {
  int64_t a = *s, r = rsimd_wrap_add_i64(a, b);
  if (((a ^ r) & (b ^ r)) < 0) *carry += b < 0 ? -1 : 1;
  *s = r;
}

/* Adds h * 2^32 + l to the total (*s, *carry) of rsimd_add_carry_i64(),
   for |h| <= 2^62 and any l. */
static inline void rsimd_add_hl_i64(int64_t *s, int64_t *carry, int64_t h, uint64_t l) {
  uint64_t hl = (uint64_t) h & 0xFFFFFFFFu, u1 = hl << 32, u = u1 + l;
  /* h = hh * 2^32 + hl exactly, so the total is c * 2^64 + u. */
  int64_t hh = (h - (int64_t) hl) / INT64_C(4294967296), c = hh + (u < u1);
  int64_t w;
  if (u >= UINT64_C(0x8000000000000000)) {
    w = (int64_t) (u - UINT64_C(0x8000000000000000)) - INT64_MAX - 1; /* u - 2^64 */
    c += 1;
  } else {
    w = (int64_t) u;
  }
  rsimd_add_carry_i64(s, carry, w);
  *carry += c;
}

/* Conversions of one element (RSIMD_CVT_* rules in kernel_types.h). */
#define RSIMD_TWO53 9007199254740992.0
#define RSIMD_TWO63 9223372036854775808.0
#define RSIMD_TWO64 18446744073709551616.0
static inline double rsimd_cvt_i64_f64_1(int64_t x, int *st) {
  if (x == RSIMD_NA_I64) return rsimd_na_real();
  if (x >= (int64_t) RSIMD_TWO53 || x <= -(int64_t) RSIMD_TWO53) *st |= RSIMD_CVT_WARN_PRECISION;
  return (double) x;
}
static inline int32_t rsimd_cvt_i64_i32_1(int64_t x, int mode, int *st) {
  if (x == RSIMD_NA_I64) return RSIMD_NA_I32;
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    if (x > INT32_MAX || x < -INT32_MAX) {
      *st |= RSIMD_CVT_WARN_I32_OVF;
      return RSIMD_NA_I32;
    }
    return (int32_t) x;
  case RSIMD_CVT_SATURATING: return x > INT32_MAX ? INT32_MAX : x < -INT32_MAX ? -INT32_MAX : (int32_t) x;
  default: {
    uint32_t w = (uint32_t) ((uint64_t) x & 0xFFFFFFFFu);
    return w <= 0x7FFFFFFFu ? (int32_t) w : (int32_t) (w - 0x80000000u) - 0x7FFFFFFF - 1;
  }
  }
}
static inline uint8_t rsimd_cvt_i64_u8_1(int64_t x, int mode, int *st) {
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    if (x == RSIMD_NA_I64 || x < 0 || x > 255) {
      *st |= RSIMD_CVT_WARN_RAW;
      return 0;
    }
    return (uint8_t) x;
  case RSIMD_CVT_SATURATING:
    if (x == RSIMD_NA_I64 || x < 0) return 0;
    return x > 255 ? 255 : (uint8_t) x;
  default: return x == RSIMD_NA_I64 ? 0 : (uint8_t) ((uint64_t) x & 0xFFu);
  }
}
static inline int64_t rsimd_cvt_f64_i64_1(double x, int mode, int *st) {
  double t, w;
  if (isnan(x)) return RSIMD_NA_I64;
  t = trunc(x);
  switch (mode) {
  case RSIMD_CVT_CHECKED:
    if (t >= RSIMD_TWO63 || t <= -RSIMD_TWO63) {
      *st |= RSIMD_CVT_WARN_I64;
      return RSIMD_NA_I64;
    }
    return (int64_t) t;
  case RSIMD_CVT_SATURATING:
    if (t >= RSIMD_TWO63) return INT64_MAX;
    if (t <= -RSIMD_TWO63) return -INT64_MAX;
    return (int64_t) t;
  default:
    if (isinf(x)) return RSIMD_NA_I64;
    if (t < RSIMD_TWO63 && t > -RSIMD_TWO63) return (int64_t) t;
    /* |t| >= 2^63 is a multiple of 2^11, so w is exact, in [0, 2^64). */
    w = fmod(t, RSIMD_TWO64);
    if (w < 0) w += RSIMD_TWO64;
    if (w < RSIMD_TWO63) return (int64_t) w;
    return w == RSIMD_TWO63 ? RSIMD_NA_I64 : -(int64_t) (RSIMD_TWO64 - w);
  }
}
static inline int rsimd_cvt_i64_1(int op, int mode, const void *x, R_xlen_t i, void *out) {
  int st = 0;
  switch (op) {
  case RSIMD_CVT_I64_F64:
    ((double *) out)[i] = rsimd_cvt_i64_f64_1(((const int64_t *) x)[i], &st);
    break;
  case RSIMD_CVT_I64_I32:
    ((int *) out)[i] = rsimd_cvt_i64_i32_1(((const int64_t *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_I64_U8:
    ((Rbyte *) out)[i] = rsimd_cvt_i64_u8_1(((const int64_t *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_I64_LGL: {
    int64_t v = ((const int64_t *) x)[i];
    ((int *) out)[i] = v == RSIMD_NA_I64 ? RSIMD_NA_I32 : v != 0;
    break;
  }
  case RSIMD_CVT_F64_I64:
    ((int64_t *) out)[i] = rsimd_cvt_f64_i64_1(((const double *) x)[i], mode, &st);
    break;
  case RSIMD_CVT_I32_I64: {
    int v = ((const int *) x)[i];
    ((int64_t *) out)[i] = v == RSIMD_NA_I32 ? RSIMD_NA_I64 : (int64_t) v;
    break;
  }
  default: ((int64_t *) out)[i] = ((const Rbyte *) x)[i]; break;
  }
  return st;
}

/* The scan of one element: cumsum (op 0), cummin (2) or cummax (3) from
   s->i64; returns 0, or 1 at an NA (or, for cumsum, an overflow, which
   sets s->overflow) where the scan stops. */
static inline int rsimd_scan_i64_1(int op, int64_t x, int64_t *out, rsimd_scan_state *s) {
  int ovf = 0;
  int64_t v;
  if (x == RSIMD_NA_I64) return 1;
  if (op == 0) {
    v = rsimd_add_i64(s->i64, x, 0, &ovf);
    if (ovf) {
      s->overflow = 1;
      return 1;
    }
  } else if (op == 2) {
    v = x <= s->i64 ? x : s->i64;
  } else {
    v = x >= s->i64 ? x : s->i64;
  }
  s->i64 = v;
  *out = v;
  return 0;
}

/* pmin, pmax and their _num forms read int32 NA as NA whatever na_check
   says, as the int32 kernels do. */
static inline int rsimd_ew_i64_check(int op, const rsimd_opts *o) {
  return op >= RSIMD_EW_PMIN && op <= RSIMD_EW_PMAX_NUM ? 1 : o->na_check;
}

#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  R_xlen_t i, removed = 0;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I64) {
      r->saw_na = 1;
      if (o->na_rm) {
        removed++;
        continue;
      }
      /* The sum is NA whatever follows. */
      r->count += n;
      return;
    }
    rsimd_add_carry_i64(&r->i64, &r->carry, x[i]);
  }
  r->count += n - removed;
}

void RSIMD_KERNEL(minmax_i64)(const int64_t *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i64)(const int64_t *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  int64_t lo = r->i64, hi = r->i64_hi;
  R_xlen_t i, missing = 0;
  for (i = 0; i < n; i++) {
    /* |NA| wraps to NA. */
    int64_t v = absval ? rsimd_abs_i64(x[i]) : x[i];
    if (check && v == RSIMD_NA_I64) {
      missing++;
      r->saw_na = 1;
      /* Without na.rm the result is NA whatever follows (the scalar
         integer64 loop is not vectorised, so the test costs nothing). */
      if (!o->na_rm) {
        r->count += n;
        return;
      }
      continue;
    }
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  r->i64 = lo;
  r->i64_hi = hi;
  r->count += n - missing;
}

R_xlen_t RSIMD_KERNEL(find_i64)(const int64_t *x, R_xlen_t n, int absval, int64_t v);
R_xlen_t RSIMD_KERNEL(find_i64)(const int64_t *x, R_xlen_t n, int absval, int64_t v) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if ((absval ? rsimd_abs_i64(x[i]) : x[i]) == v) return i;
  }
  return -1;
}

void RSIMD_KERNEL(anyall_i64)(const int64_t *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_i64)(const int64_t *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  int check = o->na_check || o->na_rm;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I64) {
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

void RSIMD_KERNEL(na_i64)(const int64_t *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_i64)(const int64_t *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (x[i] != RSIMD_NA_I64) continue;
    if (mode == RSIMD_NAMODE_ANY) {
      r->saw_na = 1;
      return;
    }
    if (mode == RSIMD_NAMODE_COUNT) r->i64++;
    else RSIMD_KERNEL(put_index_)(mode, out, r, off + i);
  }
}

R_xlen_t RSIMD_KERNEL(scan_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                                rsimd_scan_state *s);
R_xlen_t RSIMD_KERNEL(scan_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                                rsimd_scan_state *s) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (rsimd_scan_i64_1(op, x[i], out + i, s)) return i;
  }
  return -1;
}

int RSIMD_KERNEL(ew1_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew1_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                          const rsimd_opts *o) {
  R_xlen_t i;
  (void) o;
  for (i = 0; i < n; i++) {
    switch (op) {
    case RSIMD_EW_NEG: out[i] = rsimd_neg_i64(x[i]); break;
    case RSIMD_EW_ABS: out[i] = rsimd_abs_i64(x[i]); break;
    default: out[i] = rsimd_sign_i64(x[i]); break;
    }
  }
  return 0;
}

int RSIMD_KERNEL(ew2_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          int64_t *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew2_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          int64_t *out, const rsimd_opts *o) {
  const int check = rsimd_ew_i64_check(op, o);
  int st = 0;
  R_xlen_t i;
  if (op == RSIMD_EW_IDIV || op == RSIMD_EW_MOD) {
    return rsimd_ew_intdiv_i64_scalar(x, y, n, flags, op == RSIMD_EW_MOD, check, out);
  }
  for (i = 0; i < n; i++) {
    out[i] = rsimd_ew2_i64_1(op, rsimd_i64_get(x, flags, 0, i, check),
                             rsimd_i64_get(y, flags, 1, i, check), check, &st);
  }
  return st;
}

int RSIMD_KERNEL(ew3_i64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, int64_t *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_i64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, int64_t *out, const rsimd_opts *o) {
  const int check = op == RSIMD_EW_CLAMP ? 1 : o->na_check;
  int st = 0;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    out[i] = rsimd_ew3_i64_1(op, rsimd_i64_get(x, flags, 0, i, check),
                             rsimd_i64_get(y, flags, 1, i, check),
                             rsimd_i64_get(z, flags, 2, i, check), check, &st);
  }
  return st;
}

void RSIMD_KERNEL(cmp_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    out[i] = rsimd_cmp_i64_1(op, rsimd_i64_get(x, flags, 0, i, 1), rsimd_i64_get(y, flags, 1, i, 1));
  }
}

int RSIMD_KERNEL(pred_i64)(int op, const int64_t *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_i64)(int op, const int64_t *x, R_xlen_t n, int mode, int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int v = rsimd_pred_i64_1(op, x[i]);
    if (mode == RSIMD_PRED_ELT) out[i] = v;
    else if (mode == RSIMD_PRED_ANY && v) return 1;
    else if (mode == RSIMD_PRED_ALL && !v) return 0;
  }
  return mode == RSIMD_PRED_ALL;
}

void RSIMD_KERNEL(bit_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags, int k,
                           void *out, const rsimd_opts *o);
void RSIMD_KERNEL(bit_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags, int k,
                           void *out, const rsimd_opts *o) {
  const int binary = op <= RSIMD_BIT_XOR, check = o->na_check, counts = rsimd_bit_counts_i64(op);
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    int64_t a = rsimd_i64_get(x, flags, 0, i, 1), b = binary ? rsimd_i64_get(y, flags, 1, i, 1) : 0;
    int na = check && rsimd_na2_i64(a, b);
    if (counts) ((int *) out)[i] = na ? RSIMD_NA_I32 : (int) rsimd_bit_i64_1(op, a, b, k);
    else ((int64_t *) out)[i] = na ? RSIMD_NA_I64 : rsimd_bit_i64_1(op, a, b, k);
  }
}

void RSIMD_KERNEL(popcnt_sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o);
void RSIMD_KERNEL(popcnt_sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I64) {
      r->saw_na = 1;
      if (!o->na_rm) return;
      continue;
    }
    r->i64 += rsimd_popcnt64_1((uint64_t) x[i]);
  }
}

void RSIMD_KERNEL(hamming_i64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_i64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  rsimd_ham_i64_from(x, y, 0, n, flags, r, o->na_rm);
}

void RSIMD_KERNEL(hamming_bits_i64)(const int64_t *x, const int64_t *y, R_xlen_t n, int flags,
                                    rsimd_reduce_result *r);
void RSIMD_KERNEL(hamming_bits_i64)(const int64_t *x, const int64_t *y, R_xlen_t n, int flags,
                                    rsimd_reduce_result *r) {
  rsimd_ham_bits_i64_from(x, y, 0, n, flags, r);
}

/* The integer64 conversions of the convert slot (convert.inc.c). */
RSIMD_INLINE int RSIMD_KERNEL(convert_i64_)(int op, int mode, const void *x, R_xlen_t n, void *out) {
  int st = 0;
  R_xlen_t i;
  for (i = 0; i < n; i++) st |= rsimd_cvt_i64_1(op, mode, x, i, out);
  return st;
}

#else /* vector tiers */

/* cumsum, cummin and cummax are sequential: the none tier's kernel is used. */
#define RSIMD_SKIP_scan_i64 1

/* ---- Lane helpers ---------------------------------------------------------- */

RSIMD_INLINE rsimd_mi64 rsimd_vi64_na2(rsimd_vi64 x, rsimd_vi64 y) {
  return rsimd_mi64_or(rsimd_vi64_is_na(x), rsimd_vi64_is_na(y));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_set_na(rsimd_vi64 r, rsimd_mi64 m) {
  return rsimd_vi64_blend(r, rsimd_vi64_set1(RSIMD_NA_I64), m);
}
/* Logical lanes from 64-bit masks: 1 where t, 0 elsewhere, NA_integer_
   where na (stored as int32 elements with rsimd_vi64_storeu_i32). */
RSIMD_INLINE rsimd_vi64 rsimd_lgl_vi64m(rsimd_mi64 t, rsimd_mi64 na) {
  rsimd_vi64 r = rsimd_vi64_blend(rsimd_vi64_zero(), rsimd_vi64_set1(1), t);
  return rsimd_vi64_blend(r, rsimd_vi64_set1(RSIMD_NA_I32), na);
}
/* The overflow lanes of r = a + b and r = a - b (wrapped), INT64_MIN
   results included. */
RSIMD_INLINE rsimd_mi64 rsimd_vi64_add_ovf(rsimd_vi64 a, rsimd_vi64 b, rsimd_vi64 r) {
  rsimd_vi64 t = rsimd_vi64_and(rsimd_vi64_xor(a, r), rsimd_vi64_xor(b, r));
  return rsimd_mi64_or(rsimd_vi64_cmp_lt(t, rsimd_vi64_zero()), rsimd_vi64_is_na(r));
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_sub_ovf(rsimd_vi64 a, rsimd_vi64 b, rsimd_vi64 r) {
  rsimd_vi64 t = rsimd_vi64_and(rsimd_vi64_xor(a, b), rsimd_vi64_xor(a, r));
  return rsimd_mi64_or(rsimd_vi64_cmp_lt(t, rsimd_vi64_zero()), rsimd_vi64_is_na(r));
}
/* The lanes in (-2^31, 2^31), whose products cannot overflow (NA is not
   among them). */
RSIMD_INLINE rsimd_mi64 rsimd_vi64_small(rsimd_vi64 a) {
  return rsimd_mi64_and(rsimd_vi64_cmp_gt(a, rsimd_vi64_set1(-INT64_C(2147483648))),
                        rsimd_vi64_cmp_lt(a, rsimd_vi64_set1(INT64_C(2147483648))));
}

/* Population count and leading zeros of 64-bit lanes. */
#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
RSIMD_INLINE rsimd_vi64 rsimd_vi64_popcnt(rsimd_vi64 x) {
  return svreinterpret_s64_u64(svcnt_s64_x(RSIMD_PT64, x));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_clz(rsimd_vi64 x) {
  return svreinterpret_s64_u64(svclz_s64_x(RSIMD_PT64, x));
}
#else
#if RSIMD_TIER_IS(neon) && defined(SIMDE_ARM_NEON_A32V7_NATIVE)
/* Byte counts summed pairwise into 64-bit lanes. */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_popcnt(rsimd_vi64 x) {
  uint8x16_t c = vcntq_u8(vreinterpretq_u8_u64(simde__m128i_to_neon_u64(x)));
  return simde__m128i_from_neon_u64(vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(c))));
}
#else
/* SWAR population count of each byte, then the bytes of a lane summed:
   psadbw against zero on sse2 and avx2, shifts on avx512 (the tier's
   subset has neither a 512-bit psadbw nor VPOPCNTQ). */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_popcnt(rsimd_vi64 x) {
  const rsimd_vi64 m1 = rsimd_vi64_set1(INT64_C(0x5555555555555555)),
                   m2 = rsimd_vi64_set1(INT64_C(0x3333333333333333)),
                   m4 = rsimd_vi64_set1(INT64_C(0x0F0F0F0F0F0F0F0F));
  x = rsimd_vi64_sub(x, rsimd_vi64_and(rsimd_vi64_srl(x, 1), m1));
  x = rsimd_vi64_add(rsimd_vi64_and(x, m2), rsimd_vi64_and(rsimd_vi64_srl(x, 2), m2));
  x = rsimd_vi64_and(rsimd_vi64_add(x, rsimd_vi64_srl(x, 4)), m4);
#if RSIMD_TIER_IS(avx2)
  return simde_mm256_sad_epu8(x, simde_mm256_setzero_si256());
#elif RSIMD_TIER_IS(avx512)
  x = rsimd_vi64_add(x, rsimd_vi64_srl(x, 8));
  x = rsimd_vi64_add(x, rsimd_vi64_srl(x, 16));
  x = rsimd_vi64_add(x, rsimd_vi64_srl(x, 32));
  return rsimd_vi64_and(x, rsimd_vi64_set1(0x7F));
#else
  return simde_mm_sad_epu8(x, simde_mm_setzero_si128());
#endif
}
#endif
/* Smear the highest set bit down, then 64 - popcount. */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_clz(rsimd_vi64 x) {
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 1));
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 2));
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 4));
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 8));
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 16));
  x = rsimd_vi64_or(x, rsimd_vi64_srl(x, 32));
  return rsimd_vi64_sub(rsimd_vi64_set1(64), rsimd_vi64_popcnt(x));
}
#endif
/* Trailing zeros: the bits below the lowest set bit, ~x & (x - 1). */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_ctz(rsimd_vi64 x) {
  return rsimd_vi64_popcnt(rsimd_vi64_andnot(x, rsimd_vi64_sub(x, rsimd_vi64_set1(1))));
}

/* Operand k of an elementwise kernel at i, as int64 lanes: the broadcast
   bc, int64 elements, or int32 elements sign-extended with NA_integer_
   becoming NA (with check). The _p form reads the active lanes of pg and
   repeats element i in the others. */
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_i64_from_i32(rsimd_vi64 v, int check) {
  if (!check) return v;
  return rsimd_vi64_set_na(v, rsimd_vi64_cmp_eq(v, rsimd_vi64_set1(RSIMD_NA_I32)));
}
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_i64_ld(const void *p, int flags, int k, rsimd_vi64 bc,
                                            ptrdiff_t i, int check) {
  if (flags & RSIMD_EW_SCALAR(k)) return bc;
  if (flags & RSIMD_EW_I32(k)) {
    return rsimd_i64_from_i32(rsimd_vi64_loadu_i32((const int32_t *) p + i), check);
  }
  return rsimd_vi64_loadu((const int64_t *) p + i);
}
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_i64_ld_p(const void *p, int flags, int k, rsimd_vi64 bc,
                                              ptrdiff_t i, rsimd_p64 pg, int check) {
  if (flags & RSIMD_EW_SCALAR(k)) return bc;
  if (flags & RSIMD_EW_I32(k)) {
    const int32_t *q = (const int32_t *) p + i;
    return rsimd_i64_from_i32(rsimd_vi64_loadu_i32_p(pg, q, q[0]), check);
  }
  return rsimd_vi64_loadu_p(pg, (const int64_t *) p + i, ((const int64_t *) p)[i]);
}
/* Operand k broadcast from its element 0 (zeros for a NULL operand). */
static inline rsimd_vi64 rsimd_i64_bcast(const void *p, int flags, int k, int check) {
  if (p == NULL) return rsimd_vi64_zero();
  return rsimd_vi64_set1(rsimd_i64_get(p, flags, k, 0, check));
}

/* Runs `expr`, which sets the rsimd_vi64 r from a and b (and c), over the
   chunk with operand flags fl, storing r at out with ST / ST_P (int64 or
   int32 elements of type OT); `lanes` is the number of lanes computed. */
#define RSIMD_I64_LOOP_(fl, nargs, OT, ST, ST_P, expr)                                     \
  do {                                                                                     \
    ptrdiff_t i = 0;                                                                       \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                                 \
      rsimd_vi64 a = rsimd_i64_ld(x, fl, 0, bc0, i, check), b = bc1, c = bc2, r;           \
      if ((nargs) > 1) b = rsimd_i64_ld(y, fl, 1, bc1, i, check);                          \
      if ((nargs) > 2) c = rsimd_i64_ld(z, fl, 2, bc2, i, check);                          \
      (void) b;                                                                            \
      (void) c;                                                                            \
      expr;                                                                                \
      ST((OT *) out + i, r);                                                               \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                                \
      rsimd_vi64 a = rsimd_i64_ld_p(x, fl, 0, bc0, i, pg, check), b = bc1, c = bc2, r;     \
      if ((nargs) > 1) b = rsimd_i64_ld_p(y, fl, 1, bc1, i, pg, check);                    \
      if ((nargs) > 2) c = rsimd_i64_ld_p(z, fl, 2, bc2, i, pg, check);                    \
      (void) b;                                                                            \
      (void) c;                                                                            \
      expr;                                                                                \
      ST_P(pg, (OT *) out + i, r);                                                         \
    }                                                                                      \
  } while (0)
/* The loop with the flags as constants in the two common cases (all
   operands int64 vectors; an int64 scalar y). */
#define RSIMD_I64_LOOP_T(nargs, OT, ST, ST_P, expr)                                        \
  do {                                                                                     \
    if (flags == 0) {                                                                      \
      RSIMD_I64_LOOP_(0, nargs, OT, ST, ST_P, expr);                                       \
    } else if (flags == RSIMD_EW_SCALAR(1)) {                                              \
      RSIMD_I64_LOOP_(RSIMD_EW_SCALAR(1), nargs, OT, ST, ST_P, expr);                      \
    } else {                                                                               \
      RSIMD_I64_LOOP_(flags, nargs, OT, ST, ST_P, expr);                                   \
    }                                                                                      \
  } while (0)
#define RSIMD_I64_LOOP(nargs, expr)                                                        \
  RSIMD_I64_LOOP_T(nargs, int64_t, rsimd_vi64_storeu, rsimd_vi64_storeu_p, expr)
#define RSIMD_I64_LOOP32(nargs, expr)                                                      \
  RSIMD_I64_LOOP_T(nargs, int32_t, rsimd_vi64_storeu_i32, rsimd_vi64_storeu_i32_p, expr)

#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
/* Checked product of a and b from the high half of the product. */
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_i64_mul_checked(rsimd_vi64 a, rsimd_vi64 b, int check,
                                                     int *st) {
  rsimd_vi64 r = svmul_s64_x(RSIMD_PT64, a, b), hi = svmulh_s64_x(RSIMD_PT64, a, b);
  rsimd_mi64 na = check ? rsimd_vi64_na2(a, b) : rsimd_mi64_none();
  rsimd_mi64 o = rsimd_mi64_or(rsimd_mi64_not(rsimd_vi64_cmp_eq(hi, rsimd_vi64_sign(r))),
                               rsimd_vi64_is_na(r));
  if (rsimd_mi64_any(rsimd_mi64_andnot(na, o))) *st |= RSIMD_EW_OVERFLOW;
  return rsimd_vi64_set_na(r, rsimd_mi64_or(o, na));
}
#else
/* Checked products: exact in vectors where both operands of every lane are
   small; a vector that has a larger one takes the scalar product lane by
   lane, and the loop goes on in vectors, until 8 such vectors in a row,
   after which the rest of the chunk is scalar (testing every vector is
   slower than the scalar loop when most operands are large). */
RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(mul_i64_)(const void *x, const void *y, R_xlen_t n,
                                              const int flags, int64_t *out, const int check) {
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, check), bc1 = rsimd_i64_bcast(y, flags, 1, check);
  int st = 0, large = 0;
  ptrdiff_t i = 0, j;
  for (; i + RSIMD_LANES_64 <= n && large < 8; i += RSIMD_LANES_64) {
    rsimd_vi64 a = rsimd_i64_ld(x, flags, 0, bc0, i, check), b = rsimd_i64_ld(y, flags, 1, bc1, i, check);
    if (rsimd_mi64_all(rsimd_mi64_and(rsimd_vi64_small(a), rsimd_vi64_small(b)))) {
      rsimd_vi64_storeu(out + i, rsimd_vi64_mul32(a, b));
      large = 0;
      continue;
    }
    large++;
    for (j = i; j < i + RSIMD_LANES_64; j++) {
      int ovf = 0;
      out[j] = rsimd_mul_i64(rsimd_i64_get(x, flags, 0, j, check),
                             rsimd_i64_get(y, flags, 1, j, check), check, &ovf);
      if (ovf) st = RSIMD_EW_OVERFLOW;
    }
  }
  for (; i < n; i++) {
    int ovf = 0;
    out[i] = rsimd_mul_i64(rsimd_i64_get(x, flags, 0, i, check), rsimd_i64_get(y, flags, 1, i, check),
                           check, &ovf);
    if (ovf) st = RSIMD_EW_OVERFLOW;
  }
  return st;
}
#endif

/* ---- Reductions -------------------------------------------------------------- */

/* The sum splits every element v into its high word, biased to the
   unsigned (v >> 32) ^ 2^31 = hi(v) + 2^31, and its low word, unsigned, and
   adds them in separate 64-bit lanes: five native instructions per vector
   on every tier (no 64-bit compares, no arithmetic shift). A lane that has
   taken m elements holds at most m * 2^32 in each, so the lanes are moved
   into the total every RSIMD_SUM_I64_FLUSH vectors, long before they could
   overflow; a lane's elements total (hi - m * 2^31) * 2^32 + lo. The
   total (r->i64, r->carry) represents each exact value one way only, so
   every tier gives the none tier's result. */
#define RSIMD_SUM_I64_FLUSH ((ptrdiff_t) 1 << (sizeof(ptrdiff_t) > 4 ? 30 : 24))

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(sum_i64_)(const int64_t *x, R_xlen_t n, const int check,
                                               const int narm, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vi64 zero = rsimd_vi64_zero(), bias = rsimd_vi64_set1(INT64_C(2147483648)),
                   lo_mask = rsimd_vi64_set1(INT64_C(0xFFFFFFFF));
  rsimd_mi64 mna = rsimd_mi64_none();
  rsimd_vi64 rcnt = zero;
  int64_t lh[RSIMD_MAX_LANES_64], ll[RSIMD_MAX_LANES_64];
  ptrdiff_t i = 0, j;
  /* Without na.rm, blocks of RSIMD_FOLD_BLOCK elements, after each of which
     an NA ends the sum (it is then NA whatever follows). */
  const ptrdiff_t block = check && !narm ? RSIMD_FOLD_BLOCK : 2 * W * RSIMD_SUM_I64_FLUSH;
  while (i < n) {
    rsimd_vi64 h0 = zero, h1 = zero, l0 = zero, l1 = zero;
    ptrdiff_t m = 0, end = n - i > block ? i + block : n;
#define RSIMD_SUM_I64_(h, l, v)                                                            \
  do {                                                                                     \
    rsimd_vi64 v_ = (v);                                                                   \
    if (check) {                                                                           \
      rsimd_mi64 m_ = rsimd_vi64_is_na(v_);                                                \
      mna = rsimd_mi64_or(mna, m_);                                                        \
      if (narm) {                                                                          \
        v_ = rsimd_vi64_blend(v_, zero, m_);                                               \
        rcnt = rsimd_vi64_inc(rcnt, m_);                                                   \
      }                                                                                    \
    }                                                                                      \
    (h) = rsimd_vi64_add((h), rsimd_vi64_xor(rsimd_vi64_srl(v_, 32), bias));               \
    (l) = rsimd_vi64_add((l), rsimd_vi64_and(v_, lo_mask));                                \
  } while (0)
    for (; i + 2 * W <= end; i += 2 * W, m++) {
      RSIMD_SUM_I64_(h0, l0, rsimd_vi64_loadu(x + i));
      RSIMD_SUM_I64_(h1, l1, rsimd_vi64_loadu(x + i + W));
    }
    /* The tail, in h1 and l1: the fill 0 adds the bias only (and is not NA);
       a missing second vector adds only the bias too. */
    if (i < end) {
      RSIMD_SUM_I64_(h0, l0, rsimd_vi64_loadu_p(rsimd_p64_while(i, end), x + i, 0));
      RSIMD_SUM_I64_(h1, l1, rsimd_vi64_loadu_p(rsimd_p64_while(i + W, end), x + i + W, 0));
      i = end;
      m++;
    }
#undef RSIMD_SUM_I64_
    /* Each lane of h0 and h1 took m biased high words. */
    rsimd_vi64_storeu(lh, h0);
    rsimd_vi64_storeu(ll, l0);
    for (j = 0; j < W; j++) {
      rsimd_add_hl_i64(&r->i64, &r->carry, lh[j] - m * INT64_C(2147483648), (uint64_t) ll[j]);
    }
    rsimd_vi64_storeu(lh, h1);
    rsimd_vi64_storeu(ll, l1);
    for (j = 0; j < W; j++) {
      rsimd_add_hl_i64(&r->i64, &r->carry, lh[j] - m * INT64_C(2147483648), (uint64_t) ll[j]);
    }
    if (check && !narm && rsimd_mi64_any(mna)) {
      r->saw_na = 1;
      r->count += n;
      return;
    }
  }
  if (rsimd_mi64_any(mna)) r->saw_na = 1;
  r->count += n - rsimd_vi64_reduce_add(rcnt);
}
#undef RSIMD_SUM_I64_FLUSH

void RSIMD_KERNEL(sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o);
void RSIMD_KERNEL(sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                           const rsimd_opts *o) {
  if (o->na_rm) RSIMD_KERNEL(sum_i64_)(x, n, 1, 1, r);
  else if (o->na_check) RSIMD_KERNEL(sum_i64_)(x, n, 1, 0, r);
  else RSIMD_KERNEL(sum_i64_)(x, n, 0, 0, r);
}

/* The number of NA elements among the n, counted in the lanes of four
   accumulators (mask counts are slow on the 128-bit tiers). */
RSIMD_INLINE R_xlen_t RSIMD_KERNEL(count_na_i64_)(const int64_t *x, R_xlen_t n) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vi64 zero = rsimd_vi64_zero();
  rsimd_vi64 c0 = zero, c1 = zero, c2 = zero, c3 = zero;
  R_xlen_t i = 0, k;
  for (; i + 4 * W <= n; i += 4 * W) {
    c0 = rsimd_vi64_inc(c0, rsimd_vi64_is_na(rsimd_vi64_loadu(x + i)));
    c1 = rsimd_vi64_inc(c1, rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + W)));
    c2 = rsimd_vi64_inc(c2, rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + 2 * W)));
    c3 = rsimd_vi64_inc(c3, rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + 3 * W)));
  }
  for (; i + W <= n; i += W) c0 = rsimd_vi64_inc(c0, rsimd_vi64_is_na(rsimd_vi64_loadu(x + i)));
  k = rsimd_vi64_reduce_add(rsimd_vi64_add(rsimd_vi64_add(c0, c1), rsimd_vi64_add(c2, c3)));
  for (; i < n; i++) k += x[i] == RSIMD_NA_I64;
  return k;
}

/* Minimum and/or maximum (ext) of the n integer64 elements, or of their
   absolute values, into r. With stop (na_check without na.rm) a block of
   RSIMD_FOLD_BLOCK elements with an NA ends the scan, as the result is
   then NA. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(minmax_i64_)(const int64_t *x, R_xlen_t n, const int check,
                                                  const int stop, const int absval,
                                                  const int ext, rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int want_lo = ext & RSIMD_EXT_MIN, want_hi = ext & RSIMD_EXT_MAX;
  const rsimd_vi64 big = rsimd_vi64_set1(INT64_MAX), zero = rsimd_vi64_zero();
  rsimd_vi64 lo0 = big, lo1 = big, hi0 = rsimd_vi64_set1(INT64_MIN), hi1 = hi0;
  rsimd_mi64 mna = rsimd_mi64_none();
  R_xlen_t i = 0, missing = 0;
  /* NA (INT64_MIN) lanes cannot raise the maximum; for the minimum they
     are replaced by INT64_MAX. The absolute value max(v, -v) of NA wraps
     to NA. */
#define RSIMD_MINMAX_I64_(lo, hi, v0)                                                    \
  do {                                                                                   \
    rsimd_vi64 v_ = absval ? rsimd_vi64_max((v0), rsimd_vi64_sub(zero, (v0))) : (v0);    \
    rsimd_vi64 m_ = v_;                                                                  \
    if (check) {                                                                         \
      rsimd_mi64 na_ = rsimd_vi64_is_na(v_);                                             \
      mna = rsimd_mi64_or(mna, na_);                                                     \
      if (want_lo) m_ = rsimd_vi64_blend(m_, big, na_);                                  \
    }                                                                                    \
    if (want_lo) (lo) = rsimd_vi64_min((lo), m_);                                        \
    if (want_hi) (hi) = rsimd_vi64_max((hi), v_);                                        \
  } while (0)
  while (i + 2 * W <= n) {
    const R_xlen_t end = stop && n - i > RSIMD_FOLD_BLOCK ? i + RSIMD_FOLD_BLOCK : n;
    for (; i + 2 * W <= end; i += 2 * W) {
      rsimd_vi64 v0 = rsimd_vi64_loadu(x + i), v1 = rsimd_vi64_loadu(x + i + W);
      RSIMD_MINMAX_I64_(lo0, hi0, v0);
      RSIMD_MINMAX_I64_(lo1, hi1, v1);
    }
    if (stop && rsimd_mi64_any(mna)) {
      r->saw_na = 1;
      r->count += n;
      return;
    }
  }
  for (; i < n; i += W) {
    /* The tail repeats element i in the inactive lanes. */
    rsimd_vi64 v = rsimd_vi64_loadu_p(rsimd_p64_while(i, n), x + i, x[i]);
    RSIMD_MINMAX_I64_(lo0, hi0, v);
  }
#undef RSIMD_MINMAX_I64_
  if (check && rsimd_mi64_any(mna)) {
    r->saw_na = 1;
    if (stop) {
      r->count += n;
      return;
    }
    missing = RSIMD_KERNEL(count_na_i64_)(x, n);
  }
  if (missing < n) {
    if (want_lo) {
      int64_t lo = rsimd_vi64_reduce_min(rsimd_vi64_min(lo0, lo1));
      if (lo < r->i64) r->i64 = lo;
    }
    if (want_hi) {
      int64_t hi = rsimd_vi64_reduce_max(rsimd_vi64_max(hi0, hi1));
      if (hi > r->i64_hi) r->i64_hi = hi;
    }
  }
  r->count += n - missing;
}

void RSIMD_KERNEL(minmax_i64)(const int64_t *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(minmax_i64)(const int64_t *x, R_xlen_t n, int absval, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
#define RSIMD_MINMAX_I64_EXT_(c, s, a)                                                   \
  do {                                                                                   \
    if (o->extrema == RSIMD_EXT_MIN) RSIMD_KERNEL(minmax_i64_)(x, n, c, s, a, RSIMD_EXT_MIN, r); \
    else if (o->extrema == RSIMD_EXT_MAX) RSIMD_KERNEL(minmax_i64_)(x, n, c, s, a, RSIMD_EXT_MAX, r); \
    else RSIMD_KERNEL(minmax_i64_)(x, n, c, s, a, RSIMD_EXT_BOTH, r);                   \
  } while (0)
  if (!check) {
    if (absval) RSIMD_MINMAX_I64_EXT_(0, 0, 1);
    else RSIMD_MINMAX_I64_EXT_(0, 0, 0);
  } else if (o->na_rm) {
    if (absval) RSIMD_MINMAX_I64_EXT_(1, 0, 1);
    else RSIMD_MINMAX_I64_EXT_(1, 0, 0);
  } else {
    if (absval) RSIMD_MINMAX_I64_EXT_(1, 1, 1);
    else RSIMD_MINMAX_I64_EXT_(1, 1, 0);
  }
#undef RSIMD_MINMAX_I64_EXT_
}

RSIMD_ALWAYS_INLINE R_xlen_t RSIMD_KERNEL(find_i64_)(const int64_t *x, R_xlen_t n,
                                                    const int absval, int64_t v) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const rsimd_vi64 vv = rsimd_vi64_set1(v), zero = rsimd_vi64_zero();
  R_xlen_t i = 0;
#define RSIMD_FIND_I64_(j)                                                               \
  rsimd_vi64_cmp_eq(absval ? rsimd_vi64_max(rsimd_vi64_loadu(x + i + (j) * W),           \
                                            rsimd_vi64_sub(zero, rsimd_vi64_loadu(x + i + (j) * W))) \
                           : rsimd_vi64_loadu(x + i + (j) * W),                          \
                    vv)
  for (; i + 4 * W <= n; i += 4 * W) {
    rsimd_mi64 m = rsimd_mi64_or(rsimd_mi64_or(RSIMD_FIND_I64_(0), RSIMD_FIND_I64_(1)),
                                 rsimd_mi64_or(RSIMD_FIND_I64_(2), RSIMD_FIND_I64_(3)));
    if (rsimd_mi64_any(m)) break;
  }
#undef RSIMD_FIND_I64_
  for (; i < n; i++) {
    if ((absval ? rsimd_abs_i64(x[i]) : x[i]) == v) return i;
  }
  return -1;
}

R_xlen_t RSIMD_KERNEL(find_i64)(const int64_t *x, R_xlen_t n, int absval, int64_t v);
R_xlen_t RSIMD_KERNEL(find_i64)(const int64_t *x, R_xlen_t n, int absval, int64_t v) {
  if (absval) return RSIMD_KERNEL(find_i64_)(x, n, 1, v);
  return RSIMD_KERNEL(find_i64_)(x, n, 0, v);
}

/* Whole vectors give the TRUE, FALSE and NA flags; the tail is scalar.
   Early exit is checked after every vector. */
void RSIMD_KERNEL(anyall_i64)(const int64_t *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o);
void RSIMD_KERNEL(anyall_i64)(const int64_t *x, R_xlen_t n, int stop, rsimd_reduce_result *r,
                              const rsimd_opts *o) {
  const ptrdiff_t W = RSIMD_LANES_64;
  const int check = o->na_check || o->na_rm;
  const rsimd_vi64 zero = rsimd_vi64_zero();
  R_xlen_t i = 0;
  for (; i + W <= n; i += W) {
    rsimd_vi64 v = rsimd_vi64_loadu(x + i);
    rsimd_mi64 na = check ? rsimd_vi64_is_na(v) : rsimd_mi64_none();
    rsimd_mi64 z = rsimd_vi64_cmp_eq(v, zero);
    if (rsimd_mi64_any(na)) r->saw_na = 1;
    if (rsimd_mi64_any(z)) r->any_false = 1;
    if (rsimd_mi64_any(rsimd_mi64_not(rsimd_mi64_or(na, z)))) r->any_true = 1;
    if ((stop == RSIMD_STOP_TRUE && r->any_true) || (stop == RSIMD_STOP_FALSE && r->any_false)) {
      return;
    }
  }
  for (; i < n; i++) {
    if (check && x[i] == RSIMD_NA_I64) {
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

void RSIMD_KERNEL(na_i64)(const int64_t *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r);
void RSIMD_KERNEL(na_i64)(const int64_t *x, R_xlen_t n, int mode, R_xlen_t off, void *out,
                          rsimd_reduce_result *r) {
  const ptrdiff_t W = RSIMD_LANES_64;
  R_xlen_t i = 0, j;
  if (mode == RSIMD_NAMODE_ANY) {
    for (; i + 4 * W <= n; i += 4 * W) {
      rsimd_mi64 m = rsimd_mi64_or(rsimd_mi64_or(rsimd_vi64_is_na(rsimd_vi64_loadu(x + i)),
                                                 rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + W))),
                                   rsimd_mi64_or(rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + 2 * W)),
                                                 rsimd_vi64_is_na(rsimd_vi64_loadu(x + i + 3 * W))));
      if (rsimd_mi64_any(m)) break;
    }
    for (; i < n; i++) {
      if (x[i] == RSIMD_NA_I64) {
        r->saw_na = 1;
        return;
      }
    }
  } else if (mode == RSIMD_NAMODE_COUNT) {
    r->i64 += RSIMD_KERNEL(count_na_i64_)(x, n);
  } else {
    for (; i + W <= n; i += W) {
      if (!rsimd_mi64_any(rsimd_vi64_is_na(rsimd_vi64_loadu(x + i)))) continue;
      for (j = i; j < i + W; j++) {
        if (x[j] == RSIMD_NA_I64) RSIMD_KERNEL(put_index_)(mode, out, r, off + j);
      }
    }
    for (; i < n; i++) {
      if (x[i] == RSIMD_NA_I64) RSIMD_KERNEL(put_index_)(mode, out, r, off + i);
    }
  }
}

/* ---- Elementwise --------------------------------------------------------------- */

int RSIMD_KERNEL(ew1_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                          const rsimd_opts *o);
int RSIMD_KERNEL(ew1_i64)(int op, const int64_t *x, R_xlen_t n, int64_t *out,
                          const rsimd_opts *o) {
  const void *y = NULL, *z = NULL;
  const int flags = 0, check = 1;
  const rsimd_vi64 bc0 = rsimd_vi64_zero(), bc1 = bc0, bc2 = bc0, zero = bc0;
  (void) o;
  (void) y;
  (void) z;
  switch (op) {
  case RSIMD_EW_NEG: RSIMD_I64_LOOP(1, r = rsimd_vi64_sub(zero, a)); break;
  case RSIMD_EW_ABS: RSIMD_I64_LOOP(1, r = rsimd_vi64_max(a, rsimd_vi64_sub(zero, a))); break;
  default:
    RSIMD_I64_LOOP(1, {
      r = rsimd_vi64_blend(zero, rsimd_vi64_set1(1), rsimd_vi64_cmp_gt(a, zero));
      r = rsimd_vi64_blend(r, rsimd_vi64_set1(-1), rsimd_vi64_cmp_lt(a, zero));
      r = rsimd_vi64_set_na(r, rsimd_vi64_is_na(a));
    });
    break;
  }
  return 0;
}

#if (RSIMD_TIER_IS(neon) && defined(SIMDE_ARM_NEON_A64V8_NATIVE)) || RSIMD_TIER_IS(sve) || \
  RSIMD_TIER_IS(sve2)
/* Doubles truncated to int64 lanes, saturating (NaN gives 0). */
RSIMD_INLINE rsimd_vi64 rsimd_i64_cvtt_f64(rsimd_vf64 v) {
#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
  return svcvt_s64_f64_x(RSIMD_PT64, v);
#else
  return simde__m128i_from_neon_i64(vcvtq_s64_f64(simde__m128d_to_neon_f64(v)));
#endif
}
#define RSIMD_HAVE_I64_CVTT 1
#endif

#ifndef RSIMD_NO_F64_SIMD
/* %/% and %% through double lanes. For |x|, |y| <= 2^51 (y != 0) the
   rounded quotient x / y lies on the same side of every integer as the
   exact one (which is at least 1/|y| from the integers it is not, a
   relative 1/|x| > 2^-52, more than the rounding error), so its floor is
   x %/% y, and x - floor * y is exact. Each block of
   RSIMD_I64_DIV_BLOCK elements is computed so, then again with
   rsimd_intdiv_i64_loop if it had an operand outside [-2^51, 2^51), NA
   included, or a zero divisor. */
#define RSIMD_I64_DIV_BLOCK 256
/* int64 lanes in [-2^51, 2^51] to double, and integer doubles of that
   range to int64 lanes: native where the tier converts (avx512 rounds,
   which is exact for integers), else through the mantissa of 1.5 * 2^52. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_vi64_to_vf64_51(rsimd_vi64 a) {
#if defined(RSIMD_HAVE_VI64_TO_VF64) || RSIMD_TIER_IS(avx512) || RSIMD_TIER_IS(sve) || \
  RSIMD_TIER_IS(sve2)
  return rsimd_vi64_to_vf64(a);
#else
  const rsimd_vi64 k = rsimd_vi64_set1(INT64_C(0x4338000000000000));
  return rsimd_vf64_sub(rsimd_vi64_as_vf64(rsimd_vi64_add(a, k)), rsimd_vf64_set1(0x1.8p52));
#endif
}
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_vf64_to_vi64_51(rsimd_vf64 v) {
#if defined(RSIMD_HAVE_I64_CVTT)
  return rsimd_i64_cvtt_f64(v);
#elif RSIMD_TIER_IS(avx512)
  return simde_mm512_cvtpd_epi64(v);
#else
  return rsimd_vi64_sub(rsimd_vf64_as_vi64(rsimd_vf64_add(v, rsimd_vf64_set1(0x1.8p52))),
                        rsimd_vi64_set1(INT64_C(0x4338000000000000)));
#endif
}
/* x %/% y or x %% y of lanes in range (above), from their doubles. */
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_vf64_intdiv_i64(rsimd_vf64 a, rsimd_vf64 b, int mod) {
  rsimd_vf64 q = rsimd_vf64_floor(rsimd_vf64_div(a, b));
  if (mod) q = rsimd_vf64_sub(a, rsimd_vf64_mul(q, b));
  return rsimd_vf64_to_vi64_51(q);
}
/* Elements [i, e) of the operands through double lanes into out (with
   compute; else only tested). Returns 1 if the block has an operand
   outside [-2^51, 2^51) (NA included, which is INT64_MIN) or a zero
   divisor: the lanes of an OR of the operands biased by 2^51 that have a
   bit at 2^52 or above. */
RSIMD_ALWAYS_INLINE int rsimd_intdiv_i64_block(const void *x, const void *y, ptrdiff_t i, ptrdiff_t e,
                                               const int flags, const int mod, const int check,
                                               const int compute, rsimd_vi64 bc0, rsimd_vi64 bc1,
                                               rsimd_vf64 bd, int64_t *out) {
  const int ys = (flags & RSIMD_EW_SCALAR(1)) != 0;
  const rsimd_vi64 bias = rsimd_vi64_set1(INT64_C(1) << 51), zero = rsimd_vi64_zero();
  rsimd_vi64 acc = zero;
  rsimd_mi64 dz = rsimd_mi64_none();
  for (; i + RSIMD_LANES_64 <= e; i += RSIMD_LANES_64) {
    rsimd_vi64 a = rsimd_i64_ld(x, flags, 0, bc0, i, check), b = bc1;
    acc = rsimd_vi64_or(acc, rsimd_vi64_add(a, bias));
    if (!ys) {
      b = rsimd_i64_ld(y, flags, 1, bc1, i, check);
      acc = rsimd_vi64_or(acc, rsimd_vi64_add(b, bias));
      dz = rsimd_mi64_or(dz, rsimd_vi64_cmp_eq(b, zero));
    }
    if (compute) {
      rsimd_vi64_storeu(out + i, rsimd_vf64_intdiv_i64(rsimd_vi64_to_vf64_51(a),
                                                       ys ? bd : rsimd_vi64_to_vf64_51(b), mod));
    }
  }
  if (i < e) {
    rsimd_p64 pg = rsimd_p64_while(i, e);
    rsimd_vi64 a = rsimd_i64_ld_p(x, flags, 0, bc0, i, pg, check), b = bc1;
    acc = rsimd_vi64_or(acc, rsimd_vi64_add(a, bias));
    if (!ys) {
      b = rsimd_i64_ld_p(y, flags, 1, bc1, i, pg, check);
      acc = rsimd_vi64_or(acc, rsimd_vi64_add(b, bias));
      dz = rsimd_mi64_or(dz, rsimd_vi64_cmp_eq(b, zero));
    }
    if (compute) {
      rsimd_vi64_storeu_p(pg, out + i, rsimd_vf64_intdiv_i64(rsimd_vi64_to_vf64_51(a),
                                                             ys ? bd : rsimd_vi64_to_vf64_51(b), mod));
    }
  }
  return rsimd_mi64_any(rsimd_vi64_cmp_gt(rsimd_vi64_srl(acc, 52), zero)) || rsimd_mi64_any(dz);
}
/* The kernel, inlined with constant flags, mod and check; with a scalar
   y, the caller has checked that y is in range and not 0, and g is its
   divider (or NULL) for the blocks done by rsimd_intdiv_i64_loop. After
   such a block the next is tested before it is computed, so that data
   mostly out of range is not computed twice. */
RSIMD_ALWAYS_INLINE int rsimd_ew_intdiv_i64_(const void *x, const void *y, R_xlen_t n, const int flags,
                                             const int mod, const int check,
                                             const rsimd_divmagic_i64 *g, int64_t *out) {
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, check), bc1 = rsimd_i64_bcast(y, flags, 1, check);
  const rsimd_vf64 bd = rsimd_vi64_to_vf64_51(bc1);
  int st = 0, scalar = 0;
  ptrdiff_t i = 0;
  for (; i < n && scalar < 8; i += RSIMD_I64_DIV_BLOCK) {
    const ptrdiff_t e = n - i > RSIMD_I64_DIV_BLOCK ? i + RSIMD_I64_DIV_BLOCK : n;
    if (scalar && !rsimd_intdiv_i64_block(x, y, i, e, flags, mod, check, 0, bc0, bc1, bd, out)) {
      scalar = 0;
    }
    if (!scalar && !rsimd_intdiv_i64_block(x, y, i, e, flags, mod, check, 1, bc0, bc1, bd, out)) {
      continue;
    }
    scalar++;
    st |= rsimd_intdiv_i64_run(x, y, i, e, flags, mod, check, g, out);
  }
  if (i < n) st |= rsimd_intdiv_i64_run(x, y, i, n, flags, mod, check, g, out);
  return st;
}
#endif

/* %/% (mod = 0) or %% (mod = 1) of int64 operands. */
static int rsimd_ew_intdiv_i64(const void *x, const void *y, R_xlen_t n, int flags, int mod,
                               int check, int64_t *out) {
#ifdef RSIMD_NO_F64_SIMD
  return rsimd_ew_intdiv_i64_scalar(x, y, n, flags, mod, check, out);
#else
  int64_t x0 = 0, y0 = 0;
  rsimd_divmagic_i64 gs;
  const rsimd_divmagic_i64 *g;
  x = rsimd_intdiv_i64_scalar_arg(x, &flags, 0, check, &x0);
  y = rsimd_intdiv_i64_scalar_arg(y, &flags, 1, check, &y0);
  g = rsimd_intdiv_i64_magic(y, flags, check, &gs);
  if (flags & RSIMD_EW_SCALAR(1)) {
    if (y0 == 0 || y0 < -(INT64_C(1) << 51) || y0 >= (INT64_C(1) << 51)) {
      return rsimd_ew_intdiv_i64_scalar(x, y, n, flags, mod, check, out);
    }
  }
  if (flags == 0) {
    if (mod) {
      return check ? rsimd_ew_intdiv_i64_(x, y, n, 0, 1, 1, NULL, out)
                   : rsimd_ew_intdiv_i64_(x, y, n, 0, 1, 0, NULL, out);
    }
    return check ? rsimd_ew_intdiv_i64_(x, y, n, 0, 0, 1, NULL, out)
                 : rsimd_ew_intdiv_i64_(x, y, n, 0, 0, 0, NULL, out);
  }
  if (flags == RSIMD_EW_SCALAR(1)) {
    return mod ? rsimd_ew_intdiv_i64_(x, y, n, RSIMD_EW_SCALAR(1), 1, check, g, out)
               : rsimd_ew_intdiv_i64_(x, y, n, RSIMD_EW_SCALAR(1), 0, check, g, out);
  }
  if (flags == RSIMD_EW_I32(1)) {
    return rsimd_ew_intdiv_i64_(x, y, n, RSIMD_EW_I32(1), mod, check, g, out);
  }
  return rsimd_ew_intdiv_i64_(x, y, n, flags, mod, check, g, out);
#endif
}

/* Checked add or sub of a and b into r, NA where an operand is NA (with
   check) or the result is out of range; overflow outside the NA lanes
   sets the status. */
#define RSIMD_I64_CHECKED(r, a, b, opfn, ovffn)                                            \
  do {                                                                                     \
    rsimd_mi64 na_ = check ? rsimd_vi64_na2(a, b) : rsimd_mi64_none();                     \
    rsimd_mi64 o_;                                                                         \
    r = opfn(a, b);                                                                        \
    o_ = ovffn(a, b, r);                                                                   \
    if (rsimd_mi64_any(rsimd_mi64_andnot(na_, o_))) st |= RSIMD_EW_OVERFLOW;               \
    r = rsimd_vi64_set_na(r, rsimd_mi64_or(o_, na_));                                      \
  } while (0)
/* Wrapping op, NA where an operand is NA (with check). */
#define RSIMD_I64_WRAP(r, a, b, opfn)                                                      \
  do {                                                                                     \
    r = opfn(a, b);                                                                        \
    if (check) r = rsimd_vi64_set_na(r, rsimd_vi64_na2(a, b));                             \
  } while (0)

int RSIMD_KERNEL(ew2_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          int64_t *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew2_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                          int64_t *out, const rsimd_opts *o) {
  const void *z = NULL;
  const int check = rsimd_ew_i64_check(op, o);
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, check), bc1 = rsimd_i64_bcast(y, flags, 1, check),
                   bc2 = rsimd_vi64_zero();
  int st = 0;
  (void) z;
  switch (op) {
  case RSIMD_EW_ADD: RSIMD_I64_LOOP(2, RSIMD_I64_CHECKED(r, a, b, rsimd_vi64_add, rsimd_vi64_add_ovf)); break;
  case RSIMD_EW_SUB: RSIMD_I64_LOOP(2, RSIMD_I64_CHECKED(r, a, b, rsimd_vi64_sub, rsimd_vi64_sub_ovf)); break;
#if RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)
  case RSIMD_EW_MUL: RSIMD_I64_LOOP(2, r = rsimd_i64_mul_checked(a, b, check, &st)); break;
#else
  case RSIMD_EW_MUL:
    if (flags == 0) st = RSIMD_KERNEL(mul_i64_)(x, y, n, 0, out, check);
    else st = RSIMD_KERNEL(mul_i64_)(x, y, n, flags, out, check);
    break;
#endif
  case RSIMD_EW_ADD_WRAP: RSIMD_I64_LOOP(2, RSIMD_I64_WRAP(r, a, b, rsimd_vi64_add)); break;
  case RSIMD_EW_SUB_WRAP: RSIMD_I64_LOOP(2, RSIMD_I64_WRAP(r, a, b, rsimd_vi64_sub)); break;
  case RSIMD_EW_MUL_WRAP: RSIMD_I64_LOOP(2, RSIMD_I64_WRAP(r, a, b, rsimd_vi64_mul)); break;
  /* NA is INT64_MIN, so min propagates it and max ignores it by itself. */
  case RSIMD_EW_PMIN: RSIMD_I64_LOOP(2, r = rsimd_vi64_min(a, b)); break;
  case RSIMD_EW_PMAX:
    RSIMD_I64_LOOP(2, r = rsimd_vi64_set_na(rsimd_vi64_max(a, b), rsimd_vi64_na2(a, b)));
    break;
  case RSIMD_EW_PMIN_NUM:
    RSIMD_I64_LOOP(2, {
      r = rsimd_vi64_blend(rsimd_vi64_min(a, b), b, rsimd_vi64_is_na(a));
      r = rsimd_vi64_blend(r, a, rsimd_vi64_is_na(b));
    });
    break;
  case RSIMD_EW_PMAX_NUM: RSIMD_I64_LOOP(2, r = rsimd_vi64_max(a, b)); break;
  default: st = rsimd_ew_intdiv_i64(x, y, n, flags, op == RSIMD_EW_MOD, check, out); break;
  }
  return st;
}

int RSIMD_KERNEL(ew3_i64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, int64_t *out, const rsimd_opts *o);
int RSIMD_KERNEL(ew3_i64)(int op, const void *x, const void *y, const void *z, R_xlen_t n,
                          int flags, int64_t *out, const rsimd_opts *o) {
  const int check = op == RSIMD_EW_CLAMP ? 1 : o->na_check;
  int st = 0;
  R_xlen_t i;
  if (op == RSIMD_EW_CLAMP) {
    const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, 1), bc1 = rsimd_i64_bcast(y, flags, 1, 1),
                     bc2 = rsimd_i64_bcast(z, flags, 2, 1);
    RSIMD_I64_LOOP(3, {
      rsimd_mi64 bad = rsimd_mi64_andnot(rsimd_vi64_na2(b, c), rsimd_vi64_cmp_gt(b, c));
      if (rsimd_mi64_any(bad)) st |= RSIMD_EW_LO_GT_HI;
      r = rsimd_vi64_set_na(rsimd_vi64_max(a, b), rsimd_vi64_na2(a, b));
      r = rsimd_vi64_min(r, c);
    });
    return st;
  }
  /* mul_add and add_mul: scalar. */
  for (i = 0; i < n; i++) {
    out[i] = rsimd_ew3_i64_1(op, rsimd_i64_get(x, flags, 0, i, check),
                             rsimd_i64_get(y, flags, 1, i, check),
                             rsimd_i64_get(z, flags, 2, i, check), check, &st);
  }
  return st;
}

#undef RSIMD_I64_CHECKED
#undef RSIMD_I64_WRAP

/* ---- Comparisons and predicates ------------------------------------------------- */

RSIMD_ALWAYS_INLINE rsimd_mi64 rsimd_cmp_vi64(const int op, rsimd_vi64 a, rsimd_vi64 b) {
  switch (op) {
  case RSIMD_CMP_EQ: return rsimd_vi64_cmp_eq(a, b);
  case RSIMD_CMP_NE: return rsimd_mi64_not(rsimd_vi64_cmp_eq(a, b));
  case RSIMD_CMP_LT: return rsimd_vi64_cmp_lt(a, b);
  case RSIMD_CMP_LE: return rsimd_mi64_not(rsimd_vi64_cmp_gt(a, b));
  case RSIMD_CMP_GT: return rsimd_vi64_cmp_gt(a, b);
  default: return rsimd_mi64_not(rsimd_vi64_cmp_lt(a, b));
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(cmp_i64_)(const int op, const void *x, const void *y,
                                                R_xlen_t n, int flags, int *out) {
  const void *z = NULL;
  const int check = 1;
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, 1), bc1 = rsimd_i64_bcast(y, flags, 1, 1),
                   bc2 = rsimd_vi64_zero();
  (void) z;
  RSIMD_I64_LOOP32(2, r = rsimd_lgl_vi64m(rsimd_cmp_vi64(op, a, b), rsimd_vi64_na2(a, b)));
}

void RSIMD_KERNEL(cmp_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out);
void RSIMD_KERNEL(cmp_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                           int *out) {
  switch (op) {
  case RSIMD_CMP_EQ: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_EQ, x, y, n, flags, out); break;
  case RSIMD_CMP_NE: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_NE, x, y, n, flags, out); break;
  case RSIMD_CMP_LT: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_LT, x, y, n, flags, out); break;
  case RSIMD_CMP_LE: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_LE, x, y, n, flags, out); break;
  case RSIMD_CMP_GT: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_GT, x, y, n, flags, out); break;
  default: RSIMD_KERNEL(cmp_i64_)(RSIMD_CMP_GE, x, y, n, flags, out); break;
  }
}

RSIMD_ALWAYS_INLINE rsimd_mi64 rsimd_pred_vi64(const int op, rsimd_vi64 v) {
  const rsimd_vi64 zero = rsimd_vi64_zero(), one = rsimd_vi64_set1(1);
  switch (op) {
  case RSIMD_PRED_NA: return rsimd_vi64_is_na(v);
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return rsimd_mi64_not(rsimd_vi64_is_na(v));
  case RSIMD_PRED_NEGATIVE:
    return rsimd_mi64_andnot(rsimd_vi64_is_na(v), rsimd_vi64_cmp_lt(v, zero));
  case RSIMD_PRED_NORMAL:
    return rsimd_mi64_not(rsimd_mi64_or(rsimd_vi64_is_na(v), rsimd_vi64_cmp_eq(v, zero)));
  /* NA (INT64_MIN) has the low bit of an even number. */
  case RSIMD_PRED_EVEN:
    return rsimd_mi64_andnot(rsimd_vi64_is_na(v),
                             rsimd_vi64_cmp_eq(rsimd_vi64_and(v, one), zero));
  case RSIMD_PRED_ODD: return rsimd_vi64_cmp_eq(rsimd_vi64_and(v, one), one);
  case RSIMD_PRED_POW2:
    return rsimd_mi64_and(rsimd_vi64_cmp_gt(v, zero),
                          rsimd_vi64_cmp_eq(rsimd_vi64_and(v, rsimd_vi64_sub(v, one)), zero));
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return rsimd_mi64_none();
  default: return rsimd_vi64_cmp_eq(v, zero);
  }
}

/* In the last vector the inactive lanes repeat element i, so they cannot
   change an any/all answer. */
RSIMD_ALWAYS_INLINE int RSIMD_KERNEL(pred_i64_)(const int op, const int64_t *x, R_xlen_t n,
                                                int mode, int32_t *out) {
  const rsimd_mi64 none = rsimd_mi64_none();
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_mi64 m = rsimd_pred_vi64(op, rsimd_vi64_loadu(x + i));
    if (mode == RSIMD_PRED_ELT) rsimd_vi64_storeu_i32(out + i, rsimd_lgl_vi64m(m, none));
    else if (mode == RSIMD_PRED_ANY) {
      if (rsimd_mi64_any(m)) return 1;
    } else if (!rsimd_mi64_all(m)) {
      return 0;
    }
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_mi64 m = rsimd_pred_vi64(op, rsimd_vi64_loadu_p(pg, x + i, x[i]));
    if (mode == RSIMD_PRED_ELT) rsimd_vi64_storeu_i32_p(pg, out + i, rsimd_lgl_vi64m(m, none));
    else if (mode == RSIMD_PRED_ANY) return rsimd_mi64_any(m);
    else return rsimd_mi64_all(m);
  }
  return mode == RSIMD_PRED_ALL;
}

int RSIMD_KERNEL(pred_i64)(int op, const int64_t *x, R_xlen_t n, int mode, int *out);
int RSIMD_KERNEL(pred_i64)(int op, const int64_t *x, R_xlen_t n, int mode, int *out) {
  switch (op) {
  case RSIMD_PRED_NA: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_NA, x, n, mode, out);
  case RSIMD_PRED_FINITE:
  case RSIMD_PRED_WHOLE: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_FINITE, x, n, mode, out);
  case RSIMD_PRED_NEGATIVE: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_NEGATIVE, x, n, mode, out);
  case RSIMD_PRED_NORMAL: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_NORMAL, x, n, mode, out);
  case RSIMD_PRED_EVEN: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_EVEN, x, n, mode, out);
  case RSIMD_PRED_ODD: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_ODD, x, n, mode, out);
  case RSIMD_PRED_POW2: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_POW2, x, n, mode, out);
  case RSIMD_PRED_NAN:
  case RSIMD_PRED_INFINITE:
  case RSIMD_PRED_SUBNORMAL: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_NAN, x, n, mode, out);
  default: return RSIMD_KERNEL(pred_i64_)(RSIMD_PRED_ZERO, x, n, mode, out);
  }
}

/* ---- Bitwise ops -------------------------------------------------------------- */

/* op of the int64 lanes a and b (count k), without NA handling. */
RSIMD_ALWAYS_INLINE rsimd_vi64 rsimd_bit_vi64(const int op, rsimd_vi64 a, rsimd_vi64 b, int k) {
  switch (op) {
  case RSIMD_BIT_AND: return rsimd_vi64_and(a, b);
  case RSIMD_BIT_OR: return rsimd_vi64_or(a, b);
  case RSIMD_BIT_XOR: return rsimd_vi64_xor(a, b);
  case RSIMD_BIT_NOT: return rsimd_vi64_xor(a, rsimd_vi64_set1(-1));
  case RSIMD_BIT_SHL: return rsimd_vi64_sll(a, k);
  case RSIMD_BIT_SHR: return rsimd_vi64_srl(a, k);
  case RSIMD_BIT_SAR: return rsimd_vi64_sra(a, k);
  case RSIMD_BIT_ROTL:
    return k == 0 ? a : rsimd_vi64_or(rsimd_vi64_sll(a, k), rsimd_vi64_srl(a, 64 - k));
  case RSIMD_BIT_ROTR:
    return k == 0 ? a : rsimd_vi64_or(rsimd_vi64_srl(a, k), rsimd_vi64_sll(a, 64 - k));
  case RSIMD_BIT_POPCNT: return rsimd_vi64_popcnt(a);
  case RSIMD_BIT_LZCNT: return rsimd_vi64_clz(a);
  default: return rsimd_vi64_ctz(a);
  }
}

RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(bit_i64_)(const int op, const void *x, const void *y,
                                                R_xlen_t n, int flags, int k, void *out,
                                                const int nacheck) {
  const void *z = NULL;
  const int check = 1;
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, 1), bc1 = rsimd_i64_bcast(y, flags, 1, 1),
                   bc2 = rsimd_vi64_zero();
  (void) z;
  if (rsimd_bit_counts_i64(op)) {
    RSIMD_I64_LOOP32(1, {
      r = rsimd_bit_vi64(op, a, b, k);
      if (nacheck) r = rsimd_vi64_blend(r, rsimd_vi64_set1(RSIMD_NA_I32), rsimd_vi64_is_na(a));
    });
  } else {
    RSIMD_I64_LOOP(op <= RSIMD_BIT_XOR ? 2 : 1, {
      r = rsimd_bit_vi64(op, a, b, k);
      if (nacheck) {
        r = rsimd_vi64_set_na(r, op <= RSIMD_BIT_XOR ? rsimd_vi64_na2(a, b) : rsimd_vi64_is_na(a));
      }
    });
  }
}

/* Calls KERNEL(op, ...) with op as a constant. */
#define RSIMD_BIT_I64_DISPATCH(KERNEL, ...)                                                \
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

void RSIMD_KERNEL(bit_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags, int k,
                           void *out, const rsimd_opts *o);
void RSIMD_KERNEL(bit_i64)(int op, const void *x, const void *y, R_xlen_t n, int flags, int k,
                           void *out, const rsimd_opts *o) {
  if (op > RSIMD_BIT_XOR) y = NULL;
  if (o->na_check) {
    RSIMD_BIT_I64_DISPATCH(RSIMD_KERNEL(bit_i64_), x, y, n, flags, k, out, 1)
  } else {
    RSIMD_BIT_I64_DISPATCH(RSIMD_KERNEL(bit_i64_), x, y, n, flags, k, out, 0)
  }
}
#undef RSIMD_BIT_I64_DISPATCH

void RSIMD_KERNEL(popcnt_sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o);
void RSIMD_KERNEL(popcnt_sum_i64)(const int64_t *x, R_xlen_t n, rsimd_reduce_result *r,
                                  const rsimd_opts *o) {
  const int check = o->na_check || o->na_rm;
  int64_t buf[RSIMD_MAX_LANES_64];
  rsimd_vi64 acc = rsimd_vi64_zero();
  ptrdiff_t i = 0, j;
  for (; i < n; i += RSIMD_LANES_64) {
    /* The fill 0 counts nothing and is not NA. */
    rsimd_vi64 v = i + RSIMD_LANES_64 <= n ? rsimd_vi64_loadu(x + i)
                                           : rsimd_vi64_loadu_p(rsimd_p64_while(i, n), x + i, 0);
    rsimd_vi64 c = rsimd_vi64_popcnt(v);
    if (check) {
      rsimd_mi64 m = rsimd_vi64_is_na(v);
      if (rsimd_mi64_any(m)) {
        r->saw_na = 1;
        if (!o->na_rm) break;
        c = rsimd_vi64_blend(c, rsimd_vi64_zero(), m);
      }
    }
    acc = rsimd_vi64_add(acc, c);
  }
  rsimd_vi64_storeu(buf, acc);
  for (j = 0; j < RSIMD_LANES_64; j++) r->i64 += buf[j];
}

/* ---- Conversions --------------------------------------------------------------- */

RSIMD_INLINE int RSIMD_KERNEL(convert_i64_)(int op, int mode, const void *x, R_xlen_t n, void *out) {
  const int flags = 0, check = 1;
  const void *y = NULL, *z = NULL;
  const rsimd_vi64 bc0 = rsimd_vi64_zero(), bc1 = bc0, bc2 = bc0;
  int st = 0;
  R_xlen_t i;
  (void) y;
  (void) z;
  switch (op) {
  case RSIMD_CVT_I64_I32: {
    const rsimd_vi64 hi = rsimd_vi64_set1(INT32_MAX), lo = rsimd_vi64_set1(-INT32_MAX),
                     na = rsimd_vi64_set1(RSIMD_NA_I32);
    switch (mode) {
    case RSIMD_CVT_CHECKED:
      RSIMD_I64_LOOP32(1, {
        rsimd_mi64 m = rsimd_vi64_is_na(a);
        rsimd_mi64 oor = rsimd_mi64_or(rsimd_vi64_cmp_gt(a, hi), rsimd_vi64_cmp_lt(a, lo));
        if (rsimd_mi64_any(rsimd_mi64_andnot(m, oor))) st |= RSIMD_CVT_WARN_I32_OVF;
        r = rsimd_vi64_blend(a, na, rsimd_mi64_or(m, oor));
      });
      break;
    case RSIMD_CVT_SATURATING:
      RSIMD_I64_LOOP32(1, {
        r = rsimd_vi64_min(rsimd_vi64_max(a, lo), hi);
        r = rsimd_vi64_blend(r, na, rsimd_vi64_is_na(a));
      });
      break;
    default:
      /* The store keeps the low 32 bits; 0x80000000 there is NA anyway. */
      RSIMD_I64_LOOP32(1, r = rsimd_vi64_blend(a, na, rsimd_vi64_is_na(a)));
      break;
    }
    return st;
  }
  case RSIMD_CVT_I64_LGL:
    RSIMD_I64_LOOP32(1, {
      r = rsimd_vi64_blend(rsimd_vi64_set1(1), rsimd_vi64_zero(),
                           rsimd_vi64_cmp_eq(a, rsimd_vi64_zero()));
      r = rsimd_vi64_blend(r, rsimd_vi64_set1(RSIMD_NA_I32), rsimd_vi64_is_na(a));
    });
    return 0;
  case RSIMD_CVT_I32_I64: {
    const int32_t *xi = (const int32_t *) x;
    int64_t *po = (int64_t *) out;
    ptrdiff_t j = 0;
    for (; j + RSIMD_LANES_64 <= n; j += RSIMD_LANES_64) {
      rsimd_vi64_storeu(po + j, rsimd_i64_from_i32(rsimd_vi64_loadu_i32(xi + j), 1));
    }
    if (j < n) {
      rsimd_p64 pg = rsimd_p64_while(j, n);
      rsimd_vi64_storeu_p(pg, po + j, rsimd_i64_from_i32(rsimd_vi64_loadu_i32_p(pg, xi + j, 0), 1));
    }
    return 0;
  }
#ifndef RSIMD_NO_F64_SIMD
  case RSIMD_CVT_I64_F64: {
    const rsimd_vi64 big = rsimd_vi64_set1(INT64_C(9007199254740991));
    const rsimd_vf64 na = rsimd_vf64_set1(rsimd_na_real());
    double *po = (double *) out;
    rsimd_mi64 prec = rsimd_mi64_none();
    ptrdiff_t j = 0;
    /* |x| >= 2^53, outside the NA lanes, loses precision. */
#define RSIMD_I64_F64_(v, res)                                                             \
  do {                                                                                     \
    rsimd_mi64 m_ = rsimd_vi64_is_na(v);                                                   \
    rsimd_mi64 p_ = rsimd_mi64_or(rsimd_vi64_cmp_gt(v, big),                               \
                                  rsimd_vi64_cmp_lt(v, rsimd_vi64_sub(rsimd_vi64_zero(), big))); \
    prec = rsimd_mi64_or(prec, rsimd_mi64_andnot(m_, p_));                                 \
    res = rsimd_vf64_blend(rsimd_vi64_to_vf64(v), na, rsimd_mi64_to_mf64(m_));             \
  } while (0)
    for (; j + RSIMD_LANES_64 <= n; j += RSIMD_LANES_64) {
      rsimd_vi64 v = rsimd_vi64_loadu((const int64_t *) x + j);
      rsimd_vf64 res;
      RSIMD_I64_F64_(v, res);
      rsimd_vf64_storeu(po + j, res);
    }
    if (j < n) {
      rsimd_p64 pg = rsimd_p64_while(j, n);
      rsimd_vi64 v = rsimd_vi64_loadu_p(pg, (const int64_t *) x + j, 0);
      rsimd_vf64 res;
      RSIMD_I64_F64_(v, res);
      rsimd_vf64_storeu_p(pg, po + j, res);
    }
#undef RSIMD_I64_F64_
    return rsimd_mi64_any(prec) ? RSIMD_CVT_WARN_PRECISION : 0;
  }
#endif
#ifdef RSIMD_HAVE_I64_CVTT
  case RSIMD_CVT_F64_I64:
    if (mode != RSIMD_CVT_TRUNCATING) {
      const rsimd_vf64 two63 = rsimd_vf64_set1(RSIMD_TWO63);
      const rsimd_vi64 na = rsimd_vi64_set1(RSIMD_NA_I64), lo = rsimd_vi64_set1(-INT64_MAX);
      const double *xd = (const double *) x;
      int64_t *po = (int64_t *) out;
      rsimd_mf64 bad = rsimd_mf64_none();
      ptrdiff_t j = 0;
      /* Out of range: |x| >= 2^63 (not NaN); the conversion saturates. */
#define RSIMD_F64_I64_(v, res)                                                             \
  do {                                                                                     \
    rsimd_mf64 nan_ = rsimd_vf64_is_nan(v);                                                \
    rsimd_mf64 oor_ = rsimd_mf64_andnot(rsimd_vf64_cmp_lt(rsimd_vf64_abs(v), two63),       \
                                        rsimd_mf64_not(nan_));                             \
    res = rsimd_i64_cvtt_f64(v);                                                           \
    if (mode == RSIMD_CVT_CHECKED) {                                                       \
      bad = rsimd_mf64_or(bad, oor_);                                                      \
      res = rsimd_vi64_blend(res, na, rsimd_mf64_to_mi64(rsimd_mf64_or(nan_, oor_)));      \
    } else {                                                                               \
      res = rsimd_vi64_blend(rsimd_vi64_max(res, lo), na, rsimd_mf64_to_mi64(nan_));       \
    }                                                                                      \
  } while (0)
      for (; j + RSIMD_LANES_64 <= n; j += RSIMD_LANES_64) {
        rsimd_vi64 res;
        RSIMD_F64_I64_(rsimd_vf64_loadu(xd + j), res);
        rsimd_vi64_storeu(po + j, res);
      }
      if (j < n) {
        rsimd_p64 pg = rsimd_p64_while(j, n);
        rsimd_vi64 res;
        RSIMD_F64_I64_(rsimd_vf64_loadu_p(pg, xd + j, 0.0), res);
        rsimd_vi64_storeu_p(pg, po + j, res);
      }
#undef RSIMD_F64_I64_
      return rsimd_mf64_any(bad) ? RSIMD_CVT_WARN_I64 : 0;
    }
    break;
#endif
  default: break;
  }
  /* The rest: scalar. */
  for (i = 0; i < n; i++) st |= rsimd_cvt_i64_1(op, mode, x, i, out);
  (void) bc0;
  (void) bc1;
  (void) bc2;
  (void) flags;
  (void) check;
  return st;
}

/* ---- Hamming distance ------------------------------------------------------ */

/* Lane counters folded into r every RSIMD_HAMMING_BLOCK vectors; the
   last partial vector in scalar code. */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(hamming_i64_)(const void *x, const void *y, R_xlen_t n,
                                                    int flags, rsimd_reduce_result *r,
                                                    int na_rm) {
  const rsimd_vi64 bc0 = rsimd_i64_bcast(x, flags, 0, 1), bc1 = rsimd_i64_bcast(y, flags, 1, 1);
  ptrdiff_t i = 0;
  while (i + RSIMD_LANES_64 <= n) {
    const ptrdiff_t end = (n - i) / RSIMD_LANES_64 > RSIMD_HAMMING_BLOCK
                              ? i + RSIMD_HAMMING_BLOCK * RSIMD_LANES_64
                              : n;
    rsimd_vi64 d = rsimd_vi64_zero(), m = rsimd_vi64_zero();
    int64_t missing;
    for (; i + RSIMD_LANES_64 <= end; i += RSIMD_LANES_64) {
      rsimd_vi64 a = rsimd_i64_ld(x, flags, 0, bc0, i, 1), b = rsimd_i64_ld(y, flags, 1, bc1, i, 1);
      rsimd_mi64 na = rsimd_vi64_na2(a, b);
      d = rsimd_vi64_inc(d, rsimd_mi64_not(rsimd_mi64_or(rsimd_vi64_cmp_eq(a, b), na)));
      m = rsimd_vi64_inc(m, na);
    }
    missing = rsimd_vi64_reduce_add(m);
    r->i64 += rsimd_vi64_reduce_add(d);
    if (missing > 0) {
      r->saw_na = 1;
      if (!na_rm) return;
    }
  }
  rsimd_ham_i64_from(x, y, i, n, flags, r, na_rm);
}

void RSIMD_KERNEL(hamming_i64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o);
void RSIMD_KERNEL(hamming_i64)(const void *x, const void *y, R_xlen_t n, int flags,
                               rsimd_reduce_result *r, const rsimd_opts *o) {
  if (flags == 0) RSIMD_KERNEL(hamming_i64_)(x, y, n, 0, r, o->na_rm);
  else RSIMD_KERNEL(hamming_i64_)(x, y, n, flags, r, o->na_rm);
}

/* Popcounts of x ^ y summed in 64-bit lanes (which cannot overflow). */
RSIMD_ALWAYS_INLINE void RSIMD_KERNEL(hamming_bits_i64_)(const int64_t *x, const int64_t *y,
                                                         R_xlen_t n, int flags,
                                                         rsimd_reduce_result *r) {
  const int sx = (flags & RSIMD_EW_SCALAR(0)) != 0, sy = (flags & RSIMD_EW_SCALAR(1)) != 0;
  const rsimd_vi64 vx = rsimd_vi64_set1(x[0]), vy = rsimd_vi64_set1(y[0]);
  rsimd_vi64 acc = rsimd_vi64_zero();
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vi64 a = sx ? vx : rsimd_vi64_loadu(x + i), b = sy ? vy : rsimd_vi64_loadu(y + i);
    acc = rsimd_vi64_add(acc, rsimd_vi64_popcnt(rsimd_vi64_xor(a, b)));
  }
  r->i64 += rsimd_vi64_reduce_add(acc);
  rsimd_ham_bits_i64_from(x, y, i, n, flags, r);
}

void RSIMD_KERNEL(hamming_bits_i64)(const int64_t *x, const int64_t *y, R_xlen_t n, int flags,
                                    rsimd_reduce_result *r);
void RSIMD_KERNEL(hamming_bits_i64)(const int64_t *x, const int64_t *y, R_xlen_t n, int flags,
                                    rsimd_reduce_result *r) {
  if (flags == 0) RSIMD_KERNEL(hamming_bits_i64_)(x, y, n, 0, r);
  else RSIMD_KERNEL(hamming_bits_i64_)(x, y, n, flags, r);
}

#undef RSIMD_I64_LOOP
#undef RSIMD_I64_LOOP32
#undef RSIMD_I64_LOOP_T
#undef RSIMD_I64_LOOP_

#endif /* vector tiers */
