/* Parts of the vector layer shared by the 128- and 256-bit SIMDe tiers
   (sse2, neon, avx2): lane-count predicates, predicated loads and stores
   through a stack buffer, and horizontal reductions through a store.
   Included at the end of vec_simde128.inc.h and vec_simde256.inc.h, after
   the types, RSIMD_WIDTH_* and the full-width loadu/storeu (including the
   loadu_i32/storeu_i32 conversions) are defined. */

/* A predicate is the number of active lanes, counted from lane 0. */
typedef int rsimd_p64;
typedef int rsimd_p32;

RSIMD_INLINE rsimd_p64 rsimd_p64_while(ptrdiff_t i, ptrdiff_t n) {
  ptrdiff_t k = n - i;
  return k >= RSIMD_WIDTH_I64 ? RSIMD_WIDTH_I64 : (k > 0 ? (int) k : 0);
}
RSIMD_INLINE rsimd_p64 rsimd_p64_true(void) { return RSIMD_WIDTH_I64; }
RSIMD_INLINE int rsimd_p64_count(rsimd_p64 pg) { return pg; }
RSIMD_INLINE rsimd_p32 rsimd_p32_while(ptrdiff_t i, ptrdiff_t n) {
  ptrdiff_t k = n - i;
  return k >= RSIMD_WIDTH_I32 ? RSIMD_WIDTH_I32 : (k > 0 ? (int) k : 0);
}
RSIMD_INLINE rsimd_p32 rsimd_p32_true(void) { return RSIMD_WIDTH_I32; }
RSIMD_INLINE int rsimd_p32_count(rsimd_p32 pg) { return pg; }

#define RSIMD_FIXED_PARTIAL(v, elt, p, w)                                        \
  RSIMD_INLINE rsimd_##v rsimd_##v##_loadu_p(rsimd_##p pg, const elt *ptr,      \
                                             elt fill) {                         \
    elt buf[w];                                                                  \
    int j;                                                                       \
    if (pg >= (w)) return rsimd_##v##_loadu(ptr);                                \
    for (j = 0; j < (w); j++) buf[j] = j < pg ? ptr[j] : fill;                   \
    return rsimd_##v##_loadu(buf);                                               \
  }                                                                              \
  RSIMD_INLINE void rsimd_##v##_storeu_p(rsimd_##p pg, elt *ptr, rsimd_##v x) {  \
    elt buf[w];                                                                  \
    int j;                                                                       \
    if (pg >= (w)) {                                                             \
      rsimd_##v##_storeu(ptr, x);                                                \
      return;                                                                    \
    }                                                                            \
    rsimd_##v##_storeu(buf, x);                                                  \
    for (j = 0; j < pg; j++) ptr[j] = buf[j];                                    \
  }

/* Lane-order folds with the layer's min/max semantics. */
#define RSIMD_FIXED_MINMAX(v, elt, w)                                            \
  RSIMD_INLINE elt rsimd_##v##_reduce_min(rsimd_##v x) {                         \
    elt buf[w], r;                                                               \
    int j;                                                                       \
    rsimd_##v##_storeu(buf, x);                                                  \
    r = buf[0];                                                                  \
    for (j = 1; j < (w); j++) r = r < buf[j] ? r : buf[j];                       \
    return r;                                                                    \
  }                                                                              \
  RSIMD_INLINE elt rsimd_##v##_reduce_max(rsimd_##v x) {                         \
    elt buf[w], r;                                                               \
    int j;                                                                       \
    rsimd_##v##_storeu(buf, x);                                                  \
    r = buf[0];                                                                  \
    for (j = 1; j < (w); j++) r = r > buf[j] ? r : buf[j];                       \
    return r;                                                                    \
  }

#ifndef RSIMD_NO_F64_SIMD
RSIMD_FIXED_PARTIAL(vf64, double, p64, RSIMD_WIDTH_F64)
RSIMD_FIXED_MINMAX(vf64, double, RSIMD_WIDTH_F64)
RSIMD_INLINE double rsimd_vf64_reduce_add(rsimd_vf64 x) {
  double buf[RSIMD_WIDTH_F64], r;
  int j;
  rsimd_vf64_storeu(buf, x);
  r = buf[0];
  for (j = 1; j < RSIMD_WIDTH_F64; j++) r += buf[j];
  return r;
}
#endif

RSIMD_FIXED_PARTIAL(vi32, int32_t, p32, RSIMD_WIDTH_I32)
RSIMD_FIXED_MINMAX(vi32, int32_t, RSIMD_WIDTH_I32)
RSIMD_INLINE int64_t rsimd_vi32_reduce_add(rsimd_vi32 x) {
  int32_t buf[RSIMD_WIDTH_I32];
  int64_t r = 0;
  int j;
  rsimd_vi32_storeu(buf, x);
  for (j = 0; j < RSIMD_WIDTH_I32; j++) r += buf[j];
  return r;
}

RSIMD_FIXED_PARTIAL(vi64, int64_t, p64, RSIMD_WIDTH_I64)
RSIMD_FIXED_MINMAX(vi64, int64_t, RSIMD_WIDTH_I64)
RSIMD_INLINE int64_t rsimd_vi64_reduce_add(rsimd_vi64 x) {
  int64_t buf[RSIMD_WIDTH_I64];
  uint64_t r = 0;
  int j;
  rsimd_vi64_storeu(buf, x);
  for (j = 0; j < RSIMD_WIDTH_I64; j++) r += (uint64_t) buf[j];
  return (int64_t) r;
}

/* 64-bit multiply (low 64 bits, wrapping) from 32-bit partial products,
   as no tier using this file has one: a_lo * b_lo + ((a_lo * b_hi +
   a_hi * b_lo) << 32). */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_mul(rsimd_vi64 a, rsimd_vi64 b) {
  rsimd_vi64 cross = rsimd_vi64_add(rsimd_vi64_mulu32(a, rsimd_vi64_srl(b, 32)),
                                    rsimd_vi64_mulu32(rsimd_vi64_srl(a, 32), b));
  return rsimd_vi64_add(rsimd_vi64_mulu32(a, b), rsimd_vi64_sll(cross, 32));
}

/* Arithmetic shift right by k in [0, 63]: the logical shift of the
   complement of a negative lane, complemented back. */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sra(rsimd_vi64 a, int k) {
  rsimd_vi64 s = rsimd_vi64_sign(a);
  return rsimd_vi64_xor(rsimd_vi64_srl(rsimd_vi64_xor(a, s), k), s);
}

#if !defined(RSIMD_NO_F64_SIMD) && !defined(RSIMD_HAVE_VI64_TO_VF64)
/* int64 lanes to double, rounded to nearest (even), with exact magic
   numbers: the high word h, biased to h + 2^31 and placed in the mantissa
   of 2^84, is 2^84 + (h + 2^31) * 2^32; the low word l in the mantissa of
   2^52 is 2^52 + l. Subtracting 2^84 + 2^63 + 2^52 from the first is exact
   (h * 2^32 - 2^52), and adding the second then rounds once. */
RSIMD_INLINE rsimd_vf64 rsimd_vi64_to_vf64(rsimd_vi64 a) {
  const rsimd_vi64 hi_bits = rsimd_vi64_set1(INT64_C(0x4530000080000000));
  const rsimd_vi64 lo_bits = rsimd_vi64_set1(INT64_C(0x4330000000000000));
  const rsimd_vi64 lo_mask = rsimd_vi64_set1(INT64_C(0xFFFFFFFF));
  rsimd_vf64 hi = rsimd_vi64_as_vf64(rsimd_vi64_xor(rsimd_vi64_srl(a, 32), hi_bits));
  rsimd_vf64 lo = rsimd_vi64_as_vf64(rsimd_vi64_or(rsimd_vi64_and(a, lo_mask), lo_bits));
  return rsimd_vf64_add(rsimd_vf64_sub(hi, rsimd_vf64_set1(0x1.00000801p84)), lo);
}
#endif

/* Predicated forms of the int32 <-> 64-bit-lane conversions, through a
   buffer of int32 elements. */
#define RSIMD_FIXED_I32_PARTIAL(v)                                               \
  RSIMD_INLINE rsimd_##v rsimd_##v##_loadu_i32_p(rsimd_p64 pg, const int32_t *ptr, \
                                                 int32_t fill) {                 \
    int32_t buf[RSIMD_WIDTH_I64];                                                \
    int j;                                                                       \
    if (pg >= RSIMD_WIDTH_I64) return rsimd_##v##_loadu_i32(ptr);                \
    for (j = 0; j < RSIMD_WIDTH_I64; j++) buf[j] = j < pg ? ptr[j] : fill;       \
    return rsimd_##v##_loadu_i32(buf);                                           \
  }                                                                              \
  RSIMD_INLINE void rsimd_##v##_storeu_i32_p(rsimd_p64 pg, int32_t *ptr, rsimd_##v x) { \
    int32_t buf[RSIMD_WIDTH_I64];                                                \
    int j;                                                                       \
    if (pg >= RSIMD_WIDTH_I64) {                                                 \
      rsimd_##v##_storeu_i32(ptr, x);                                            \
      return;                                                                    \
    }                                                                            \
    rsimd_##v##_storeu_i32(buf, x);                                              \
    for (j = 0; j < pg; j++) ptr[j] = buf[j];                                    \
  }

/* Predicated forms of the byte <-> 32-bit-lane conversions. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_u8_p(rsimd_p32 pg, const uint8_t *p, uint8_t fill) {
  uint8_t buf[RSIMD_WIDTH_I32];
  int j;
  if (pg >= RSIMD_WIDTH_I32) return rsimd_vi32_loadu_u8(p);
  for (j = 0; j < RSIMD_WIDTH_I32; j++) buf[j] = j < pg ? p[j] : fill;
  return rsimd_vi32_loadu_u8(buf);
}
RSIMD_INLINE void rsimd_vi32_storeu_u8_p(rsimd_p32 pg, uint8_t *p, rsimd_vi32 v) {
  uint8_t buf[RSIMD_WIDTH_I32];
  int j;
  if (pg >= RSIMD_WIDTH_I32) {
    rsimd_vi32_storeu_u8(p, v);
    return;
  }
  rsimd_vi32_storeu_u8(buf, v);
  for (j = 0; j < pg; j++) p[j] = buf[j];
}

RSIMD_FIXED_I32_PARTIAL(vi64)
#ifndef RSIMD_NO_F64_SIMD
RSIMD_FIXED_I32_PARTIAL(vf64)
#endif

#undef RSIMD_FIXED_PARTIAL
#undef RSIMD_FIXED_MINMAX
#undef RSIMD_FIXED_I32_PARTIAL
