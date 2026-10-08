/* Helpers shared by the kernels of predicates.inc.c, compare.inc.c and
 * bitwise.inc.c, whose results are logical vectors, int32 counts or bytes.
 * Included by those files after arith.inc.c (whose rsimd_ew_ld helpers
 * read mixed double and int32 operands); the include guard makes the
 * second inclusion empty.
 *
 * Vector tiers only:
 *   - rsimd_lgl_vi32(t, na) and rsimd_lgl_vi64(t, na): logical lanes, 1
 *     where t, 0 elsewhere and NA (INT32_MIN) where na. The 64-bit form is
 *     stored as int32 elements with rsimd_vi64_storeu_i32.
 *   - RSIMD_LANE32_LOOP: a loop over int32 or byte operands x and y (y may
 *     be NULL) with RSIMD_EW_SCALAR(k) broadcast flags.
 *   - RSIMD_LANE64_LOOP: a loop over double or int32 operands (flags
 *     RSIMD_EW_SCALAR(k), RSIMD_EW_I32(k)) writing int32 elements from
 *     64-bit lanes.
 * In the predicated last vector the inactive lanes repeat the first active
 * element of each operand, as in arith.inc.c.
 */

#ifndef RSIMD_KERNELS_LOGICAL_INC_H
#define RSIMD_KERNELS_LOGICAL_INC_H

#define RSIMD_LGL_SCALAR(k) ((flags & RSIMD_EW_SCALAR(k)) != 0)

/* Vectors between the folds of the Hamming kernels' lane counters into
   the result (and their checks for missing pairs): few enough for 32-bit
   lane counters. Tests set it smaller. */
#ifndef RSIMD_HAMMING_BLOCK
#define RSIMD_HAMMING_BLOCK ((ptrdiff_t) 1 << 12)
#endif

#if !RSIMD_TIER_IS(none)

RSIMD_INLINE rsimd_vi32 rsimd_lgl_vi32(rsimd_mi32 t, rsimd_mi32 na) {
  rsimd_vi32 r = rsimd_vi32_blend(rsimd_vi32_zero(), rsimd_vi32_set1(1), t);
  return rsimd_vi32_blend(r, rsimd_vi32_set1(RSIMD_NA_I32), na);
}

#ifndef RSIMD_NO_F64_SIMD
RSIMD_INLINE rsimd_vi64 rsimd_lgl_vi64(rsimd_mf64 t, rsimd_mf64 na) {
  rsimd_vi64 r = rsimd_vi64_blend(rsimd_vi64_zero(), rsimd_vi64_set1(1), rsimd_mf64_to_mi64(t));
  return rsimd_vi64_blend(r, rsimd_vi64_set1(RSIMD_NA_I32), rsimd_mf64_to_mi64(na));
}
#endif

/* Runs the statements in `...`, which set the rsimd_vi32 r from the operands a and b,
   over the chunk of n elements of x and y (pointers to elements of type
   T, read with LD / LD_P: rsimd_vi32_loadu(_p) for int32_t,
   rsimd_vi32_loadu_u8(_p) for uint8_t), and stores r at out as elements
   of type OT with ST / ST_P. A NULL y reads as zeros. */
#define RSIMD_LANE32_LOOP(T, LD, LD_P, OT, ST, ST_P, ...)                                   \
  do {                                                                                     \
    const T *xp_ = (const T *) x, *yp_ = (const T *) y;                                    \
    const rsimd_vi32 bx_ = rsimd_vi32_set1((int32_t) xp_[0]),                              \
                     by_ = yp_ != NULL ? rsimd_vi32_set1((int32_t) yp_[0]) : rsimd_vi32_zero(); \
    const int sx_ = RSIMD_LGL_SCALAR(0), sy_ = yp_ == NULL || RSIMD_LGL_SCALAR(1);         \
    ptrdiff_t i = 0;                                                                       \
    for (; i + RSIMD_LANES_32 <= n; i += RSIMD_LANES_32) {                                 \
      rsimd_vi32 a = sx_ ? bx_ : LD(xp_ + i), b = sy_ ? by_ : LD(yp_ + i), r;              \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      ST((OT *) out + i, r);                                                               \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p32 pg = rsimd_p32_while(i, n);                                                \
      rsimd_vi32 a = sx_ ? bx_ : LD_P(pg, xp_ + i, xp_[i]);                                \
      rsimd_vi32 b = sy_ ? by_ : LD_P(pg, yp_ + i, yp_[i]), r;                             \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      ST_P(pg, (OT *) out + i, r);                                                         \
    }                                                                                      \
  } while (0)

#ifndef RSIMD_NO_F64_SIMD
/* Runs the statements in `...`, which set the rsimd_vi64 r (logical lanes) from the
   rsimd_vf64 operands a and b, over the chunk, with operand flags fl (int32
   operands become doubles, NA becoming NA_real_). A NULL y reads as
   zeros. Needs bc0 and bc1 from rsimd_ew_bcast(). */
#define RSIMD_LANE64_LOOP(fl, ...)                                                         \
  do {                                                                                     \
    ptrdiff_t i = 0;                                                                       \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                                 \
      rsimd_vf64 a = rsimd_ew_ld(x, fl, 0, bc0, i, 1);                                     \
      rsimd_vf64 b = y == NULL ? bc1 : rsimd_ew_ld(y, fl, 1, bc1, i, 1);                   \
      rsimd_vi64 r;                                                                        \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      rsimd_vi64_storeu_i32((int32_t *) out + i, r);                                       \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                                \
      rsimd_vf64 a = rsimd_ew_ld_p(x, fl, 0, bc0, i, pg, 1);                               \
      rsimd_vf64 b = y == NULL ? bc1 : rsimd_ew_ld_p(y, fl, 1, bc1, i, pg, 1);             \
      rsimd_vi64 r;                                                                        \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      rsimd_vi64_storeu_i32_p(pg, (int32_t *) out + i, r);                                 \
    }                                                                                      \
  } while (0)
/* As RSIMD_LANE64_LOOP, reading the operands through rsimd_ew_dptr(). */
#define RSIMD_LANE64_LOOP_D(...)                                                           \
  do {                                                                                     \
    ptrdiff_t i = 0;                                                                       \
    const rsimd_vf64 bc2 = bc1;                                                            \
    RSIMD_EW_DPTRS(x, y, NULL);                                                            \
    (void) p1_;                                                                            \
    (void) p2_;                                                                            \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                                 \
      rsimd_vf64 a = RSIMD_EW_DLD(0), b = RSIMD_EW_DLD(1);                                 \
      rsimd_vi64 r;                                                                        \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      rsimd_vi64_storeu_i32((int32_t *) out + i, r);                                       \
      RSIMD_EW_DNEXT();                                                                    \
    }                                                                                      \
    if (i < n) {                                                                           \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                                \
      rsimd_vf64 a = RSIMD_EW_DLD_P(0), b = RSIMD_EW_DLD_P(1);                             \
      rsimd_vi64 r;                                                                        \
      (void) b;                                                                            \
      __VA_ARGS__;                                                                          \
      rsimd_vi64_storeu_i32_p(pg, (int32_t *) out + i, r);                                 \
    }                                                                                      \
  } while (0)
/* The loop with all-vector double operands as a constant case, the loop
   over double operands, vectors or scalars, without per-operand tests
   otherwise, and the general one for int32 vectors. */
#define RSIMD_LANE64_LOOP_FLAGS(...)                                                       \
  do {                                                                                     \
    if (flags == 0) {                                                                      \
      RSIMD_LANE64_LOOP(0, __VA_ARGS__);                                                      \
    } else if (RSIMD_EW_NO_I32_VECTOR(flags)) {                                            \
      RSIMD_LANE64_LOOP_D(__VA_ARGS__);                                                    \
    } else {                                                                               \
      RSIMD_LANE64_LOOP(flags, __VA_ARGS__);                                                  \
    }                                                                                      \
  } while (0)
#endif /* RSIMD_NO_F64_SIMD */

#endif /* vector tiers */

#endif /* RSIMD_KERNELS_LOGICAL_INC_H */
