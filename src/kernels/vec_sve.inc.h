/* Vector layer of the sve and sve2 tiers, written with the Arm C Language
   Extensions (SIMDe does not map x86 intrinsics to SVE). The vector length is
   only known at run time: RSIMD_VLA is defined, RSIMD_LANES_64/32 call
   svcntd()/svcntw(), and there are no RSIMD_WIDTH_* constants. Masks and
   predicates are both svbool_t. Included by kernels/common.inc.h only. */

#include <arm_sve.h>

#define RSIMD_VLA 1
#define RSIMD_LANES_64 ((ptrdiff_t) svcntd())
#define RSIMD_LANES_32 ((ptrdiff_t) svcntw())

typedef svfloat64_t rsimd_vf64;
typedef svint32_t rsimd_vi32;
typedef svint64_t rsimd_vi64;
typedef svbool_t rsimd_mf64;
typedef svbool_t rsimd_mi32;
typedef svbool_t rsimd_mi64;
typedef svbool_t rsimd_p64;
typedef svbool_t rsimd_p32;

#define RSIMD_PT64 svptrue_b64()
#define RSIMD_PT32 svptrue_b32()

/* Predicates */
RSIMD_INLINE rsimd_p64 rsimd_p64_while(ptrdiff_t i, ptrdiff_t n) {
  return svwhilelt_b64((int64_t) i, (int64_t) n);
}
RSIMD_INLINE rsimd_p64 rsimd_p64_true(void) { return RSIMD_PT64; }
RSIMD_INLINE int rsimd_p64_count(rsimd_p64 pg) { return (int) svcntp_b64(RSIMD_PT64, pg); }
RSIMD_INLINE rsimd_p32 rsimd_p32_while(ptrdiff_t i, ptrdiff_t n) {
  return svwhilelt_b32((int64_t) i, (int64_t) n);
}
RSIMD_INLINE rsimd_p32 rsimd_p32_true(void) { return RSIMD_PT32; }
RSIMD_INLINE int rsimd_p32_count(rsimd_p32 pg) { return (int) svcntp_b32(RSIMD_PT32, pg); }

/* Masks */
#define RSIMD_SVE_MASK(m, pt, b)                                                  \
  RSIMD_INLINE rsimd_##m rsimd_##m##_and(rsimd_##m a, rsimd_##m c) {              \
    return svand_b_z(pt, a, c);                                                  \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_or(rsimd_##m a, rsimd_##m c) {               \
    return svorr_b_z(pt, a, c);                                                  \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##m##_not(rsimd_##m a) { return svnot_b_z(pt, a); } \
  RSIMD_INLINE rsimd_##m rsimd_##m##_andnot(rsimd_##m a, rsimd_##m c) {           \
    return svbic_b_z(pt, c, a);                                                  \
  }                                                                              \
  RSIMD_INLINE int rsimd_##m##_any(rsimd_##m a) { return svptest_any(pt, a); }    \
  RSIMD_INLINE int rsimd_##m##_all(rsimd_##m a) {                                 \
    return !svptest_any(pt, svnot_b_z(pt, a));                                   \
  }                                                                              \
  RSIMD_INLINE int rsimd_##m##_count(rsimd_##m a) { return (int) svcntp_##b(pt, a); }
RSIMD_SVE_MASK(mf64, RSIMD_PT64, b64)
RSIMD_SVE_MASK(mi32, RSIMD_PT32, b32)
RSIMD_SVE_MASK(mi64, RSIMD_PT64, b64)
#undef RSIMD_SVE_MASK

RSIMD_INLINE rsimd_mf64 rsimd_mi64_to_mf64(rsimd_mi64 m) { return m; }
RSIMD_INLINE rsimd_mi64 rsimd_mf64_to_mi64(rsimd_mf64 m) { return m; }

/* Doubles */
RSIMD_INLINE svuint64_t rsimd_sve_bits(rsimd_vf64 a) { return svreinterpret_u64_f64(a); }
RSIMD_INLINE rsimd_vf64 rsimd_sve_f64(svuint64_t a) { return svreinterpret_f64_u64(a); }

RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu(const double *p) { return svld1_f64(RSIMD_PT64, p); }
RSIMD_INLINE void rsimd_vf64_storeu(double *p, rsimd_vf64 v) { svst1_f64(RSIMD_PT64, p, v); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_p(rsimd_p64 pg, const double *p, double fill) {
  return svsel_f64(pg, svld1_f64(pg, p), svdup_n_f64(fill));
}
RSIMD_INLINE void rsimd_vf64_storeu_p(rsimd_p64 pg, double *p, rsimd_vf64 v) {
  svst1_f64(pg, p, v);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_set1(double x) { return svdup_n_f64(x); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_zero(void) { return svdup_n_f64(0.0); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_add(rsimd_vf64 a, rsimd_vf64 b) { return svadd_f64_x(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sub(rsimd_vf64 a, rsimd_vf64 b) { return svsub_f64_x(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_mul(rsimd_vf64 a, rsimd_vf64 b) { return svmul_f64_x(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_div(rsimd_vf64 a, rsimd_vf64 b) { return svdiv_f64_x(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_fma(rsimd_vf64 a, rsimd_vf64 b, rsimd_vf64 c) {
  return svmla_f64_x(RSIMD_PT64, c, a, b);
}
/* x86 semantics (b when either is NaN), not SVE's NaN-propagating FMIN. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_min(rsimd_vf64 a, rsimd_vf64 b) {
  return svsel_f64(svcmplt_f64(RSIMD_PT64, a, b), a, b);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_max(rsimd_vf64 a, rsimd_vf64 b) {
  return svsel_f64(svcmpgt_f64(RSIMD_PT64, a, b), a, b);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_abs(rsimd_vf64 a) { return svabs_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_neg(rsimd_vf64 a) { return svneg_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_sqrt(rsimd_vf64 a) { return svsqrt_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_and(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_sve_f64(svand_u64_x(RSIMD_PT64, rsimd_sve_bits(a), rsimd_sve_bits(b)));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_or(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_sve_f64(svorr_u64_x(RSIMD_PT64, rsimd_sve_bits(a), rsimd_sve_bits(b)));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_xor(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_sve_f64(sveor_u64_x(RSIMD_PT64, rsimd_sve_bits(a), rsimd_sve_bits(b)));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_andnot(rsimd_vf64 a, rsimd_vf64 b) {
  return rsimd_sve_f64(svbic_u64_x(RSIMD_PT64, rsimd_sve_bits(b), rsimd_sve_bits(a)));
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_eq(rsimd_vf64 a, rsimd_vf64 b) { return svcmpeq_f64(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ne(rsimd_vf64 a, rsimd_vf64 b) {
  return svnot_b_z(RSIMD_PT64, svcmpeq_f64(RSIMD_PT64, a, b));
}
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_lt(rsimd_vf64 a, rsimd_vf64 b) { return svcmplt_f64(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_le(rsimd_vf64 a, rsimd_vf64 b) { return svcmple_f64(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_gt(rsimd_vf64 a, rsimd_vf64 b) { return svcmpgt_f64(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_cmp_ge(rsimd_vf64 a, rsimd_vf64 b) { return svcmpge_f64(RSIMD_PT64, a, b); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_nan(rsimd_vf64 a) { return svcmpuo_f64(RSIMD_PT64, a, a); }
RSIMD_INLINE rsimd_mf64 rsimd_vf64_is_na(rsimd_vf64 a) {
  svuint64_t lo = svand_n_u64_x(RSIMD_PT64, rsimd_sve_bits(a), 0xFFFFFFFFu);
  return svand_b_z(RSIMD_PT64, svcmpeq_n_u64(RSIMD_PT64, lo, RSIMD_NA_LOW_WORD),
                   svcmpuo_f64(RSIMD_PT64, a, a));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_blend(rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m) {
  return svsel_f64(m, b, a);
}
/* Deinterleave: lanes 0, 2, 4 ... (uzp_even) or 1, 3, 5 ... (uzp_odd) of
   the concatenation a:b, so the real and the imaginary parts of the
   complex numbers held in a and then b. */
RSIMD_INLINE rsimd_vf64 rsimd_vf64_uzp_even(rsimd_vf64 a, rsimd_vf64 b) { return svuzp1_f64(a, b); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_uzp_odd(rsimd_vf64 a, rsimd_vf64 b) { return svuzp2_f64(a, b); }
RSIMD_INLINE double rsimd_vf64_reduce_add(rsimd_vf64 a) { return svaddv_f64(RSIMD_PT64, a); }
RSIMD_INLINE double rsimd_vf64_reduce_min(rsimd_vf64 a) { return svminv_f64(RSIMD_PT64, a); }
RSIMD_INLINE double rsimd_vf64_reduce_max(rsimd_vf64 a) { return svmaxv_f64(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vi64 rsimd_vf64_as_vi64(rsimd_vf64 a) { return svreinterpret_s64_f64(a); }
RSIMD_INLINE rsimd_vf64 rsimd_vi64_as_vf64(rsimd_vi64 a) { return svreinterpret_f64_s64(a); }

/* Integers: SVE integer arithmetic wraps. Arguments: vector suffix, element
   type, lane suffix, all-true predicate, ACLE type suffix. */
#define RSIMD_SVE_INT(v, m, p, elt, pt, s)                                        \
  RSIMD_INLINE rsimd_##v rsimd_##v##_loadu(const elt *ptr) { return svld1_##s(pt, ptr); } \
  RSIMD_INLINE void rsimd_##v##_storeu(elt *ptr, rsimd_##v x) { svst1_##s(pt, ptr, x); } \
  RSIMD_INLINE rsimd_##v rsimd_##v##_loadu_p(rsimd_##p pg, const elt *ptr, elt fill) { \
    return svsel_##s(pg, svld1_##s(pg, ptr), svdup_n_##s(fill));                 \
  }                                                                              \
  RSIMD_INLINE void rsimd_##v##_storeu_p(rsimd_##p pg, elt *ptr, rsimd_##v x) {   \
    svst1_##s(pg, ptr, x);                                                       \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_set1(elt x) { return svdup_n_##s(x); }       \
  RSIMD_INLINE rsimd_##v rsimd_##v##_zero(void) { return svdup_n_##s(0); }        \
  RSIMD_INLINE rsimd_##v rsimd_##v##_add(rsimd_##v a, rsimd_##v b) {              \
    return svadd_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_sub(rsimd_##v a, rsimd_##v b) {              \
    return svsub_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_mul(rsimd_##v a, rsimd_##v b) {              \
    return svmul_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_min(rsimd_##v a, rsimd_##v b) {              \
    return svmin_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_max(rsimd_##v a, rsimd_##v b) {              \
    return svmax_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_and(rsimd_##v a, rsimd_##v b) {              \
    return svand_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_or(rsimd_##v a, rsimd_##v b) {               \
    return svorr_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_xor(rsimd_##v a, rsimd_##v b) {              \
    return sveor_##s##_x(pt, a, b);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_andnot(rsimd_##v a, rsimd_##v b) {           \
    return svbic_##s##_x(pt, b, a);                                              \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##v##_cmp_eq(rsimd_##v a, rsimd_##v b) {           \
    return svcmpeq_##s(pt, a, b);                                                \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##v##_cmp_gt(rsimd_##v a, rsimd_##v b) {           \
    return svcmpgt_##s(pt, a, b);                                                \
  }                                                                              \
  RSIMD_INLINE rsimd_##m rsimd_##v##_cmp_lt(rsimd_##v a, rsimd_##v b) {           \
    return svcmplt_##s(pt, a, b);                                                \
  }                                                                              \
  RSIMD_INLINE rsimd_##v rsimd_##v##_blend(rsimd_##v a, rsimd_##v b, rsimd_##m k) { \
    return svsel_##s(k, b, a);                                                   \
  }                                                                              \
  RSIMD_INLINE elt rsimd_##v##_reduce_min(rsimd_##v a) { return svminv_##s(pt, a); } \
  RSIMD_INLINE elt rsimd_##v##_reduce_max(rsimd_##v a) { return svmaxv_##s(pt, a); }
RSIMD_SVE_INT(vi32, mi32, p32, int32_t, RSIMD_PT32, s32)
RSIMD_SVE_INT(vi64, mi64, p64, int64_t, RSIMD_PT64, s64)
#undef RSIMD_SVE_INT

RSIMD_INLINE rsimd_mi32 rsimd_vi32_is_na(rsimd_vi32 a) {
  return svcmpeq_n_s32(RSIMD_PT32, a, INT32_MIN);
}
RSIMD_INLINE rsimd_mi64 rsimd_vi64_is_na(rsimd_vi64 a) {
  return svcmpeq_n_s64(RSIMD_PT64, a, INT64_MIN);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_inc(rsimd_vi32 acc, rsimd_mi32 m) {
  return svadd_n_s32_m(m, acc, 1);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_inc(rsimd_vi64 acc, rsimd_mi64 m) {
  return svadd_n_s64_m(m, acc, 1);
}
/* SADDV widens 32-bit lanes to a 64-bit sum, so this is exact. */
RSIMD_INLINE int64_t rsimd_vi32_reduce_add(rsimd_vi32 a) { return svaddv_s32(RSIMD_PT32, a); }
RSIMD_INLINE int64_t rsimd_vi64_reduce_add(rsimd_vi64 a) { return svaddv_s64(RSIMD_PT64, a); }

/* Width conversions and extras. The 64-bit-lane conversions move
   svcntd() int32 elements. */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_mulhi(rsimd_vi32 a, rsimd_vi32 b) {
  return svmulh_s32_x(RSIMD_PT32, a, b);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sll(rsimd_vi32 a, int k) {
  return svlsl_n_s32_x(RSIMD_PT32, a, (uint32_t) k);
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_srl(rsimd_vi32 a, int k) {
  return svreinterpret_s32_u32(svlsr_n_u32_x(RSIMD_PT32, svreinterpret_u32_s32(a), (uint32_t) k));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_sra(rsimd_vi32 a, int k) {
  return svasr_n_s32_x(RSIMD_PT32, a, (uint32_t) k);
}
/* svcntw() bytes zero-extended to 32-bit lanes, and back (low 8 bits). */
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_u8(const uint8_t *p) {
  return svreinterpret_s32_u32(svld1ub_u32(RSIMD_PT32, p));
}
RSIMD_INLINE rsimd_vi32 rsimd_vi32_loadu_u8_p(rsimd_p32 pg, const uint8_t *p, uint8_t fill) {
  return svsel_s32(pg, svreinterpret_s32_u32(svld1ub_u32(pg, p)), svdup_n_s32(fill));
}
RSIMD_INLINE void rsimd_vi32_storeu_u8(uint8_t *p, rsimd_vi32 v) {
  svst1b_u32(RSIMD_PT32, p, svreinterpret_u32_s32(v));
}
RSIMD_INLINE void rsimd_vi32_storeu_u8_p(rsimd_p32 pg, uint8_t *p, rsimd_vi32 v) {
  svst1b_u32(pg, p, svreinterpret_u32_s32(v));
}
/* 64-bit lanes: shifts by a run-time count k in [0, 63], the sign mask,
   the exact products of the low 32 bits (unsigned, signed) and the
   conversion to double (rounded to nearest). */
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sll(rsimd_vi64 a, int k) {
  return svlsl_n_s64_x(RSIMD_PT64, a, (uint64_t) k);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_srl(rsimd_vi64 a, int k) {
  return svreinterpret_s64_u64(svlsr_n_u64_x(RSIMD_PT64, svreinterpret_u64_s64(a), (uint64_t) k));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sra(rsimd_vi64 a, int k) {
  return svasr_n_s64_x(RSIMD_PT64, a, (uint64_t) k);
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_sign(rsimd_vi64 a) { return svasr_n_s64_x(RSIMD_PT64, a, 63); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_mulu32(rsimd_vi64 a, rsimd_vi64 b) {
  svuint64_t ua = svextw_u64_x(RSIMD_PT64, svreinterpret_u64_s64(a));
  svuint64_t ub = svextw_u64_x(RSIMD_PT64, svreinterpret_u64_s64(b));
  return svreinterpret_s64_u64(svmul_u64_x(RSIMD_PT64, ua, ub));
}
RSIMD_INLINE rsimd_vi64 rsimd_vi64_mul32(rsimd_vi64 a, rsimd_vi64 b) {
  return svmul_s64_x(RSIMD_PT64, svextw_s64_x(RSIMD_PT64, a), svextw_s64_x(RSIMD_PT64, b));
}
RSIMD_INLINE rsimd_vf64 rsimd_vi64_to_vf64(rsimd_vi64 a) { return svcvt_f64_s64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32(const int32_t *p) { return svld1sw_s64(RSIMD_PT64, p); }
RSIMD_INLINE rsimd_vi64 rsimd_vi64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return svsel_s64(pg, svld1sw_s64(pg, p), svdup_n_s64(fill));
}
RSIMD_INLINE void rsimd_vi64_storeu_i32(int32_t *p, rsimd_vi64 v) { svst1w_s64(RSIMD_PT64, p, v); }
RSIMD_INLINE void rsimd_vi64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vi64 v) {
  svst1w_s64(pg, p, v);
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_floor(rsimd_vf64 a) { return svrintm_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_ceil(rsimd_vf64 a) { return svrintp_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_trunc(rsimd_vf64 a) { return svrintz_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_rint(rsimd_vf64 a) { return svrintn_f64_x(RSIMD_PT64, a); }
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32(const int32_t *p) {
  return svcvt_f64_s64_x(RSIMD_PT64, rsimd_vi64_loadu_i32(p));
}
RSIMD_INLINE rsimd_vf64 rsimd_vf64_loadu_i32_p(rsimd_p64 pg, const int32_t *p, int32_t fill) {
  return svcvt_f64_s64_x(RSIMD_PT64, rsimd_vi64_loadu_i32_p(pg, p, fill));
}
RSIMD_INLINE void rsimd_vf64_storeu_i32(int32_t *p, rsimd_vf64 v) {
  svst1w_s64(RSIMD_PT64, p, svcvt_s64_f64_x(RSIMD_PT64, v));
}
RSIMD_INLINE void rsimd_vf64_storeu_i32_p(rsimd_p64 pg, int32_t *p, rsimd_vf64 v) {
  svst1w_s64(pg, p, svcvt_s64_f64_x(RSIMD_PT64, v));
}
