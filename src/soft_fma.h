#ifndef RSIMD_SOFT_FMA_H
#define RSIMD_SOFT_FMA_H

/* A correctly rounded fused multiply-add for every scalar fma in the
 * package: rsimd_fma(x, y, z) is x * y + z rounded once.
 *
 * The C library's fma is used where it is exact, and compiles to one
 * instruction where the tier has FMA (avx2, avx512, arm64). The Windows C
 * library's fma rounds some results twice, so without hardware FMA (the
 * none and sse2 tiers, and the files compiled with R's flags) Windows
 * builds use rsimd_soft_fma below. Define RSIMD_SOFT_FMA=1 to use it
 * everywhere, for testing.
 *
 * rsimd_soft_fma is musl's fma (src/math/fma.c, MIT licence, see
 * inst/COPYRIGHTS), unchanged except for names, the bit casts and count of
 * leading zeros written portably, and the underflow exception not raised. */

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define RSIMD_SOFT_FMA_ZEROINFNAN (0x7ff - 0x3ff - 52 - 1)

struct rsimd_soft_fma_num {
  uint64_t m;
  int e;
  int sign;
};

static inline uint64_t rsimd_soft_fma_bits(double x) {
  uint64_t i;
  memcpy(&i, &x, sizeof i);
  return i;
}

static inline struct rsimd_soft_fma_num rsimd_soft_fma_normalize(double x) {
  uint64_t ix = rsimd_soft_fma_bits(x);
  int e = (int) (ix >> 52);
  int sign = e & 0x800;
  e &= 0x7ff;
  if (!e) {
    ix = rsimd_soft_fma_bits(x * 0x1p63);
    e = (int) (ix >> 52 & 0x7ff);
    e = e ? e - 63 : 0x800;
  }
  ix &= (1ull << 52) - 1;
  ix |= 1ull << 52;
  ix <<= 1;
  e -= 0x3ff + 52 + 1;
  struct rsimd_soft_fma_num n = {ix, e, sign};
  return n;
}

static inline void rsimd_soft_fma_mul(uint64_t *hi, uint64_t *lo, uint64_t x, uint64_t y) {
  uint64_t t1, t2, t3;
  uint64_t xlo = (uint32_t) x, xhi = x >> 32;
  uint64_t ylo = (uint32_t) y, yhi = y >> 32;

  t1 = xlo * ylo;
  t2 = xlo * yhi + xhi * ylo;
  t3 = xhi * yhi;
  *lo = t1 + (t2 << 32);
  *hi = t3 + (t2 >> 32) + (t1 > *lo);
}

static inline int rsimd_soft_fma_clz64(uint64_t x) {
  int n = 0;
  if (!(x >> 32)) n += 32, x <<= 32;
  if (!(x >> 48)) n += 16, x <<= 16;
  if (!(x >> 56)) n += 8, x <<= 8;
  if (!(x >> 60)) n += 4, x <<= 4;
  if (!(x >> 62)) n += 2, x <<= 2;
  if (!(x >> 63)) n += 1;
  return n;
}

static inline double rsimd_soft_fma(double x, double y, double z) {
  /* normalize so top 10bits and last bit are 0 */
  struct rsimd_soft_fma_num nx, ny, nz;
  nx = rsimd_soft_fma_normalize(x);
  ny = rsimd_soft_fma_normalize(y);
  nz = rsimd_soft_fma_normalize(z);

  if (nx.e >= RSIMD_SOFT_FMA_ZEROINFNAN || ny.e >= RSIMD_SOFT_FMA_ZEROINFNAN) return x * y + z;
  if (nz.e >= RSIMD_SOFT_FMA_ZEROINFNAN) {
    if (nz.e > RSIMD_SOFT_FMA_ZEROINFNAN) /* z==0 */
      return x * y;
    return z;
  }

  /* mul: r = x*y */
  uint64_t rhi, rlo, zhi, zlo;
  rsimd_soft_fma_mul(&rhi, &rlo, nx.m, ny.m);
  /* either top 20 or 21 bits of rhi and last 2 bits of rlo are 0 */

  /* align exponents */
  int e = nx.e + ny.e;
  int d = nz.e - e;
  /* shift bits z<<=kz, r>>=kr, so kz+kr == d, set e = e+kr (== ez-kz) */
  if (d > 0) {
    if (d < 64) {
      zlo = nz.m << d;
      zhi = nz.m >> (64 - d);
    } else {
      zlo = 0;
      zhi = nz.m;
      e = nz.e - 64;
      d -= 64;
      if (d == 0) {
      } else if (d < 64) {
        rlo = rhi << (64 - d) | rlo >> d | !!(rlo << (64 - d));
        rhi = rhi >> d;
      } else {
        rlo = 1;
        rhi = 0;
      }
    }
  } else {
    zhi = 0;
    d = -d;
    if (d == 0) {
      zlo = nz.m;
    } else if (d < 64) {
      zlo = nz.m >> d | !!(nz.m << (64 - d));
    } else {
      zlo = 1;
    }
  }

  /* add */
  int sign = nx.sign ^ ny.sign;
  int samesign = !(sign ^ nz.sign);
  int nonzero = 1;
  if (samesign) {
    /* r += z */
    rlo += zlo;
    rhi += zhi + (rlo < zlo);
  } else {
    /* r -= z */
    uint64_t t = rlo;
    rlo -= zlo;
    rhi = rhi - zhi - (t < rlo);
    if (rhi >> 63) {
      rlo = -rlo;
      rhi = -rhi - !!rlo;
      sign = !sign;
    }
    nonzero = !!rhi;
  }

  /* set rhi to top 63bit of the result (last bit is sticky) */
  if (nonzero) {
    e += 64;
    d = rsimd_soft_fma_clz64(rhi) - 1;
    /* note: d > 0 */
    rhi = rhi << d | rlo >> (64 - d) | !!(rlo << d);
  } else if (rlo) {
    d = rsimd_soft_fma_clz64(rlo) - 1;
    if (d < 0)
      rhi = rlo >> 1 | (rlo & 1);
    else
      rhi = rlo << d;
  } else {
    /* exact +-0 */
    return x * y + z;
  }
  e -= d;

  /* convert to double */
  int64_t i = (int64_t) rhi; /* i is in [1<<62,(1<<63)-1] */
  if (sign) i = -i;
  double r = (double) i; /* |r| is in [0x1p62,0x1p63] */

  if (e < -1022 - 62) {
    /* result is subnormal before rounding */
    if (e == -1022 - 63) {
      double c = 0x1p63;
      if (sign) c = -c;
      if (r == c) {
        /* min normal after rounding, underflow depends
           on arch behaviour which can be imitated by
           a double to float conversion */
        float fltmin = (float) (0x0.ffffff8p-63 * FLT_MIN * r);
        return DBL_MIN / FLT_MIN * fltmin;
      }
      /* one bit is lost when scaled, add another top bit to
         only round once at conversion if it is inexact */
      if (rhi << 53) {
        i = (int64_t) (rhi >> 1 | (rhi & 1) | 1ull << 62);
        if (sign) i = -i;
        r = (double) i;
        r = 2 * r - c; /* remove top bit */
      }
    } else {
      /* only round once when scaled */
      d = 10;
      i = (int64_t) ((rhi >> d | !!(rhi << (64 - d))) << d);
      if (sign) i = -i;
      r = (double) i;
    }
  }
  return ldexp(r, e);
}

#if defined(RSIMD_SOFT_FMA) || (defined(_WIN32) && !defined(__FMA__))
#define rsimd_fma rsimd_soft_fma
#else
#define rsimd_fma fma
#endif

#endif
