# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-11T22:03:29Z |
| rsimd | 0.1.0 (git d95c7863d200e98e30b6f9679a28d91f9d7d1615) |
| rsimd library | `/Users/runner/work/_temp/Library/rsimd` (built R 4.6.1; aarch64-apple-darwin23; 2026-10-11 21:53:30 UTC; unix) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `clang -arch arm64` |
| CFLAGS | `-falign-functions=64 -Wall -g -O2` |
| OS | Darwin 25.6.0 (arm64) |
| CPU | Apple M1 (Virtual), 3 cores |
| CPU features | neon, fp16, dotprod |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 550 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.943 µs (3.5× n, 1.9× b) | 3.32 µs (0.54× b) | 1.8 µs |
| 1e5 | 17.7 µs (28.1× n, 13.7× b) | 496 µs (0.49× b) | 243 µs |
| 1e7 | 2.77 ms (11.7× n, 7.5× b) | 32.5 ms (0.64× b) | 20.8 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.03 µs (1.2× n, 0.80× b) | 1.27 µs (0.65× b) | 0.82 µs |
| 1e5 | 26.2 µs (2.3× n, 4.5× b) | 59.1 µs (2.0× b) | 118 µs |
| 1e7 | 3.03 ms (2.0× n, 2.6× b) | 5.91 ms (1.3× b) | 7.78 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.943 µs (3.4× n, 4.0× b) | 3.2 µs (1.2× b) | 3.77 µs |
| 1e5 | 17.8 µs (17.6× n, 12.8× b) | 314 µs (0.72× b) | 227 µs |
| 1e7 | 3.22 ms (10.3× n, 8.1× b) | 33.3 ms (0.79× b) | 26.2 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 1.27 µs (3.3× n, 1.7× b) | 4.14 µs (0.52× b) | 2.17 µs | 2.3 µs (0.95× b) |
| 1e5 | 22.5 µs (14.0× n, 10.2× b) | 315 µs (0.73× b) | 229 µs | 227 µs (1.0× b) |
| 1e7 | 5.11 ms (8.0× n, 8.5× b) | 41 ms (1.1× b) | 43.5 ms | 31 ms (1.4× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.6 µs (1.2× n, 0.28× b) | 1.97 µs (0.23× b) | 0.451 µs |
| 1e5 | 39.8 µs (2.7× n, 0.80× b) | 107 µs (0.30× b) | 31.9 µs |
| 1e7 | 77.8 ms (0.19× n, 0.11× b) | 14.9 ms (0.55× b) | 8.22 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.64 µs (1.8× n, 1.8× b) | 3.03 µs (1.0× b) | 2.99 µs |
| 1e5 | 38.3 µs (4.8× n, 7.8× b) | 184 µs (1.6× b) | 298 µs |
| 1e7 | 28.2 ms (4.0× n, 2.2× b) | 114 ms (0.54× b) | 60.9 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.8 µs (1.5× n, 0.61× b) | 2.66 µs (0.42× b) | 1.11 µs |
| 1e5 | 46.2 µs (3.7× n, 1.8× b) | 169 µs (0.50× b) | 84.2 µs |
| 1e7 | 28.4 ms (1.4× n, 0.86× b) | 39.7 ms (0.62× b) | 24.5 ms |

### scalar_sub (double)

`simd_sub(1, x)` vs base `1 - x`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.68 µs (1.3× n, 0.24× b) | 2.26 µs (0.18× b) | 0.41 µs |
| 1e5 | 44 µs (3.0× n, 0.47× b) | 133 µs (0.16× b) | 20.7 µs |
| 1e7 | 8.62 ms (2.5× n, 0.60× b) | 21.1 ms (0.24× b) | 5.15 ms |

### fma_scalar (double)

`simd_fma(x, 2, 1)` vs base `x * 2 + 1`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.76 µs (1.5× n, 0.47× b) | 2.62 µs (0.31× b) | 0.82 µs |
| 1e5 | 46.3 µs (3.4× n, 1.4× b) | 158 µs (0.42× b) | 65.8 µs |
| 1e7 | 6.47 ms (3.6× n, 1.9× b) | 23.4 ms (0.52× b) | 12.1 ms |

### lt_scalar (double)

`simd_lt(x, 0)` vs base `x < 0`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.48 µs (1.7× n, 0.94× b) | 2.5 µs (0.56× b) | 1.39 µs |
| 1e5 | 34.1 µs (4.6× n, 3.9× b) | 156 µs (0.85× b) | 133 µs |
| 1e7 | 7.1 ms (3.3× n, 2.8× b) | 23.6 ms (0.85× b) | 20.1 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.72 µs (1.1× n, 3.1× b) | 1.89 µs (2.8× b) | 5.37 µs |
| 1e5 | 37.2 µs (1.0× n, 18.7× b) ⚠ | 38.5 µs (18.1× b) | 697 µs |
| 1e7 | 23 ms (0.32× n, 3.4× b) ⚠ | 7.41 ms (10.6× b) | 78.2 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.67 µs (1.4× n, 1.4× b) | 3.69 µs (1.0× b) | 3.81 µs |
| 1e5 | 217 µs (3.1× n, 3.6× b) | 669 µs (1.2× b) | 784 µs |
| 1e7 | 25.1 ms (3.1× n, 3.4× b) | 77.3 ms (1.1× b) | 85.3 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.533 µs (1.5× n, 0.85× b) | 0.779 µs (0.58× b) | 0.451 µs |
| 1e5 | 12.1 µs (3.0× n, 2.8× b) | 36.4 µs (0.94× b) | 34.3 µs |
| 1e7 | 2.15 ms (1.8× n, 1.9× b) | 3.79 ms (1.1× b) | 4.08 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.492 µs (1.7× n, 1.2× b) | 0.82 µs (0.70× b) | 0.574 µs |
| 1e5 | 6.11 µs (6.1× n, 8.0× b) | 37.4 µs (1.3× b) | 48.7 µs |
| 1e7 | 1.2 ms (3.8× n, 4.4× b) | 4.57 ms (1.2× b) | 5.3 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (2.4× n, 0.71× b) | 2.38 µs (0.29× b) | 0.697 µs |
| 1e5 | 17.6 µs (11.1× n, 3.6× b) | 195 µs (0.32× b) | 63.3 µs |
| 1e7 | 13.4 ms (2.1× n, 0.59× b) | 27.6 ms (0.28× b) | 7.86 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.52 µs (1.1× n, 1.1× b) | 1.6 µs (1.0× b) | 1.64 µs |
| 1e5 | 94.5 µs (1.2× n, 4.1× b) | 117 µs (3.3× b) | 383 µs |
| 1e7 | 15.8 ms (1.2× n, 4.2× b) | 19.6 ms (3.4× b) | 65.8 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.779 µs (1.5× n, 1.1× b) | 1.19 µs (0.72× b) | 0.861 µs |
| 1e5 | 12 µs (4.6× n, 6.3× b) | 55.1 µs (1.4× b) | 75.8 µs |
| 1e7 | 2.47 ms (2.3× n, 3.2× b) | 5.65 ms (1.4× b) | 7.83 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.11 µs (2.1× n, 1.2× b) | 2.34 µs (0.58× b) | 1.35 µs |
| 1e5 | 20.2 µs (7.1× n, 4.7× b) | 144 µs (0.66× b) | 95.2 µs |
| 1e7 | 2.83 ms (6.4× n, 10.8× b) | 18 ms (1.7× b) | 30.6 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.861 µs (2.5× n, 5.6× b) | 2.13 µs (2.3× b) | 4.84 µs |
| 1e5 | 19 µs (7.7× n, 25.2× b) | 147 µs (3.3× b) | 478 µs |
| 1e7 | 3.17 ms (5.5× n, 17.8× b) | 17.3 ms (3.2× b) | 56.2 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.24 µs (1.5× n, 0.28× b) | 4.76 µs (0.19× b) | 0.902 µs |
| 1e5 | 74.4 µs (2.5× n, 1.1× b) | 189 µs (0.42× b) | 80.1 µs |
| 1e7 | 23.9 ms (1.3× n, 0.58× b) | 31.1 ms (0.45× b) | 13.9 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.68 µs (1.6× n, 1.3× b) | 2.71 µs (0.79× b) | 2.13 µs |
| 1e5 | 38.2 µs (3.7× n, 5.4× b) | 142 µs (1.5× b) | 208 µs |
| 1e7 | 4.15 ms (4.3× n, 7.5× b) | 17.9 ms (1.7× b) | 31.2 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.6 µs (2.2× n, 1.6× b) | 3.49 µs (0.72× b) | 2.5 µs |
| 1e5 | 25.3 µs (9.1× n, 9.1× b) | 230 µs (1.0× b) | 231 µs |
| 1e7 | 2.75 ms (8.4× n, 10.3× b) | 23.2 ms (1.2× b) | 28.4 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.11 µs (2.1× n, 2.3× b) | 2.38 µs (1.1× b) | 2.58 µs |
| 1e5 | 37.1 µs (4.7× n, 6.8× b) | 174 µs (1.5× b) | 253 µs |
| 1e7 | 6.13 ms (2.7× n, 4.2× b) | 16.7 ms (1.6× b) | 26 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.23 µs (3.6× n, 12.0× b) | 4.43 µs (3.3× b) | 14.7 µs |
| 1e5 | 42.4 µs (18.8× n, 39.1× b) | 798 µs (2.1× b) | 1.66 ms |
| 1e7 | 4.97 ms (16.8× n, 42.4× b) | 83.6 ms (2.5× b) | 211 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.27 µs (0.84× n, 0.35× b) | 1.07 µs (0.42× b) | 0.451 µs |
| 1e5 | 45.6 µs (0.56× n, 0.51× b) | 25.3 µs (0.92× b) | 23.2 µs |
| 1e7 | 9.46 ms (3.3× n, 0.66× b) | 31.5 ms (0.20× b) | 6.26 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.35 µs (1.0× n, 4.0× b) | 1.31 µs (4.1× b) | 5.37 µs |
| 1e5 | 66.2 µs (1.1× n, 9.3× b) | 74.3 µs (8.3× b) | 614 µs |
| 1e7 | 9.04 ms (0.77× n, 8.7× b) | 6.99 ms (11.3× b) | 78.9 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.31 µs (1.0× n, 4.1× b) | 1.31 µs (4.1× b) | 5.41 µs |
| 1e5 | 71.3 µs (1.0× n, 8.8× b) ⚠ | 72.1 µs (8.7× b) | 629 µs |
| 1e7 | 6.83 ms (0.95× n, 10.0× b) ⚠ | 6.47 ms (10.5× b) | 68 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 31.6 µs (1.1× n, 0.37× b) | 33.9 µs (0.35× b) | 11.8 µs |
| 1e5 | 2.9 ms (1.2× n, 0.44× b) | 3.37 ms (0.38× b) | 1.27 ms |
| 1e7 | 303 ms (1.1× n, 0.43× b) | 335 ms (0.39× b) | 129 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.38 µs (0.83× n, 0.64× b) | 1.97 µs (0.77× b) | 1.52 µs |
| 1e5 | 91.2 µs (1.0× n, 1.6× b) | 91 µs (1.6× b) | 148 µs |
| 1e7 | 139 ms (0.11× n, 0.12× b) | 14.8 ms (1.1× b) | 16.2 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 15 µs (1.0× n, 0.95× b) | 15.5 µs (0.92× b) | 14.2 µs |
| 1e5 | 1.73 ms (1.5× n, 1.0× b) | 2.57 ms (0.66× b) | 1.69 ms |
| 1e7 | 221 ms (1.0× n, 0.88× b) ⚠ | 227 ms (0.86× b) | 196 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.3 µs (4.0× n, 1.2× b) | 9.14 µs (0.30× b) | 2.79 µs |
| 1e5 | 143 µs (6.2× n, 1.9× b) | 880 µs (0.31× b) | 275 µs |
| 1e7 | 13.4 ms (5.9× n, 2.0× b) | 79.5 ms (0.34× b) | 27 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.45 µs (0.50× n, 0.38× b) | 2.71 µs (0.77× b) | 2.09 µs |
| 1e5 | 517 µs (0.42× n, 0.40× b) ⚠ | 216 µs (1.0× b) | 207 µs |
| 1e7 | 86.4 ms (0.39× n, 0.28× b) ⚠ | 33.9 ms (0.70× b) | 23.8 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.85 µs (1.6× n, 1.3× b) | 6.23 µs (0.80× b) | 4.96 µs |
| 1e5 | 210 µs (2.1× n, 2.3× b) | 449 µs (1.1× b) | 487 µs |
| 1e7 | 89.3 ms (0.49× n, 0.52× b) | 43.5 ms (1.1× b) | 46.4 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.23 µs (2.5× n, 1.4× b) | 3.03 µs (0.57× b) | 1.72 µs |
| 1e5 | 39.3 µs (6.3× n, 4.5× b) | 246 µs (0.71× b) | 175 µs |
| 1e7 | 3.57 ms (6.3× n, 4.4× b) | 22.4 ms (0.70× b) | 15.7 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.902 µs (2.1× n, 0.86× b) | 1.93 µs (0.40× b) | 0.779 µs |
| 1e5 | 22.7 µs (6.3× n, 3.0× b) | 144 µs (0.48× b) | 68.7 µs |
| 1e7 | 1.84 ms (7.5× n, 3.7× b) | 13.8 ms (0.49× b) | 6.72 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.82 µs (3.1× n, 1.4× b) | 2.58 µs (0.43× b) | 1.11 µs |
| 1e5 | 14.9 µs (14.4× n, 7.4× b) | 214 µs (0.52× b) | 110 µs |
| 1e7 | 1.36 ms (13.9× n, 7.8× b) | 18.8 ms (0.56× b) | 10.6 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.615 µs (4.2× n, 4.0× b) | 2.58 µs (1.0× b) | 2.46 µs |
| 1e5 | 28.2 µs (10.3× n, 8.4× b) | 290 µs (0.82× b) | 237 µs |
| 1e7 | 2.3 ms (11.5× n, 16.5× b) | 26.4 ms (1.4× b) | 38 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.902 µs (2.0× n, 1.5× b) | 1.85 µs (0.73× b) | 1.35 µs |
| 1e5 | 5.99 µs (18.9× n, 25.8× b) | 113 µs (1.4× b) | 155 µs |
| 1e7 | 2.31 ms (4.7× n, 6.9× b) | 10.9 ms (1.5× b) | 15.9 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.779 µs (1.1× n, 20.2× b) | 0.861 µs (18.2× b) | 15.7 µs |
| 1e5 | 8.24 µs (1.2× n, 204.9× b) | 9.72 µs (173.8× b) | 1.69 ms |
| 1e7 | 604 µs (1.2× n, 186.2× b) | 711 µs (158.1× b) | 112 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.03 µs (3.2× n, 1.8× b) | 3.28 µs (0.55× b) | 1.8 µs |
| 1e5 | 1.56 µs (8.0× n, 117.4× b) | 12.5 µs (14.7× b) | 183 µs |
| 1e7 | 1.31 µs (8.2× n, 12474.0× b) | 10.7 µs (1523.5× b) | 16.4 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (0.83× n, 0.21× b) | 0.82 µs (0.25× b) | 0.205 µs |
| 1e5 | 1.8 µs (0.43× n, 0.11× b) | 0.779 µs (0.26× b) | 0.205 µs |
| 1e7 | 1.76 µs (0.47× n, 0.14× b) | 0.82 µs (0.30× b) | 0.246 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.492 µs (1.0× n, 0.25× b) | 0.492 µs (0.25× b) | 0.123 µs |
| 1e5 | 0.451 µs (1.2× n, 0.18× b) | 0.533 µs (0.15× b) | 0.082 µs |
| 1e7 | 0.451 µs (1.0× n, 0.27× b) | 0.451 µs (0.27× b) | 0.123 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.451 µs (1.0× n, 0.27× b) | 0.451 µs (0.27× b) | 0.123 µs |
| 1e5 | 0.41 µs (1.1× n, 0.20× b) | 0.451 µs (0.18× b) | 0.082 µs |
| 1e7 | 0.41 µs (1.1× n, 0.30× b) | 0.451 µs (0.27× b) | 0.123 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.15 µs (2.9× n, 1.5× b) | 3.28 µs (0.54× b) | 1.76 µs |
| 1e5 | 36.3 µs (7.7× n, 4.9× b) | 279 µs (0.63× b) | 177 µs |
| 1e7 | 3.63 ms (6.8× n, 4.6× b) | 24.7 ms (0.67× b) | 16.6 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.07 µs (2.0× n, 0.81× b) | 2.13 µs (0.40× b) | 0.861 µs |
| 1e5 | 1.56 µs (3.9× n, 50.6× b) | 6.15 µs (12.8× b) | 78.8 µs |
| 1e7 | 1.56 µs (3.9× n, 4783.6× b) | 6.15 µs (1211.8× b) | 7.45 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (2.8× n, 0.25× b) | 2.79 µs (0.09× b) | 0.246 µs |
| 1e5 | 1.44 µs (6.3× n, 0.11× b) | 8.98 µs (0.02× b) | 0.164 µs |
| 1e7 | 1.23 µs (6.8× n, 0.17× b) | 8.36 µs (0.02× b) | 0.205 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.35 ms (16.1× n, 11.5× b) | 21.8 ms (0.71× b) | 15.5 ms |
| pairwise | 1.41 ms (16.8× n, 11.2× b) | 23.6 ms (0.67× b) | 15.7 ms |
| compensated | 3.94 ms (5.0× n, 3.9× b) | 19.5 ms (0.80× b) | 15.5 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 3.59 ms (7.8× n, 9.4× b) | 28.1 ms (1.2× b) | 33.6 ms | 19.7 ms (1.7× b) |
| pairwise | 2.91 ms (10.0× n, 7.6× b) | 29 ms (0.76× b) | 22.2 ms | 21.6 ms (1.0× b) |
| compensated | 4.58 ms (5.6× n, 4.6× b) | 25.4 ms (0.83× b) | 21.2 ms | 21.2 ms (1.0× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 76.1 ms (1.1× n, 1.1× b) | 80.4 ms (1.0× b) | 83.7 ms |
| fast | 44.1 ms (1.8× n, 1.9× b) | 79.5 ms (1.1× b) | 84.4 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 35.5 ms (0.88× n, 1.0× b) | 31.4 ms (1.1× b) | 35.1 ms |
| fast | 28.7 ms (1.6× n, 1.4× b) | 45.4 ms (0.89× b) | 40.4 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 89.1 ms (0.93× n, 0.91× b) | 82.7 ms (1.0× b) | 81.5 ms |
| fast | 34.2 ms (2.3× n, 2.6× b) | 78.1 ms (1.1× b) | 87.5 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 73.4 ms (3.7× n, 3.8× b) | 269 ms (1.0× b) | 279 ms |
| fast | 34.9 ms (7.4× n, 8.0× b) | 257 ms (1.1× b) | 280 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 52.9 ms (0.62× n, 1.0× b) | 33 ms (1.7× b) | 55 ms |
| fast | 14.1 ms (2.7× n, 4.3× b) | 37.9 ms (1.6× b) | 60.8 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 122 ms (1.1× n, 1.5× b) | 136 ms (1.4× b) | 189 ms |
| fast | 166 ms (0.89× n, 1.1× b) | 147 ms (1.2× b) | 181 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 61.5 ms (1.1× n, 0.87× b) | 64.7 ms (0.82× b) | 53.3 ms |
| fast | 54 ms (1.0× n, 1.0× b) | 56.1 ms (0.94× b) | 52.6 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 129 ms (0.32× n, 0.32× b) | 40.9 ms (1.0× b) | 41.8 ms |
| fast | 32 ms (1.2× n, 1.0× b) | 37.6 ms (0.86× b) | 32.4 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 81 ms (2.6× n, 2.7× b) | 214 ms (1.0× b) | 216 ms |
| fast | 59.9 ms (4.6× n, 5.1× b) | 277 ms (1.1× b) | 306 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 342 ms (1.4× n, 1.5× b) | 469 ms (1.1× b) | 519 ms |
| fast | 154 ms (2.7× n, 3.0× b) | 408 ms (1.1× b) | 457 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 138 ms (1.4× n, 1.3× b) | 190 ms (0.92× b) | 175 ms |
| fast | 115 ms (1.6× n, 1.5× b) | 183 ms (1.0× b) | 176 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 268 ms (2.2× n, 2.2× b) | 594 ms (1.0× b) | 602 ms |
| fast | 237 ms (2.8× n, 2.5× b) | 655 ms (0.90× b) | 588 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 137 ms (1.0× n, 1.1× b) | 143 ms (1.1× b) | 152 ms |
| fast | 119 ms (1.2× n, 1.3× b) | 143 ms (1.1× b) | 157 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 451 ms (1.6× n, 1.8× b) | 707 ms (1.1× b) | 796 ms |
| fast | 321 ms (2.3× n, 3.0× b) | 742 ms (1.3× b) | 963 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.164 µs | 0.779 µs (0.21× b) | 0.533 µs (0.31× b) |
| 1e1 | 0.164 µs | 0.779 µs (0.21× b) | 0.615 µs (0.27× b) |
| 1e2 | 0.328 µs | 0.861 µs (0.38× b) | 0.615 µs (0.53× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.041 µs | 0.902 µs (0.05× b) |
| 1e1 | 0.082 µs | 1.07 µs (0.08× b) |
| 1e2 | 0.164 µs | 1.19 µs (0.14× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.041 µs | 0.697 µs (0.06× b) |
| 1e1 | 0.082 µs | 0.697 µs (0.12× b) |
| 1e2 | 0.492 µs | 0.984 µs (0.50× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.041 µs | 0.984 µs (0.04× b) |
| 1e1 | 0.123 µs | 0.984 µs (0.13× b) |
| 1e2 | 0.287 µs | 1.15 µs (0.25× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.041 µs | 0.943 µs (0.04× b) |
| 1e1 | 0.082 µs | 1.11 µs (0.07× b) |
| 1e2 | 0.205 µs | 1.31 µs (0.16× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.205 µs | 1.11 µs (0.19× b) |
| 1e1 | 0.287 µs | 1.23 µs (0.23× b) |
| 1e2 | 0.574 µs | 1.31 µs (0.44× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.082 µs | 0.984 µs (0.08× b) |
| 1e1 | 0.123 µs | 1.07 µs (0.12× b) |
| 1e2 | 0.287 µs | 2.79 µs (0.10× b) |

### add (integer64)

`simd_add(xi64, yi64)` vs base `xi64 + yi64`

| n | base | auto |
|---|---|---|
| 1e0 | 31.7 µs | 2.99 µs (10.6× b) |
| 1e1 | 32.2 µs | 3.16 µs (10.2× b) |
| 1e2 | 41.4 µs | 3.73 µs (11.1× b) |

