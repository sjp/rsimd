# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-07T06:49:13Z |
| rsimd | 0.1.0 (git 8a2d3e839cdb544674eb5fd938b645cac69469fd) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `clang -arch arm64` |
| CFLAGS | `-falign-functions=64 -Wall -g -O2` |
| OS | Darwin 25.6.0 (arm64) |
| CPU | Apple M1 (Virtual), 3 cores |
| CPU features | neon, fp16, dotprod |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 467 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (2.9× n, 2.3× b) | 2.83 µs (0.78× b) | 2.21 µs |
| 1e5 | 15.2 µs (13.7× n, 10.9× b) | 208 µs (0.79× b) | 165 µs |
| 1e7 | 1.87 ms (11.5× n, 9.2× b) | 21.5 ms (0.80× b) | 17.2 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.943 µs (1.5× n, 0.87× b) | 1.39 µs (0.59× b) | 0.82 µs |
| 1e5 | 13.4 µs (3.7× n, 4.9× b) | 50 µs (1.3× b) | 66.1 µs |
| 1e7 | 1.36 ms (3.8× n, 4.9× b) | 5.19 ms (1.3× b) | 6.74 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (3.0× n, 3.8× b) | 2.99 µs (1.3× b) | 3.77 µs |
| 1e5 | 14.9 µs (14.0× n, 13.5× b) | 209 µs (1.0× b) | 201 µs |
| 1e7 | 1.92 ms (11.1× n, 11.3× b) | 21.4 ms (1.0× b) | 21.8 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 1.31 µs (2.8× n, 1.8× b) | 3.69 µs (0.63× b) | 2.34 µs | 2.34 µs (1.0× b) |
| 1e5 | 23.8 µs (11.3× n, 8.5× b) | 268 µs (0.76× b) | 203 µs | 211 µs (1.0× b) |
| 1e7 | 3.37 ms (8.4× n, 9.8× b) | 28.2 ms (1.2× b) | 33.2 ms | 22.9 ms (1.4× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.93 µs (1.5× n, 0.46× b) | 2.93 µs (0.30× b) | 0.882 µs |
| 1e5 | 33.5 µs (4.0× n, 0.83× b) | 135 µs (0.21× b) | 27.8 µs |
| 1e7 | 5.41 ms (2.7× n, 1.0× b) | 14.4 ms (0.37× b) | 5.27 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.44 µs (1.8× n, 2.0× b) | 2.58 µs (1.1× b) | 2.87 µs |
| 1e5 | 36.8 µs (4.2× n, 7.4× b) | 153 µs (1.8× b) | 272 µs |
| 1e7 | 4.71 ms (3.4× n, 5.4× b) | 16.1 ms (1.6× b) | 25.6 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.6 µs (2.0× n, 0.72× b) | 3.2 µs (0.36× b) | 1.15 µs |
| 1e5 | 39.9 µs (5.1× n, 1.7× b) | 202 µs (0.34× b) | 68.8 µs |
| 1e7 | 7.88 ms (2.6× n, 1.4× b) | 20.5 ms (0.53× b) | 10.8 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.84 µs (1.6× n, 3.1× b) | 2.93 µs (1.9× b) | 5.66 µs |
| 1e5 | 30.2 µs (3.9× n, 19.8× b) | 117 µs (5.1× b) | 600 µs |
| 1e7 | 4.78 ms (2.5× n, 22.8× b) | 12 ms (9.1× b) | 109 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.87 µs (1.1× n, 1.2× b) | 3.28 µs (1.1× b) | 3.44 µs |
| 1e5 | 213 µs (3.0× n, 3.3× b) | 646 µs (1.1× b) | 702 µs |
| 1e7 | 27.2 ms (2.7× n, 2.8× b) | 72.3 ms (1.1× b) | 76.6 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.697 µs (1.3× n, 0.59× b) | 0.902 µs (0.45× b) | 0.41 µs |
| 1e5 | 12.3 µs (2.8× n, 2.6× b) | 34.2 µs (0.95× b) | 32.3 µs |
| 1e7 | 1.9 ms (1.8× n, 2.0× b) | 3.34 ms (1.1× b) | 3.71 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.574 µs (1.6× n, 0.86× b) | 0.902 µs (0.55× b) | 0.492 µs |
| 1e5 | 5.86 µs (6.1× n, 7.3× b) | 36 µs (1.2× b) | 42.7 µs |
| 1e7 | 1.17 ms (3.6× n, 3.8× b) | 4.26 ms (1.1× b) | 4.48 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.11 µs (2.1× n, 0.56× b) | 2.3 µs (0.27× b) | 0.615 µs |
| 1e5 | 31 µs (5.5× n, 1.7× b) | 170 µs (0.30× b) | 51.7 µs |
| 1e7 | 27 ms (0.83× n, 0.21× b) | 22.4 ms (0.26× b) | 5.77 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.44 µs (1.1× n, 1.0× b) | 1.56 µs (0.93× b) | 1.46 µs |
| 1e5 | 74.5 µs (1.1× n, 3.6× b) | 82.9 µs (3.3× b) | 271 µs |
| 1e7 | 9.06 ms (4.4× n, 5.6× b) | 39.9 ms (1.3× b) | 51.1 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.738 µs (1.6× n, 1.1× b) | 1.15 µs (0.71× b) | 0.82 µs |
| 1e5 | 10.7 µs (4.7× n, 5.9× b) | 50.3 µs (1.3× b) | 63.8 µs |
| 1e7 | 1.97 ms (3.2× n, 4.0× b) | 6.23 ms (1.3× b) | 7.91 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.11 µs (2.0× n, 1.1× b) | 2.21 µs (0.53× b) | 1.17 µs |
| 1e5 | 22.7 µs (5.7× n, 3.6× b) | 129 µs (0.64× b) | 82.2 µs |
| 1e7 | 2.45 ms (5.8× n, 15.8× b) | 14.3 ms (2.7× b) | 38.6 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.902 µs (2.5× n, 4.8× b) | 2.25 µs (1.9× b) | 4.35 µs |
| 1e5 | 22.3 µs (6.0× n, 18.8× b) | 134 µs (3.1× b) | 420 µs |
| 1e7 | 2.53 ms (5.5× n, 18.5× b) | 14 ms (3.3× b) | 46.7 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.3 µs (1.4× n, 0.38× b) | 3.32 µs (0.26× b) | 0.861 µs |
| 1e5 | 66.2 µs (2.5× n, 1.1× b) | 169 µs (0.42× b) | 71.2 µs |
| 1e7 | 31.3 ms (0.61× n, 0.24× b) | 19.2 ms (0.38× b) | 7.38 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.76 µs (1.6× n, 1.3× b) | 2.83 µs (0.83× b) | 2.34 µs |
| 1e5 | 37.4 µs (3.6× n, 5.0× b) | 135 µs (1.4× b) | 189 µs |
| 1e7 | 4.12 ms (3.3× n, 5.3× b) | 13.5 ms (1.6× b) | 21.8 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.64 µs (2.1× n, 1.6× b) | 3.49 µs (0.75× b) | 2.62 µs |
| 1e5 | 24.7 µs (8.8× n, 8.9× b) | 219 µs (1.0× b) | 219 µs |
| 1e7 | 2.43 ms (9.2× n, 11.7× b) | 22.2 ms (1.3× b) | 28.5 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.23 µs (2.1× n, 2.0× b) | 2.54 µs (1.0× b) | 2.46 µs |
| 1e5 | 42.7 µs (4.0× n, 5.2× b) | 172 µs (1.3× b) | 224 µs |
| 1e7 | 9.68 ms (1.9× n, 3.0× b) | 18.6 ms (1.6× b) | 28.9 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.23 µs (3.6× n, 11.6× b) | 4.47 µs (3.2× b) | 14.3 µs |
| 1e5 | 44.8 µs (17.2× n, 33.2× b) | 773 µs (1.9× b) | 1.49 ms |
| 1e7 | 4.7 ms (17.7× n, 33.3× b) | 83.4 ms (1.9× b) | 157 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.27 µs (1.0× n, 0.55× b) | 1.31 µs (0.53× b) | 0.697 µs |
| 1e5 | 45.9 µs (0.50× n, 0.44× b) | 22.8 µs (0.89× b) | 20.3 µs |
| 1e7 | 4.38 ms (0.71× n, 0.72× b) | 3.13 ms (1.0× b) | 3.17 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.03 µs (1.0× n, 2.7× b) | 1.97 µs (2.8× b) | 5.45 µs |
| 1e5 | 60 µs (1.0× n, 8.6× b) | 61.7 µs (8.3× b) | 513 µs |
| 1e7 | 5.89 ms (1.1× n, 8.8× b) | 6.2 ms (8.4× b) | 52.1 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.35 µs (1.0× n, 3.8× b) | 1.39 µs (3.7× b) | 5.21 µs |
| 1e5 | 59.8 µs (1.1× n, 9.0× b) ⚠ | 63 µs (8.5× b) | 537 µs |
| 1e7 | 6.38 ms (1.0× n, 8.9× b) ⚠ | 6.69 ms (8.5× b) | 57 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 30.7 µs (1.2× n, 0.38× b) | 35.9 µs (0.32× b) | 11.5 µs |
| 1e5 | 2.68 ms (1.1× n, 0.47× b) | 2.98 ms (0.43× b) | 1.27 ms |
| 1e7 | 290 ms (1.1× n, 0.43× b) ⚠ | 310 ms (0.40× b) | 125 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.63 µs (1.3× n, 0.50× b) | 6.11 µs (0.38× b) | 2.34 µs |
| 1e5 | 87.7 µs (2.4× n, 1.6× b) | 211 µs (0.65× b) | 138 µs |
| 1e7 | 98.9 ms (0.23× n, 0.14× b) | 22.6 ms (0.60× b) | 13.6 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 17.6 µs (1.0× n, 0.82× b) | 18.1 µs (0.79× b) | 14.4 µs |
| 1e5 | 1.39 ms (1.0× n, 1.0× b) ⚠ | 1.46 ms (1.0× b) | 1.43 ms |
| 1e7 | 143 ms (1.1× n, 1.1× b) ⚠ | 150 ms (1.0× b) | 156 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.38 µs (3.8× n, 1.2× b) | 9.06 µs (0.31× b) | 2.79 µs |
| 1e5 | 140 µs (6.0× n, 1.9× b) | 844 µs (0.31× b) | 266 µs |
| 1e7 | 17 ms (5.0× n, 1.6× b) | 85.2 ms (0.32× b) | 27.2 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.78 µs (0.50× n, 0.40× b) | 2.91 µs (0.79× b) | 2.3 µs |
| 1e5 | 492 µs (0.41× n, 0.41× b) ⚠ | 200 µs (1.0× b) | 200 µs |
| 1e7 | 139 ms (0.17× n, 0.16× b) ⚠ | 24.3 ms (0.90× b) | 22 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.22 µs (2.2× n, 1.1× b) | 9.31 µs (0.52× b) | 4.84 µs |
| 1e5 | 183 µs (4.0× n, 2.6× b) | 725 µs (0.65× b) | 474 µs |
| 1e7 | 24.8 ms (3.2× n, 2.1× b) | 79.3 ms (0.66× b) | 52.4 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.15 µs (2.5× n, 1.6× b) | 2.83 µs (0.65× b) | 1.84 µs |
| 1e5 | 34.8 µs (5.8× n, 4.8× b) | 202 µs (0.82× b) | 166 µs |
| 1e7 | 3.65 ms (5.7× n, 5.0× b) | 20.8 ms (0.88× b) | 18.3 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (2.1× n, 0.79× b) | 2.09 µs (0.37× b) | 0.779 µs |
| 1e5 | 18.8 µs (7.1× n, 3.5× b) | 134 µs (0.50× b) | 66.3 µs |
| 1e7 | 2.59 ms (5.6× n, 2.7× b) | 14.4 ms (0.49× b) | 7.07 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.861 µs (3.1× n, 1.2× b) | 2.67 µs (0.40× b) | 1.07 µs |
| 1e5 | 13.7 µs (14.6× n, 7.4× b) | 199 µs (0.51× b) | 102 µs |
| 1e7 | 1.42 ms (15.3× n, 7.4× b) | 21.6 ms (0.48× b) | 10.5 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.779 µs (3.7× n, 2.6× b) | 2.91 µs (0.70× b) | 2.05 µs |
| 1e5 | 22.5 µs (10.4× n, 8.4× b) | 233 µs (0.81× b) | 189 µs |
| 1e7 | 2.5 ms (11.3× n, 13.9× b) | 28.2 ms (1.2× b) | 34.8 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.44 µs (1.5× n, 0.86× b) | 2.21 µs (0.56× b) | 1.23 µs |
| 1e5 | 5.7 µs (17.9× n, 21.3× b) | 102 µs (1.2× b) | 121 µs |
| 1e7 | 3.35 ms (3.9× n, 4.2× b) | 13.1 ms (1.1× b) | 13.9 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (1.0× n, 13.5× b) | 0.943 µs (14.0× b) | 13.2 µs |
| 1e5 | 7.38 µs (1.1× n, 165.9× b) | 8.12 µs (150.8× b) | 1.22 ms |
| 1e7 | 638 µs (1.1× n, 219.6× b) ⚠ | 700 µs (200.0× b) | 140 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.11 µs (2.5× n, 1.6× b) | 2.79 µs (0.65× b) | 1.8 µs |
| 1e5 | 1.6 µs (7.2× n, 103.9× b) | 11.6 µs (14.4× b) | 166 µs |
| 1e7 | 1.48 µs (6.3× n, 13307.7× b) | 9.31 µs (2110.5× b) | 19.6 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.984 µs (0.87× n, 0.25× b) | 0.861 µs (0.29× b) | 0.246 µs |
| 1e5 | 1.31 µs (0.62× n, 0.16× b) | 0.82 µs (0.25× b) | 0.205 µs |
| 1e7 | 1.43 µs (0.66× n, 0.17× b) | 0.943 µs (0.26× b) | 0.246 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.574 µs (1.1× n, 0.21× b) | 0.615 µs (0.20× b) | 0.123 µs |
| 1e5 | 0.574 µs (1.1× n, 0.14× b) | 0.615 µs (0.13× b) | 0.082 µs |
| 1e7 | 0.533 µs (1.1× n, 0.23× b) | 0.574 µs (0.21× b) | 0.123 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.492 µs (1.0× n, 0.25× b) | 0.492 µs (0.25× b) | 0.123 µs |
| 1e5 | 0.533 µs (1.1× n, 0.15× b) | 0.574 µs (0.14× b) | 0.082 µs |
| 1e7 | 0.533 µs (1.2× n, 0.23× b) | 0.615 µs (0.20× b) | 0.123 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.25 µs (2.3× n, 1.4× b) | 2.83 µs (0.61× b) | 1.72 µs |
| 1e5 | 34.3 µs (6.1× n, 4.8× b) | 208 µs (0.80× b) | 166 µs |
| 1e7 | 3.61 ms (6.7× n, 5.0× b) | 24.3 ms (0.73× b) | 17.9 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.07 µs (1.9× n, 0.77× b) | 2.01 µs (0.41× b) | 0.82 µs |
| 1e5 | 1.56 µs (3.9× n, 48.2× b) | 6.11 µs (12.3× b) | 75 µs |
| 1e7 | 1.56 µs (4.1× n, 5327.1× b) | 6.42 µs (1293.5× b) | 8.3 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.861 µs (3.1× n, 0.24× b) | 2.67 µs (0.08× b) | 0.205 µs |
| 1e5 | 1.27 µs (6.9× n, 0.13× b) | 8.82 µs (0.02× b) | 0.164 µs |
| 1e7 | 1.27 µs (7.1× n, 0.19× b) | 8.98 µs (0.03× b) | 0.246 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.83 ms (12.3× n, 10.1× b) | 22.5 ms (0.82× b) | 18.5 ms |
| pairwise | 2.23 ms (12.1× n, 7.8× b) | 26.9 ms (0.65× b) | 17.5 ms |
| compensated | 5.62 ms (4.7× n, 3.2× b) | 26.5 ms (0.67× b) | 17.9 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 4.2 ms (7.8× n, 14.4× b) | 33 ms (1.8× b) | 60.5 ms | 24.1 ms (2.5× b) |
| pairwise | 3.97 ms (7.7× n, 12.9× b) | 30.3 ms (1.7× b) | 51 ms | 22.1 ms (2.3× b) |
| compensated | 6.76 ms (4.7× n, 3.5× b) | 31.7 ms (0.74× b) | 23.4 ms | 22.4 ms (1.0× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 80.3 ms (1.1× n, 1.2× b) | 89.9 ms (1.1× b) | 95.1 ms |
| fast | 50 ms (1.8× n, 1.9× b) | 88.5 ms (1.1× b) | 93.5 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 36.6 ms (1.0× n, 0.95× b) | 35.9 ms (1.0× b) | 34.7 ms |
| fast | 24.3 ms (1.4× n, 1.5× b) | 34.4 ms (1.1× b) | 36.6 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 98.9 ms (1.0× n, 0.88× b) | 101 ms (0.86× b) | 87.3 ms |
| fast | 36.2 ms (2.2× n, 2.4× b) | 80.7 ms (1.1× b) | 86.3 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 115 ms (2.6× n, 2.8× b) | 295 ms (1.1× b) | 321 ms |
| fast | 55.9 ms (5.0× n, 4.9× b) | 280 ms (1.0× b) | 275 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 52 ms (0.69× n, 1.1× b) | 35.8 ms (1.7× b) | 59.3 ms |
| fast | 12.8 ms (2.8× n, 4.8× b) | 36 ms (1.7× b) | 61.8 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 174 ms (1.0× n, 1.0× b) | 173 ms (1.0× b) | 169 ms |
| fast | 174 ms (1.0× n, 1.0× b) | 174 ms (1.0× b) | 169 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 60.9 ms (0.94× n, 0.93× b) | 57.5 ms (1.0× b) | 56.5 ms |
| fast | 57.5 ms (1.0× n, 1.0× b) | 57.4 ms (1.0× b) | 56.5 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 110 ms (0.43× n, 0.32× b) | 47.6 ms (0.73× b) | 34.9 ms |
| fast | 31.7 ms (1.1× n, 1.1× b) | 35.9 ms (0.94× b) | 33.8 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 86.4 ms (2.7× n, 2.7× b) | 230 ms (1.0× b) | 230 ms |
| fast | 58.2 ms (3.9× n, 3.9× b) | 230 ms (1.0× b) | 230 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 253 ms (1.6× n, 1.6× b) | 396 ms (1.0× b) | 405 ms |
| fast | 121 ms (3.3× n, 3.3× b) | 395 ms (1.0× b) | 402 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 135 ms (1.3× n, 1.3× b) | 181 ms (1.0× b) | 173 ms |
| fast | 113 ms (1.6× n, 1.7× b) | 179 ms (1.0× b) | 186 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 260 ms (2.2× n, 2.3× b) | 582 ms (1.0× b) | 607 ms |
| fast | 231 ms (2.6× n, 2.6× b) | 593 ms (1.0× b) | 598 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 131 ms (1.1× n, 1.0× b) | 141 ms (1.0× b) | 135 ms |
| fast | 107 ms (1.3× n, 1.3× b) | 141 ms (1.0× b) | 135 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 427 ms (1.7× n, 1.6× b) | 713 ms (0.93× b) | 667 ms |
| fast | 237 ms (3.5× n, 3.8× b) | 838 ms (1.1× b) | 909 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.164 µs | 0.902 µs (0.18× b) | 0.779 µs (0.21× b) |
| 1e1 | 0.164 µs | 0.861 µs (0.19× b) | 0.697 µs (0.24× b) |
| 1e2 | 0.943 µs | 0.902 µs (1.0× b) | 0.82 µs (1.2× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.082 µs | 1.15 µs (0.07× b) |
| 1e1 | 0.082 µs | 0.943 µs (0.09× b) |
| 1e2 | 0.164 µs | 1.44 µs (0.11× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.041 µs | 0.738 µs (0.06× b) |
| 1e1 | 0.082 µs | 0.82 µs (0.10× b) |
| 1e2 | 0.41 µs | 2.23 µs (0.18× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.246 µs | 0.943 µs (0.26× b) |
| 1e1 | 0.123 µs | 0.902 µs (0.14× b) |
| 1e2 | 0.246 µs | 1.15 µs (0.21× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.082 µs | 1.76 µs (0.05× b) |
| 1e1 | 0.082 µs | 1.21 µs (0.07× b) |
| 1e2 | 0.205 µs | 1.27 µs (0.16× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.246 µs | 1.19 µs (0.21× b) |
| 1e1 | 1.07 µs | 3.85 µs (0.28× b) |
| 1e2 | 0.656 µs | 3.89 µs (0.17× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.082 µs | 3.18 µs (0.03× b) |
| 1e1 | 0.164 µs | 9.92 µs (0.02× b) |
| 1e2 | 0.943 µs | 10.6 µs (0.09× b) |

