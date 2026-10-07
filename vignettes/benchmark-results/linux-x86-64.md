# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-07T06:57:38Z |
| rsimd | 0.1.0 (git 8a2d3e839cdb544674eb5fd938b645cac69469fd) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `gcc -std=gnu2x` |
| CFLAGS | `-g -O2` |
| OS | Linux 6.17.0-1022-azure (x86_64) |
| CPU | AMD EPYC 7763 64-Core Processor, 4 cores |
| CPU features | sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, fma, os_avx |
| Tiers | avx2, sse2, none (auto: avx2) |
| bench | 1.1.4 |
| Run time | 830 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.53 µs (2.6× n, 1.4× b) | 1.6 µs (2.4× n, 1.3× b) | 3.93 µs (0.55× b) | 2.16 µs |
| 1e5 | 20.2 µs (12.4× n, 9.2× b) | 24.8 µs (10.1× n, 7.5× b) | 250 µs (0.74× b) | 186 µs |
| 1e7 | 3.33 ms (7.5× n, 5.6× b) | 3.9 ms (6.4× n, 4.8× b) | 25.1 ms (0.75× b) | 18.8 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.44 µs (1.4× n, 0.63× b) | 1.51 µs (1.3× n, 0.60× b) | 1.95 µs (0.46× b) | 0.902 µs |
| 1e5 | 9.49 µs (6.0× n, 6.6× b) | 15.8 µs (3.6× n, 3.9× b) | 56.9 µs (1.1× b) | 62.2 µs |
| 1e7 | 1.85 ms (3.0× n, 3.4× b) | 2.13 ms (2.6× n, 3.0× b) | 5.48 ms (1.1× b) | 6.28 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.57 µs (2.5× n, 4.7× b) | 1.66 µs (2.4× n, 4.4× b) | 3.97 µs (1.9× b) | 7.38 µs |
| 1e5 | 17.1 µs (14.6× n, 25.7× b) | 24.4 µs (10.2× n, 18.0× b) | 250 µs (1.8× b) | 440 µs |
| 1e7 | 3.31 ms (7.6× n, 13.2× b) | 3.86 ms (6.5× n, 11.3× b) | 25.1 ms (1.7× b) | 43.7 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | avx2 | sse2 | none | base | blas |
|---|---|---|---|---|---|
| 1e3 | 2.15 µs (1.8× n, 2.8× b) | 2.3 µs (1.7× n, 2.6× b) | 3.88 µs (1.5× b) | 5.97 µs | 1.29 µs (4.6× b) |
| 1e5 | 29 µs (6.8× n, 17.9× b) | 33.2 µs (6.0× n, 15.6× b) | 197 µs (2.6× b) | 517 µs | 68.8 µs (7.5× b) |
| 1e7 | 5.09 ms (3.9× n, 6.1× b) | 5.63 ms (3.5× n, 5.5× b) | 19.7 ms (1.6× b) | 30.9 ms | 13.6 ms (2.3× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.32 µs (2.8× n, 1.7× b) | 5.62 µs (1.1× n, 0.71× b) | 6.4 µs (0.62× b) | 3.99 µs |
| 1e5 | 327 µs (1.3× n, 1.1× b) | 335 µs (1.2× n, 1.0× b) | 413 µs (0.84× b) | 347 µs |
| 1e7 | 11.3 ms (1.6× n, 1.1× b) | 11.3 ms (1.6× n, 1.1× b) | 17.9 ms (0.68× b) | 12.1 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.34 µs (1.4× n, 0.73× b) | 4.05 µs (0.81× n, 0.42× b) | 3.27 µs (0.52× b) | 1.71 µs |
| 1e5 | 45.7 µs (6.2× n, 13.9× b) | 164 µs (1.7× n, 3.9× b) | 285 µs (2.2× b) | 635 µs |
| 1e7 | 6.17 ms (2.8× n, 8.8× b) | 7.85 ms (2.2× n, 6.9× b) | 17.2 ms (3.1× b) | 54.1 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.74 µs (2.1× n, 0.61× b) | 4.34 µs (1.3× n, 0.38× b) | 5.84 µs (0.28× b) | 1.66 µs |
| 1e5 | 332 µs (2.0× n, 1.2× b) | 513 µs (1.3× n, 0.76× b) | 668 µs (0.59× b) | 392 µs |
| 1e7 | 13.7 ms (3.1× n, 1.3× b) | 27.2 ms (1.6× n, 0.67× b) | 42.9 ms (0.43× b) | 18.3 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 6.21 µs (0.64× n, 1.4× b) | 3.17 µs (1.2× n, 2.8× b) | 3.95 µs (2.2× b) | 8.85 µs |
| 1e5 | 325 µs (1.3× n, 3.0× b) | 328 µs (1.3× n, 3.0× b) | 415 µs (2.4× b) | 978 µs |
| 1e7 | 11 ms (1.6× n, 6.6× b) | 11.2 ms (1.6× n, 6.5× b) | 17.8 ms (4.1× b) | 72.8 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 3.29 µs (2.1× n, 1.8× b) | 5.99 µs (1.1× n, 1.0× b) | 6.74 µs (0.90× b) | 6.07 µs |
| 1e5 | 464 µs (1.7× n, 1.8× b) | 744 µs (1.1× n, 1.1× b) ⚠ | 791 µs (1.1× b) | 851 µs |
| 1e7 | 22.5 ms (2.5× n, 2.7× b) | 50.4 ms (1.1× n, 1.2× b) | 55.4 ms (1.1× b) | 61.7 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.17 µs (1.2× n, 0.34× b) | 1.17 µs (1.2× n, 0.34× b) | 1.38 µs (0.29× b) | 0.401 µs |
| 1e5 | 18.5 µs (1.8× n, 1.7× b) | 17.3 µs (1.9× n, 1.8× b) | 33.1 µs (0.95× b) | 31.4 µs |
| 1e7 | 3.05 ms (1.4× n, 1.4× b) | 3.64 ms (1.2× n, 1.2× b) | 4.35 ms (1.0× b) | 4.29 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.02 µs (1.7× n, 0.40× b) | 1.11 µs (1.6× n, 0.37× b) | 1.78 µs (0.23× b) | 0.411 µs |
| 1e5 | 7.81 µs (8.1× n, 4.0× b) | 8.85 µs (7.2× n, 3.5× b) | 63.3 µs (0.49× b) | 31 µs |
| 1e7 | 1.25 ms (5.0× n, 2.6× b) | 1.79 ms (3.5× n, 1.8× b) | 6.29 ms (0.51× b) | 3.23 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.64 µs (1.3× n, 0.49× b) | 1.58 µs (1.3× n, 0.51× b) | 2.07 µs (0.39× b) | 0.802 µs |
| 1e5 | 39.2 µs (4.9× n, 4.5× b) | 151 µs (1.3× n, 1.2× b) | 192 µs (0.91× b) | 175 µs |
| 1e7 | 7.09 ms (1.3× n, 1.3× b) | 7.32 ms (1.3× n, 1.3× b) | 9.43 ms (1.0× b) | 9.41 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.48 µs (1.3× n, 1.2× b) | 2.44 µs (1.3× n, 1.2× b) | 3.14 µs (0.92× b) | 2.88 µs |
| 1e5 | 105 µs (1.7× n, 4.0× b) | 115 µs (1.6× n, 3.7× b) | 181 µs (2.3× b) | 426 µs |
| 1e7 | 12.5 ms (1.5× n, 3.8× b) | 12.8 ms (1.4× n, 3.7× b) | 18.2 ms (2.6× b) | 47.6 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.44 µs (1.8× n, 0.65× b) | 1.45 µs (1.8× n, 0.64× b) | 2.59 µs (0.36× b) | 0.932 µs |
| 1e5 | 10.9 µs (11.5× n, 5.7× b) | 16.4 µs (7.6× n, 3.8× b) | 126 µs (0.49× b) | 62.1 µs |
| 1e7 | 2.47 ms (5.1× n, 2.5× b) | 3.59 ms (3.5× n, 1.8× b) | 12.5 ms (0.50× b) | 6.29 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.75 µs (1.3× n, 1.0× b) | 2.24 µs (1.0× n, 0.78× b) | 2.32 µs (0.75× b) | 1.75 µs |
| 1e5 | 17.3 µs (4.1× n, 6.6× b) | 26.4 µs (2.7× n, 4.3× b) | 71 µs (1.6× b) | 113 µs |
| 1e7 | 3.31 ms (2.1× n, 5.0× b) | 4.04 ms (1.7× n, 4.1× b) | 7.07 ms (2.4× b) | 16.7 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.4 µs (1.6× n, 1.3× b) | 1.65 µs (1.4× n, 1.1× b) | 2.26 µs (0.84× b) | 1.89 µs |
| 1e5 | 15.7 µs (4.5× n, 7.1× b) | 25.8 µs (2.7× n, 4.3× b) | 71 µs (1.6× b) | 111 µs |
| 1e7 | 3.66 ms (2.2× n, 4.5× b) | 4.57 ms (1.7× n, 3.6× b) | 7.92 ms (2.1× b) | 16.4 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 4.22 µs (2.5× n, 0.28× b) | 6 µs (1.8× n, 0.20× b) | 10.6 µs (0.11× b) | 1.19 µs |
| 1e5 | 46.2 µs (41.4× n, 2.0× b) | 294 µs (6.5× n, 0.32× b) | 1.91 ms (0.05× b) | 94.1 µs |
| 1e7 | 7.7 ms (25.4× n, 1.6× b) | 31.9 ms (6.1× n, 0.39× b) | 195 ms (0.06× b) | 12.5 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 3.19 µs (1.3× n, 0.78× b) | 3.1 µs (1.3× n, 0.80× b) | 4 µs (0.62× b) | 2.47 µs |
| 1e5 | 29 µs (4.4× n, 6.1× b) | 35.7 µs (3.6× n, 5.0× b) | 128 µs (1.4× b) | 178 µs |
| 1e7 | 5.24 ms (2.4× n, 4.0× b) | 5.53 ms (2.3× n, 3.8× b) | 12.7 ms (1.7× b) | 21 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 3.01 µs (1.2× n, 0.82× b) | 3.1 µs (1.2× n, 0.79× b) | 3.7 µs (0.66× b) | 2.46 µs |
| 1e5 | 18.1 µs (4.6× n, 17.2× b) | 27.4 µs (3.0× n, 11.4× b) | 83.4 µs (3.7× b) | 312 µs |
| 1e7 | 2.53 ms (3.2× n, 8.4× b) | 2.85 ms (2.8× n, 7.5× b) | 8.06 ms (2.6× b) | 21.3 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.78 µs (2.9× n, 2.0× b) | 2.44 µs (2.1× n, 1.5× b) | 5.11 µs (0.70× b) | 3.6 µs |
| 1e5 | 148 µs (6.6× n, 2.7× b) | 239 µs (4.1× n, 1.7× b) | 973 µs (0.42× b) | 406 µs |
| 1e7 | 7.47 ms (11.7× n, 5.2× b) | 13.7 ms (6.4× n, 2.8× b) | 87.5 ms (0.44× b) | 38.9 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.85 µs (2.5× n, 14.5× b) | 2.81 µs (1.6× n, 9.6× b) | 4.58 µs (5.9× b) | 26.8 µs |
| 1e5 | 43.4 µs (19.3× n, 60.9× b) | 135 µs (6.2× n, 19.5× b) | 839 µs (3.2× b) | 2.64 ms |
| 1e7 | 7.86 ms (11.2× n, 36.9× b) | 16.4 ms (5.4× n, 17.7× b) | 88.1 ms (3.3× b) | 290 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 5.5 µs (0.54× n, 0.30× b) | 2.8 µs (1.1× n, 0.58× b) | 2.98 µs (0.55× b) | 1.62 µs |
| 1e5 | 72.8 µs (5.9× n, 5.8× b) | 122 µs (3.5× n, 3.5× b) | 429 µs (1.0× b) | 422 µs |
| 1e7 | 12.1 ms (1.6× n, 1.6× b) | 17 ms (1.1× n, 1.1× b) | 19.2 ms (1.0× b) | 19.1 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.56 µs (2.8× n, 5.1× b) | 3.6 µs (2.0× n, 3.6× b) | 7.28 µs (1.8× b) | 12.9 µs |
| 1e5 | 393 µs (1.8× n, 3.2× b) | 495 µs (1.4× n, 2.5× b) | 711 µs (1.8× b) | 1.25 ms |
| 1e7 | 15.5 ms (3.0× n, 6.5× b) | 25.9 ms (1.8× n, 3.9× b) | 47.2 ms (2.1× b) | 101 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.66 µs (2.1× n, 3.7× b) | 3.22 µs (1.8× n, 3.0× b) | 5.69 µs (1.7× b) | 9.81 µs |
| 1e5 | 397 µs (1.8× n, 3.2× b) | 457 µs (1.6× n, 2.7× b) | 724 µs (1.7× b) | 1.25 ms |
| 1e7 | 15.3 ms (3.1× n, 6.6× b) | 21.3 ms (2.2× n, 4.7× b) | 47.3 ms (2.1× b) | 101 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 24.2 µs (3.3× n, 0.81× b) | 61.4 µs (1.3× n, 0.32× b) | 80 µs (0.24× b) | 19.5 µs |
| 1e5 | 1.87 ms (4.0× n, 1.2× b) | 5.37 ms (1.4× n, 0.43× b) | 7.51 ms (0.31× b) | 2.33 ms |
| 1e7 | 161 ms (4.5× n, 1.3× b) | 529 ms (1.4× n, 0.40× b) | 721 ms (0.29× b) | 211 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 6.44 µs (2.1× n, 0.32× b) | 6.48 µs (2.1× n, 0.31× b) | 13.8 µs (0.15× b) | 2.03 µs |
| 1e5 | 697 µs (1.9× n, 1.1× b) | 704 µs (1.9× n, 1.1× b) | 1.35 ms (0.56× b) | 755 µs |
| 1e7 | 21.7 ms (3.9× n, 1.3× b) | 22.3 ms (3.8× n, 1.2× b) | 84.9 ms (0.32× b) | 27.2 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 8.24 µs (2.5× n, 1.5× b) | 10.9 µs (1.9× n, 1.1× b) | 20.4 µs (0.59× b) | 12.1 µs |
| 1e5 | 891 µs (1.6× n, 1.3× b) | 1.14 ms (1.3× n, 1.0× b) | 1.47 ms (0.82× b) | 1.2 ms |
| 1e7 | 38.5 ms (2.5× n, 1.8× b) | 63.6 ms (1.5× n, 1.1× b) | 96.6 ms (0.73× b) | 70.5 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.5 µs (6.3× n, 2.1× b) | 3.94 µs (4.0× n, 1.3× b) | 15.6 µs (0.33× b) | 5.2 µs |
| 1e5 | 59 µs (24.0× n, 8.5× b) | 235 µs (6.0× n, 2.1× b) | 1.42 ms (0.35× b) | 500 µs |
| 1e7 | 8.31 ms (17.0× n, 6.0× b) | 23.6 ms (6.0× n, 2.1× b) | 141 ms (0.35× b) | 49.9 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 5.84 µs (1.7× n, 1.1× b) | 20.3 µs (0.48× n, 0.32× b) | 9.76 µs (0.66× b) | 6.48 µs |
| 1e5 | 747 µs (3.1× n, 2.8× b) | 2.21 ms (1.1× n, 1.0× b) ⚠ | 2.32 ms (0.91× b) | 2.12 ms |
| 1e7 | 50.2 ms (4.2× n, 3.8× b) | 197 ms (1.1× n, 1.0× b) ⚠ | 211 ms (0.90× b) | 190 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 5.28 µs (5.4× n, 2.3× b) | 6.53 µs (4.3× n, 1.9× b) | 28.3 µs (0.43× b) | 12.2 µs |
| 1e5 | 736 µs (4.2× n, 1.6× b) | 851 µs (3.6× n, 1.4× b) | 3.08 ms (0.38× b) | 1.18 ms |
| 1e7 | 24 ms (10.7× n, 2.9× b) | 36.2 ms (7.1× n, 1.9× b) | 257 ms (0.27× b) | 69.1 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.62 µs (2.5× n, 1.3× b) | 1.75 µs (2.3× n, 1.2× b) | 4.04 µs (0.54× b) | 2.17 µs |
| 1e5 | 24.8 µs (10.1× n, 7.6× b) | 28.8 µs (8.7× n, 6.5× b) | 250 µs (0.75× b) | 187 µs |
| 1e7 | 3.66 ms (6.9× n, 5.2× b) | 3.96 ms (6.3× n, 4.8× b) | 25.1 ms (0.75× b) | 18.9 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.62 µs (1.5× n, 0.59× b) | 1.78 µs (1.4× n, 0.54× b) | 2.42 µs (0.40× b) | 0.962 µs |
| 1e5 | 20.3 µs (4.7× n, 3.3× b) | 21.2 µs (4.5× n, 3.1× b) | 95.6 µs (0.70× b) | 66.8 µs |
| 1e7 | 3.2 ms (3.0× n, 2.1× b) | 3.68 ms (2.6× n, 1.9× b) | 9.53 ms (0.72× b) | 6.87 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.48 µs (1.5× n, 0.60× b) | 1.71 µs (1.3× n, 0.52× b) | 2.25 µs (0.40× b) | 0.892 µs |
| 1e5 | 9.53 µs (8.3× n, 6.5× b) | 20.6 µs (3.8× n, 3.0× b) | 79.4 µs (0.78× b) | 62.1 µs |
| 1e7 | 1.71 ms (4.6× n, 3.7× b) | 2.22 ms (3.5× n, 2.8× b) | 7.86 ms (0.80× b) | 6.29 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.11 µs (2.0× n, 1.9× b) | 1.18 µs (1.9× n, 1.8× b) | 2.24 µs (0.94× b) | 2.11 µs |
| 1e5 | 17.4 µs (7.2× n, 9.5× b) | 18.9 µs (6.6× n, 8.7× b) | 126 µs (1.3× b) | 166 µs |
| 1e7 | 2.93 ms (4.3× n, 6.4× b) | 3.76 ms (3.3× n, 5.0× b) | 12.6 ms (1.5× b) | 18.8 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 2.58 µs (1.4× n, 0.65× b) | 2.52 µs (1.5× n, 0.67× b) | 3.71 µs (0.45× b) | 1.68 µs |
| 1e5 | 9.95 µs (12.8× n, 14.6× b) | 9.42 µs (13.5× n, 15.4× b) | 127 µs (1.1× b) | 145 µs |
| 1e7 | 1.81 ms (7.5× n, 8.7× b) | 1.81 ms (7.6× n, 8.7× b) | 13.6 ms (1.2× b) | 15.7 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.66 µs (2.1× n, 10.5× b) | 1.75 µs (2.0× n, 10.0× b) | 3.42 µs (5.1× b) | 17.5 µs |
| 1e5 | 6.62 µs (21.2× n, 250.8× b) | 11.7 µs (12.0× n, 141.9× b) | 141 µs (11.8× b) | 1.66 ms |
| 1e7 | 515 µs (27.3× n, 361.0× b) | 999 µs (14.1× n, 186.2× b) | 14.1 ms (13.2× b) | 186 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.56 µs (2.6× n, 1.4× b) | 1.68 µs (2.4× n, 1.3× b) | 4.06 µs (0.54× b) | 2.17 µs |
| 1e5 | 1.93 µs (6.2× n, 96.2× b) | 2.37 µs (5.0× n, 78.3× b) | 11.9 µs (15.6× b) | 186 µs |
| 1e7 | 1.86 µs (6.3× n, 10084.1× b) | 2.35 µs (5.0× n, 7980.6× b) | 11.8 µs (1589.1× b) | 18.8 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.45 µs (1.0× n, 0.26× b) | 1.5 µs (1.0× n, 0.25× b) | 1.43 µs (0.26× b) | 0.376 µs |
| 1e5 | 1.95 µs (0.67× n, 0.15× b) | 1.95 µs (0.67× n, 0.15× b) | 1.3 µs (0.22× b) | 0.291 µs |
| 1e7 | 1.88 µs (0.70× n, 0.17× b) | 1.91 µs (0.69× n, 0.17× b) | 1.31 µs (0.24× b) | 0.321 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.01 µs (1.0× n, 0.15× b) | 1.03 µs (1.0× n, 0.15× b) | 1.06 µs (0.14× b) | 0.151 µs |
| 1e5 | 0.932 µs (1.0× n, 0.14× b) | 0.952 µs (1.0× n, 0.14× b) | 0.972 µs (0.13× b) | 0.131 µs |
| 1e7 | 0.932 µs (1.1× n, 0.21× b) | 0.922 µs (1.1× n, 0.22× b) | 0.982 µs (0.20× b) | 0.2 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 0.992 µs (1.0× n, 0.12× b) | 1.01 µs (1.0× n, 0.12× b) | 1.02 µs (0.12× b) | 0.121 µs |
| 1e5 | 0.931 µs (1.1× n, 0.12× b) | 0.942 µs (1.0× n, 0.12× b) | 0.982 µs (0.11× b) | 0.111 µs |
| 1e7 | 0.921 µs (1.1× n, 0.14× b) | 0.922 µs (1.1× n, 0.14× b) | 0.982 µs (0.13× b) | 0.13 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.81 µs (2.3× n, 1.2× b) | 1.79 µs (2.3× n, 1.2× b) | 4.08 µs (0.53× b) | 2.18 µs |
| 1e5 | 19 µs (13.7× n, 9.8× b) | 28.2 µs (9.2× n, 6.6× b) | 259 µs (0.72× b) | 186 µs |
| 1e7 | 3.64 ms (7.1× n, 5.1× b) | 3.96 ms (6.6× n, 4.7× b) | 26 ms (0.72× b) | 18.6 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.68 µs (1.5× n, 0.83× b) | 1.73 µs (1.5× n, 0.81× b) | 2.59 µs (0.54× b) | 1.4 µs |
| 1e5 | 1.91 µs (3.0× n, 54.8× b) | 2.22 µs (2.6× n, 47.2× b) | 5.71 µs (18.4× b) | 105 µs |
| 1e7 | 1.82 µs (3.1× n, 5812.5× b) | 2.21 µs (2.5× n, 4808.4× b) | 5.58 µs (1900.1× b) | 10.6 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | avx2 | sse2 | none | base |
|---|---|---|---|---|
| 1e3 | 1.55 µs (1.5× n, 0.26× b) | 1.66 µs (1.4× n, 0.25× b) | 2.28 µs (0.18× b) | 0.411 µs |
| 1e5 | 1.65 µs (2.9× n, 0.18× b) | 2.2 µs (2.2× n, 0.13× b) | 4.83 µs (0.06× b) | 0.291 µs |
| 1e7 | 1.63 µs (2.9× n, 0.21× b) | 2.18 µs (2.2× n, 0.16× b) | 4.8 µs (0.07× b) | 0.341 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| fast | 3.22 ms (7.8× n, 5.8× b) | 3.83 ms (6.6× n, 4.9× b) | 25.1 ms (0.75× b) | 18.8 ms |
| pairwise | 3.49 ms (7.8× n, 5.4× b) | 3.93 ms (6.9× n, 4.8× b) | 27 ms (0.69× b) | 18.8 ms |
| compensated | 4.04 ms (3.5× n, 4.6× b) | 7.09 ms (2.0× n, 2.6× b) | 14.1 ms (1.3× b) | 18.8 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | avx2 | sse2 | none | base | blas |
|---|---|---|---|---|---|
| fast | 4.89 ms (4.0× n, 6.1× b) | 5.3 ms (3.7× n, 5.7× b) | 19.7 ms (1.5× b) | 30 ms | 12.7 ms (2.4× b) |
| pairwise | 4.92 ms (4.4× n, 6.1× b) | 5.39 ms (4.0× n, 5.6× b) | 21.6 ms (1.4× b) | 30.1 ms | 12.6 ms (2.4× b) |
| compensated | 5.39 ms (3.5× n, 5.6× b) | 9.02 ms (2.1× n, 3.3× b) | 18.9 ms (1.6× b) | 30.2 ms | 12.9 ms (2.3× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 56.2 ms (3.6× n, 3.8× b) | 157 ms (1.3× n, 1.4× b) | 205 ms (1.0× b) | 214 ms |
| fast | 36.2 ms (5.7× n, 5.9× b) | 80.7 ms (2.5× n, 2.6× b) | 205 ms (1.0× b) | 214 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 42.3 ms (1.3× n, 1.4× b) | 113 ms (0.49× n, 0.51× b) | 55.2 ms (1.1× b) | 58 ms |
| fast | 27.2 ms (2.0× n, 2.1× b) | 60.1 ms (0.92× n, 1.0× b) | 55.2 ms (1.1× b) | 58 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 71 ms (2.8× n, 2.8× b) | 222 ms (0.89× n, 0.89× b) | 196 ms (1.0× b) | 197 ms |
| fast | 33.2 ms (5.9× n, 5.9× b) | 73.4 ms (2.7× n, 2.7× b) | 196 ms (1.0× b) | 197 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 74.3 ms (4.2× n, 4.3× b) | 231 ms (1.4× n, 1.4× b) | 314 ms (1.0× b) | 320 ms |
| fast | 35.6 ms (8.8× n, 9.0× b) | 100 ms (3.1× n, 3.2× b) | 314 ms (1.0× b) | 320 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 48.7 ms (5.0× n, 2.2× b) | 194 ms (1.3× n, 0.56× b) | 244 ms (0.45× b) | 109 ms |
| fast | 17.4 ms (14.0× n, 6.3× b) | 31.5 ms (7.7× n, 3.5× b) | 244 ms (0.45× b) | 109 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 139 ms (1.8× n, 1.5× b) | 410 ms (0.59× n, 0.52× b) | 244 ms (0.87× b) | 212 ms |
| fast | 139 ms (1.8× n, 1.5× b) | 410 ms (0.59× n, 0.52× b) | 244 ms (0.87× b) | 212 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 131 ms (0.89× n, 0.92× b) | 370 ms (0.31× n, 0.33× b) | 116 ms (1.0× b) | 121 ms |
| fast | 131 ms (0.88× n, 0.92× b) | 369 ms (0.31× n, 0.33× b) | 116 ms (1.0× b) | 121 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 92.1 ms (3.8× n, 3.7× b) | 218 ms (1.6× n, 1.6× b) | 352 ms (1.0× b) | 345 ms |
| fast | 44.9 ms (7.9× n, 7.7× b) | 82.2 ms (4.3× n, 4.2× b) | 353 ms (1.0× b) | 345 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 82.5 ms (4.4× n, 4.3× b) | 197 ms (1.8× n, 1.8× b) | 359 ms (1.0× b) | 353 ms |
| fast | 65 ms (5.5× n, 5.4× b) | 157 ms (2.3× n, 2.2× b) | 359 ms (1.0× b) | 353 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 225 ms (2.6× n, 2.6× b) | 607 ms (1.0× n, 0.95× b) | 589 ms (1.0× b) | 577 ms |
| fast | 118 ms (5.0× n, 4.9× b) | 263 ms (2.2× n, 2.2× b) | 588 ms (1.0× b) | 576 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 124 ms (5.1× n, 5.0× b) | 321 ms (2.0× n, 1.9× b) | 628 ms (1.0× b) | 614 ms |
| fast | 104 ms (6.0× n, 5.9× b) | 285 ms (2.2× n, 2.2× b) | 627 ms (1.0× b) | 614 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 261 ms (3.8× n, 3.8× b) | 670 ms (1.5× n, 1.5× b) | 988 ms (1.0× b) | 989 ms |
| fast | 198 ms (5.0× n, 5.0× b) | 425 ms (2.3× n, 2.3× b) | 988 ms (1.0× b) | 989 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 127 ms (1.9× n, 1.8× b) | 251 ms (0.94× n, 0.90× b) | 235 ms (1.0× b) | 227 ms |
| fast | 103 ms (2.3× n, 2.2× b) | 173 ms (1.4× n, 1.3× b) | 235 ms (1.0× b) | 227 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | avx2 | sse2 | none | base |
|---|---|---|---|---|
| accurate | 329 ms (2.8× n, 2.7× b) | 846 ms (1.1× n, 1.1× b) | 916 ms (1.0× b) | 899 ms |
| fast | 190 ms (4.8× n, 4.7× b) | 450 ms (2.0× n, 2.0× b) | 916 ms (1.0× b) | 900 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.251 µs | 1.24 µs (0.20× b) | 0.851 µs (0.30× b) |
| 1e1 | 0.261 µs | 1.25 µs (0.21× b) | 0.842 µs (0.31× b) |
| 1e2 | 0.421 µs | 1.26 µs (0.33× b) | 0.851 µs (0.49× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.1 µs | 1.5 µs (0.07× b) |
| 1e1 | 0.391 µs | 1.57 µs (0.25× b) |
| 1e2 | 0.231 µs | 1.66 µs (0.14× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.11 µs | 1.14 µs (0.10× b) |
| 1e1 | 0.17 µs | 1.19 µs (0.14× b) |
| 1e2 | 0.711 µs | 1.4 µs (0.51× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.0811 µs | 1.55 µs (0.05× b) |
| 1e1 | 0.171 µs | 1.55 µs (0.11× b) |
| 1e2 | 0.301 µs | 1.77 µs (0.17× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.1 µs | 1.7 µs (0.06× b) |
| 1e1 | 0.11 µs | 1.69 µs (0.07× b) |
| 1e2 | 0.24 µs | 1.77 µs (0.14× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.34 µs | 1.73 µs (0.20× b) |
| 1e1 | 0.41 µs | 1.72 µs (0.24× b) |
| 1e2 | 0.651 µs | 1.77 µs (0.37× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.13 µs | 4.36 µs (0.03× b) |
| 1e1 | 0.201 µs | 4.44 µs (0.05× b) |
| 1e2 | 0.361 µs | 4.63 µs (0.08× b) |

