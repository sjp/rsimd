# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-11T22:09:47Z |
| rsimd | 0.1.0 (git d95c7863d200e98e30b6f9679a28d91f9d7d1615) |
| rsimd library | `/home/runner/work/_temp/Library/rsimd` (built R 4.6.1; x86_64-pc-linux-gnu; 2026-10-11 21:53:52 UTC; unix) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `gcc -std=gnu2x` |
| CFLAGS | `-g -O2` |
| OS | Linux 6.17.0-1022-azure (x86_64) |
| CPU | AMD EPYC 9V74 80-Core Processor, 4 cores |
| CPU features | sse2, sse3, ssse3, sse4_1, sse4_2, avx, avx2, fma, avx512f, avx512bw, avx512dq, avx512vl, os_avx, os_avx512 |
| Tiers | avx512, avx2, sse2, none (auto: avx512) |
| bench | 1.1.4 |
| Run time | 803 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.38 µs (2.3× n, 1.4× b) | 1.41 µs (2.2× n, 1.4× b) | 1.42 µs (2.2× n, 1.4× b) | 3.15 µs (0.62× b) | 1.96 µs |
| 1e5 | 8.7 µs (22.0× n, 20.9× b) | 8.74 µs (21.9× n, 20.8× b) | 12.2 µs (15.7× n, 14.8× b) | 192 µs (0.95× b) | 182 µs |
| 1e7 | 2.32 ms (8.3× n, 7.9× b) | 1.94 ms (9.9× n, 9.4× b) | 2.11 ms (9.1× n, 8.7× b) | 19.2 ms (1.0× b) | 18.3 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.24 µs (1.4× n, 0.59× b) | 1.23 µs (1.4× n, 0.59× b) | 1.3 µs (1.3× n, 0.56× b) | 1.73 µs (0.42× b) | 0.731 µs |
| 1e5 | 8.31 µs (6.7× n, 6.5× b) | 6.99 µs (7.9× n, 7.8× b) | 12.9 µs (4.3× n, 4.2× b) | 55.6 µs (1.0× b) | 54.4 µs |
| 1e7 | 924 µs (5.9× n, 5.9× b) | 944 µs (5.8× n, 5.8× b) | 1.25 ms (4.4× n, 4.4× b) | 5.5 ms (1.0× b) | 5.49 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.36 µs (2.4× n, 4.5× b) | 1.38 µs (2.3× n, 4.4× b) | 1.41 µs (2.3× n, 4.3× b) | 3.23 µs (1.9× b) | 6.11 µs |
| 1e5 | 8.73 µs (22.0× n, 43.9× b) | 8.7 µs (22.0× n, 44.0× b) | 12.2 µs (15.7× n, 31.3× b) | 192 µs (2.0× b) | 383 µs |
| 1e7 | 2.16 ms (8.9× n, 17.5× b) | 1.88 ms (10.2× n, 20.1× b) | 2.01 ms (9.6× n, 18.9× b) | 19.2 ms (2.0× b) | 37.8 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | avx512 | avx2 | sse2 | none | base | blas |
|---|---|---|---|---|---|---|
| 1e3 | 1.78 µs (2.5× n, 3.1× b) | 1.75 µs (2.6× n, 3.2× b) | 1.85 µs (2.4× n, 3.0× b) | 4.52 µs (1.2× b) | 5.56 µs | 0.982 µs (5.7× b) |
| 1e5 | 19.1 µs (7.4× n, 25.7× b) | 19 µs (7.4× n, 25.9× b) | 19.1 µs (7.4× n, 25.7× b) | 141 µs (3.5× b) | 491 µs | 58.4 µs (8.4× b) |
| 1e7 | 4.56 ms (3.1× n, 6.1× b) | 4.48 ms (3.2× n, 6.2× b) | 4.58 ms (3.1× n, 6.1× b) | 14.2 ms (1.9× b) | 27.8 ms | 8.84 ms (3.1× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.82 µs (0.56× n, 0.73× b) | 2.46 µs (1.1× n, 1.4× b) | 4.97 µs (0.55× n, 0.71× b) | 2.71 µs (1.3× b) | 3.53 µs |
| 1e5 | 319 µs (1.1× n, 1.0× b) | 299 µs (1.2× n, 1.1× b) | 312 µs (1.1× n, 1.0× b) | 350 µs (0.91× b) | 320 µs |
| 1e7 | 8.55 ms (1.4× n, 1.1× b) | 8.76 ms (1.3× n, 1.1× b) | 9.08 ms (1.3× n, 1.1× b) | 11.7 ms (0.83× b) | 9.71 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.12 µs (1.4× n, 0.65× b) | 2.1 µs (1.4× n, 0.66× b) | 2.2 µs (1.3× n, 0.63× b) | 2.93 µs (0.47× b) | 1.38 µs |
| 1e5 | 157 µs (1.6× n, 3.5× b) | 153 µs (1.7× n, 3.6× b) | 178 µs (1.4× n, 3.1× b) | 254 µs (2.2× b) | 553 µs |
| 1e7 | 4.95 ms (2.8× n, 9.3× b) | 5.05 ms (2.7× n, 9.1× b) | 6.34 ms (2.2× n, 7.3× b) | 13.8 ms (3.3× b) | 46.1 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.51 µs (1.5× n, 0.68× b) | 2.47 µs (1.6× n, 0.69× b) | 3.82 µs (1.0× n, 0.45× b) | 3.84 µs (0.45× b) | 1.71 µs |
| 1e5 | 357 µs (1.3× n, 1.0× b) | 316 µs (1.5× n, 1.2× b) | 468 µs (1.0× n, 0.78× b) | 463 µs (0.79× b) | 365 µs |
| 1e7 | 10.3 ms (2.2× n, 1.5× b) | 10.5 ms (2.2× n, 1.4× b) | 22.7 ms (1.0× n, 0.66× b) | 22.7 ms (0.66× b) | 15 ms |

### scalar_sub (double)

`simd_sub(1, x)` vs base `1 - x`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.27 µs (2.3× n, 0.40× b) | 2.33 µs (2.2× n, 0.39× b) | 2.35 µs (2.2× n, 0.39× b) | 5.23 µs (0.17× b) | 0.911 µs |
| 1e5 | 373 µs (1.7× n, 0.78× b) | 332 µs (1.9× n, 0.88× b) | 308 µs (2.0× n, 0.95× b) | 620 µs (0.47× b) | 292 µs |
| 1e7 | 7.05 ms (5.5× n, 1.0× b) | 7.14 ms (5.5× n, 1.0× b) | 7.28 ms (5.3× n, 1.0× b) | 38.9 ms (0.18× b) | 7.09 ms |

### fma_scalar (double)

`simd_fma(x, 2, 1)` vs base `x * 2 + 1`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.47 µs (1.7× n, 0.63× b) | 2.46 µs (1.8× n, 0.63× b) | 3.79 µs (1.1× n, 0.41× b) | 4.32 µs (0.36× b) | 1.55 µs |
| 1e5 | 363 µs (1.3× n, 1.0× b) | 323 µs (1.5× n, 1.1× b) | 456 µs (1.1× n, 0.77× b) | 486 µs (0.72× b) | 350 µs |
| 1e7 | 7 ms (3.6× n, 1.8× b) | 7.04 ms (3.6× n, 1.8× b) | 22.6 ms (1.1× n, 0.56× b) | 25.4 ms (0.50× b) | 12.6 ms |

### lt_scalar (double)

`simd_lt(x, 0)` vs base `x < 0`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.96 µs (1.5× n, 0.59× b) | 1.96 µs (1.5× n, 0.59× b) | 2.05 µs (1.4× n, 0.57× b) | 2.85 µs (0.41× b) | 1.16 µs |
| 1e5 | 160 µs (1.6× n, 1.6× b) | 157 µs (1.7× n, 1.6× b) | 174 µs (1.5× n, 1.4× b) | 259 µs (1.0× b) | 248 µs |
| 1e7 | 5.44 ms (2.7× n, 2.2× b) | 5.48 ms (2.7× n, 2.2× b) | 6.22 ms (2.4× n, 2.0× b) | 14.7 ms (0.83× b) | 12.2 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.62 µs (1.1× n, 2.6× b) | 2.61 µs (1.1× n, 2.6× b) | 2.66 µs (1.0× n, 2.6× b) | 2.77 µs (2.5× b) | 6.83 µs |
| 1e5 | 357 µs (0.91× n, 2.4× b) ⚠ | 306 µs (1.1× n, 2.8× b) ⚠ | 295 µs (1.1× n, 2.9× b) | 326 µs (2.6× b) | 848 µs |
| 1e7 | 8.14 ms (1.2× n, 7.7× b) | 8.47 ms (1.2× n, 7.4× b) | 9.05 ms (1.1× n, 6.9× b) | 10.2 ms (6.1× b) | 62.5 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.22 µs (2.6× n, 2.6× b) | 2.67 µs (2.1× n, 2.1× b) | 4.92 µs (1.2× n, 1.2× b) | 5.72 µs (1.0× b) | 5.68 µs |
| 1e5 | 349 µs (2.0× n, 2.0× b) | 397 µs (1.8× n, 1.8× b) | 622 µs (1.1× n, 1.1× b) | 702 µs (1.0× b) | 706 µs |
| 1e7 | 12.1 ms (3.9× n, 3.9× b) | 16.5 ms (2.8× n, 2.9× b) | 39.3 ms (1.2× n, 1.2× b) | 47.1 ms (1.0× b) | 47.3 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.902 µs (1.2× n, 0.38× b) | 0.921 µs (1.2× n, 0.37× b) | 0.941 µs (1.2× n, 0.36× b) | 1.12 µs (0.30× b) | 0.34 µs |
| 1e5 | 8.34 µs (3.4× n, 3.3× b) | 8.67 µs (3.3× n, 3.1× b) | 11.7 µs (2.4× n, 2.3× b) | 28.5 µs (1.0× b) | 27.2 µs |
| 1e7 | 1.97 ms (1.4× n, 1.4× b) | 2.03 ms (1.4× n, 1.4× b) | 2.04 ms (1.4× n, 1.4× b) | 2.81 ms (1.0× b) | 2.8 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.902 µs (1.2× n, 0.38× b) | 0.862 µs (1.2× n, 0.40× b) | 0.891 µs (1.2× n, 0.38× b) | 1.07 µs (0.32× b) | 0.341 µs |
| 1e5 | 4.72 µs (6.1× n, 5.8× b) | 4.31 µs (6.7× n, 6.3× b) | 5.63 µs (5.1× n, 4.8× b) | 28.9 µs (0.94× b) | 27.2 µs |
| 1e7 | 818 µs (3.4× n, 3.4× b) | 836 µs (3.3× n, 3.3× b) | 911 µs (3.1× n, 3.0× b) | 2.79 ms (1.0× b) | 2.76 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.42 µs (1.2× n, 0.58× b) | 1.36 µs (1.2× n, 0.61× b) | 1.38 µs (1.2× n, 0.60× b) | 1.64 µs (0.51× b) | 0.831 µs |
| 1e5 | 159 µs (1.2× n, 1.2× b) | 146 µs (1.3× n, 1.3× b) | 149 µs (1.2× n, 1.2× b) | 185 µs (1.0× b) | 185 µs |
| 1e7 | 4.73 ms (1.6× n, 1.6× b) | 5.25 ms (1.4× n, 1.4× b) | 5.44 ms (1.4× n, 1.4× b) | 7.59 ms (1.0× b) | 7.55 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.98 µs (1.2× n, 1.1× b) | 1.95 µs (1.3× n, 1.1× b) | 2 µs (1.2× n, 1.1× b) | 2.45 µs (0.88× b) | 2.16 µs |
| 1e5 | 144 µs (1.4× n, 3.4× b) | 143 µs (1.4× n, 3.4× b) | 150 µs (1.3× n, 3.2× b) | 196 µs (2.5× b) | 485 µs |
| 1e7 | 9.81 ms (1.5× n, 4.1× b) | 9.74 ms (1.5× n, 4.2× b) | 10.2 ms (1.4× n, 4.0× b) | 14.6 ms (2.8× b) | 40.6 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.24 µs (1.5× n, 0.58× b) | 1.24 µs (1.5× n, 0.58× b) | 1.25 µs (1.5× n, 0.58× b) | 1.92 µs (0.37× b) | 0.721 µs |
| 1e5 | 8.23 µs (10.1× n, 6.6× b) | 8.01 µs (10.4× n, 6.8× b) | 12 µs (6.9× n, 4.5× b) | 83.4 µs (0.65× b) | 54.4 µs |
| 1e7 | 1.92 ms (4.3× n, 2.9× b) | 1.87 ms (4.4× n, 2.9× b) | 1.99 ms (4.1× n, 2.7× b) | 8.22 ms (0.67× b) | 5.48 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.53 µs (1.4× n, 1.2× b) | 1.54 µs (1.3× n, 1.2× b) | 1.64 µs (1.3× n, 1.2× b) | 2.07 µs (0.92× b) | 1.91 µs |
| 1e5 | 10.4 µs (6.0× n, 35.4× b) | 10.8 µs (5.8× n, 34.2× b) | 20.3 µs (3.0× n, 18.1× b) | 61.9 µs (5.9× b) | 368 µs |
| 1e7 | 1.96 ms (3.1× n, 7.8× b) | 2.09 ms (2.9× n, 7.3× b) | 2.2 ms (2.8× n, 7.0× b) | 6.06 ms (2.5× b) | 15.3 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.15 µs (1.7× n, 1.6× b) | 1.16 µs (1.7× n, 1.6× b) | 1.29 µs (1.5× n, 1.4× b) | 1.92 µs (1.0× b) | 1.86 µs |
| 1e5 | 9.84 µs (6.3× n, 35.9× b) | 10.2 µs (6.1× n, 34.8× b) | 19.7 µs (3.1× n, 17.9× b) | 61.6 µs (5.7× b) | 353 µs |
| 1e7 | 2.25 ms (3.0× n, 5.6× b) | 2.38 ms (2.9× n, 5.3× b) | 2.43 ms (2.8× n, 5.2× b) | 6.79 ms (1.9× b) | 12.6 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.45 µs (2.8× n, 0.31× b) | 3.56 µs (2.7× n, 0.30× b) | 5.61 µs (1.7× n, 0.19× b) | 9.62 µs (0.11× b) | 1.07 µs |
| 1e5 | 173 µs (10.7× n, 1.2× b) | 173 µs (10.7× n, 1.2× b) | 385 µs (4.8× n, 0.56× b) | 1.85 ms (0.12× b) | 215 µs |
| 1e7 | 5.73 ms (30.9× n, 1.8× b) | 5.93 ms (29.9× n, 1.8× b) | 26.7 ms (6.6× n, 0.40× b) | 177 ms (0.06× b) | 10.6 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.52 µs (1.3× n, 0.86× b) | 2.55 µs (1.3× n, 0.85× b) | 2.65 µs (1.3× n, 0.82× b) | 3.4 µs (0.64× b) | 2.18 µs |
| 1e5 | 20.7 µs (5.0× n, 14.8× b) | 20 µs (5.2× n, 15.3× b) | 31.3 µs (3.3× n, 9.8× b) | 104 µs (3.0× b) | 306 µs |
| 1e7 | 4.15 ms (2.4× n, 4.3× b) | 4.15 ms (2.4× n, 4.3× b) | 4.37 ms (2.3× n, 4.1× b) | 10.1 ms (1.8× b) | 18 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.47 µs (1.3× n, 0.87× b) | 2.48 µs (1.2× n, 0.87× b) | 2.59 µs (1.2× n, 0.83× b) | 3.1 µs (0.70× b) | 2.16 µs |
| 1e5 | 14.2 µs (5.2× n, 21.2× b) | 14.8 µs (4.9× n, 20.3× b) | 23.8 µs (3.1× n, 12.6× b) | 73.1 µs (4.1× b) | 300 µs |
| 1e7 | 1.9 ms (3.6× n, 9.7× b) | 1.88 ms (3.7× n, 9.8× b) | 2.27 ms (3.0× n, 8.1× b) | 6.91 ms (2.7× b) | 18.4 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.48 µs (2.8× n, 2.2× b) | 1.51 µs (2.7× n, 2.2× b) | 1.98 µs (2.1× n, 1.7× b) | 4.08 µs (0.81× b) | 3.3 µs |
| 1e5 | 159 µs (5.2× n, 1.7× b) | 153 µs (5.3× n, 1.8× b) | 214 µs (3.8× n, 1.3× b) | 818 µs (0.34× b) | 276 µs |
| 1e7 | 5.17 ms (14.4× n, 6.2× b) | 5.53 ms (13.4× n, 5.8× b) | 10.4 ms (7.1× n, 3.1× b) | 74.2 ms (0.43× b) | 32.2 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.52 µs (2.4× n, 13.8× b) | 1.55 µs (2.3× n, 13.6× b) | 2.22 µs (1.6× n, 9.5× b) | 3.62 µs (5.8× b) | 21.1 µs |
| 1e5 | 28.6 µs (30.8× n, 66.1× b) | 165 µs (5.3× n, 11.4× b) | 240 µs (3.7× n, 7.9× b) | 881 µs (2.1× b) | 1.89 ms |
| 1e7 | 5.56 ms (14.1× n, 38.5× b) | 5.86 ms (13.4× n, 36.5× b) | 12.6 ms (6.2× n, 17.0× b) | 78.5 ms (2.7× b) | 214 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.91 µs (1.4× n, 0.82× b) | 1.95 µs (1.3× n, 0.80× b) | 2.31 µs (1.1× n, 0.68× b) | 2.6 µs (0.60× b) | 1.57 µs |
| 1e5 | 33 µs (12.3× n, 12.1× b) | 319 µs (1.3× n, 1.3× b) | 361 µs (1.1× n, 1.1× b) | 404 µs (1.0× b) | 400 µs |
| 1e7 | 6.92 ms (2.4× n, 2.4× b) | 8.81 ms (1.9× n, 1.9× b) | 13.1 ms (1.3× n, 1.3× b) | 16.5 ms (1.0× b) | 16.9 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.08 µs (2.4× n, 4.2× b) | 2.18 µs (2.3× n, 4.0× b) | 3.1 µs (1.6× n, 2.8× b) | 4.95 µs (1.8× b) | 8.71 µs |
| 1e5 | 370 µs (1.7× n, 3.0× b) | 363 µs (1.7× n, 3.1× b) | 449 µs (1.4× n, 2.5× b) | 632 µs (1.8× b) | 1.11 ms |
| 1e7 | 11.9 ms (3.4× n, 7.4× b) | 12.6 ms (3.2× n, 7.0× b) | 21.9 ms (1.9× n, 4.0× b) | 40.5 ms (2.2× b) | 88 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.92 µs (2.6× n, 4.5× b) | 2.13 µs (2.3× n, 4.1× b) | 2.7 µs (1.8× n, 3.2× b) | 4.96 µs (1.8× b) | 8.74 µs |
| 1e5 | 299 µs (2.1× n, 3.7× b) | 346 µs (1.8× n, 3.2× b) | 400 µs (1.6× n, 2.8× b) | 636 µs (1.8× b) | 1.12 ms |
| 1e7 | 7 ms (5.8× n, 12.6× b) | 11.1 ms (3.7× n, 8.0× b) | 17.2 ms (2.4× n, 5.1× b) | 40.5 ms (2.2× b) | 88.1 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 12.8 µs (4.6× n, 1.1× b) | 19.1 µs (3.1× n, 0.73× b) | 50.5 µs (1.2× n, 0.27× b) | 58.7 µs (0.24× b) | 13.8 µs |
| 1e5 | 964 µs (5.7× n, 1.9× b) | 1.56 ms (3.5× n, 1.2× b) | 4.71 ms (1.2× n, 0.39× b) | 5.52 ms (0.33× b) | 1.82 ms |
| 1e7 | 70 ms (7.5× n, 2.3× b) | 131 ms (4.0× n, 1.2× b) | 446 ms (1.2× n, 0.36× b) | 527 ms (0.30× b) | 160 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.37 µs (1.1× n, 0.60× b) | 3.37 µs (1.1× n, 0.60× b) | 3.36 µs (1.1× n, 0.60× b) | 3.71 µs (0.54× b) | 2.01 µs |
| 1e5 | 665 µs (1.0× n, 1.1× b) | 611 µs (1.1× n, 1.2× b) | 634 µs (1.1× n, 1.2× b) | 693 µs (1.1× b) | 732 µs |
| 1e7 | 18.7 ms (1.1× n, 1.2× b) | 17.3 ms (1.2× n, 1.3× b) | 18.1 ms (1.2× n, 1.2× b) | 21.4 ms (1.0× b) | 22 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.92 µs (2.5× n, 1.0× b) | 4.16 µs (3.0× n, 1.1× b) | 6.49 µs (1.9× n, 0.72× b) | 12.4 µs (0.38× b) | 4.7 µs |
| 1e5 | 829 µs (1.9× n, 1.2× b) | 741 µs (2.2× n, 1.4× b) | 1 ms (1.6× n, 1.0× b) | 1.6 ms (0.64× b) | 1.01 ms |
| 1e7 | 34.5 ms (3.2× n, 1.5× b) | 26.8 ms (4.2× n, 1.9× b) | 52.1 ms (2.1× n, 1.0× b) | 112 ms (0.47× b) | 51.9 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.81 µs (4.7× n, 2.5× b) | 1.91 µs (4.5× n, 2.3× b) | 3.43 µs (2.5× n, 1.3× b) | 8.58 µs (0.52× b) | 4.49 µs |
| 1e5 | 42.8 µs (17.1× n, 10.1× b) | 51.2 µs (14.3× n, 8.5× b) | 207 µs (3.5× n, 2.1× b) | 733 µs (0.59× b) | 434 µs |
| 1e7 | 4.51 ms (16.2× n, 9.7× b) | 5.21 ms (14.0× n, 8.4× b) | 20.6 ms (3.5× n, 2.1× b) | 72.9 ms (0.60× b) | 43.7 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.46 µs (2.1× n, 1.5× b) | 4.89 µs (1.5× n, 1.0× b) | 17.9 µs (0.41× n, 0.28× b) | 7.28 µs (0.70× b) | 5.09 µs |
| 1e5 | 500 µs (4.1× n, 3.8× b) | 644 µs (3.2× n, 2.9× b) | 1.97 ms (1.0× n, 1.0× b) ⚠ | 2.06 ms (0.92× b) | 1.89 ms |
| 1e7 | 26.9 ms (6.8× n, 6.3× b) | 41.3 ms (4.4× n, 4.1× b) | 173 ms (1.1× n, 1.0× b) ⚠ | 182 ms (0.92× b) | 168 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 5.94 µs (1.4× n, 0.87× b) | 4.54 µs (1.8× n, 1.1× b) | 5.08 µs (1.6× n, 1.0× b) | 8.3 µs (0.62× b) | 5.14 µs |
| 1e5 | 813 µs (1.3× n, 1.3× b) | 651 µs (1.6× n, 1.6× b) | 710 µs (1.5× n, 1.5× b) | 1.06 ms (1.0× b) | 1.04 ms |
| 1e7 | 34.8 ms (1.7× n, 1.6× b) | 20.1 ms (2.9× n, 2.7× b) | 26.5 ms (2.2× n, 2.0× b) | 58.4 ms (0.93× b) | 54.2 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.34 µs (2.4× n, 1.4× b) | 1.36 µs (2.3× n, 1.4× b) | 1.48 µs (2.1× n, 1.2× b) | 3.16 µs (0.58× b) | 1.84 µs |
| 1e5 | 11.3 µs (16.9× n, 15.5× b) | 14.9 µs (12.9× n, 11.8× b) | 24 µs (8.0× n, 7.3× b) | 192 µs (0.92× b) | 176 µs |
| 1e7 | 2.01 ms (9.6× n, 8.4× b) | 2.09 ms (9.2× n, 8.1× b) | 2.51 ms (7.7× n, 6.8× b) | 19.2 ms (0.88× b) | 16.9 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.45 µs (1.5× n, 0.50× b) | 1.46 µs (1.5× n, 0.49× b) | 1.52 µs (1.4× n, 0.47× b) | 2.19 µs (0.33× b) | 0.721 µs |
| 1e5 | 10.3 µs (8.1× n, 5.3× b) | 10.7 µs (7.8× n, 5.1× b) | 17.1 µs (4.9× n, 3.2× b) | 83.6 µs (0.65× b) | 54.6 µs |
| 1e7 | 1.99 ms (4.1× n, 4.1× b) | 2.19 ms (3.8× n, 3.8× b) | 2.15 ms (3.8× n, 3.8× b) | 8.25 ms (1.0× b) | 8.22 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.39 µs (1.4× n, 0.52× b) | 1.39 µs (1.4× n, 0.52× b) | 1.51 µs (1.3× n, 0.48× b) | 1.97 µs (0.37× b) | 0.722 µs |
| 1e5 | 5.14 µs (11.9× n, 10.6× b) | 7.09 µs (8.6× n, 7.7× b) | 17.7 µs (3.5× n, 3.1× b) | 61.3 µs (0.89× b) | 54.4 µs |
| 1e7 | 829 µs (7.3× n, 6.6× b) | 911 µs (6.7× n, 6.0× b) | 1.68 ms (3.6× n, 3.3× b) | 6.07 ms (0.90× b) | 5.48 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.862 µs (1.9× n, 2.4× b) | 0.872 µs (1.9× n, 2.3× b) | 0.962 µs (1.7× n, 2.1× b) | 1.62 µs (1.3× b) | 2.04 µs |
| 1e5 | 8.53 µs (9.7× n, 31.9× b) | 8.77 µs (9.4× n, 31.0× b) | 16.4 µs (5.0× n, 16.6× b) | 82.5 µs (3.3× b) | 272 µs |
| 1e7 | 1.99 ms (4.1× n, 8.0× b) | 2.06 ms (4.0× n, 7.7× b) | 2.26 ms (3.6× n, 7.0× b) | 8.24 ms (1.9× b) | 15.9 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.73 µs (1.6× n, 0.77× b) | 1.69 µs (1.6× n, 0.79× b) | 1.68 µs (1.6× n, 0.80× b) | 2.69 µs (0.50× b) | 1.34 µs |
| 1e5 | 42.5 µs (3.4× n, 3.6× b) | 38.5 µs (3.8× n, 3.9× b) | 38.5 µs (3.8× n, 3.9× b) | 145 µs (1.0× b) | 151 µs |
| 1e7 | 1.53 ms (7.7× n, 8.3× b) | 1.52 ms (7.7× n, 8.3× b) | 1.57 ms (7.5× n, 8.1× b) | 11.7 ms (1.1× b) | 12.6 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.2 µs (2.0× n, 11.5× b) | 1.18 µs (2.0× n, 11.7× b) | 1.25 µs (1.9× n, 11.0× b) | 2.41 µs (5.7× b) | 13.8 µs |
| 1e5 | 5.25 µs (24.1× n, 239.8× b) | 5.48 µs (23.1× n, 229.8× b) | 9.83 µs (12.9× n, 128.0× b) | 127 µs (9.9× b) | 1.26 ms |
| 1e7 | 403 µs (31.3× n, 348.2× b) | 431 µs (29.3× n, 325.4× b) | 873 µs (14.4× n, 160.6× b) | 12.6 ms (11.1× b) | 140 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.31 µs (2.4× n, 1.5× b) | 1.3 µs (2.4× n, 1.5× b) | 1.37 µs (2.3× n, 1.4× b) | 3.19 µs (0.62× b) | 1.96 µs |
| 1e5 | 1.57 µs (5.9× n, 115.4× b) | 1.56 µs (6.0× n, 116.2× b) | 1.71 µs (5.4× n, 105.9× b) | 9.32 µs (19.5× b) | 181 µs |
| 1e7 | 1.5 µs (6.1× n, 12164.8× b) | 1.51 µs (6.1× n, 12091.8× b) | 1.67 µs (5.5× n, 10928.2× b) | 9.23 µs (1980.1× b) | 18.3 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.25 µs (1.0× n, 0.18× b) | 1.26 µs (1.0× n, 0.17× b) | 1.3 µs (0.93× n, 0.17× b) | 1.21 µs (0.18× b) | 0.22 µs |
| 1e5 | 2.39 µs (0.50× n, 0.08× b) | 1.73 µs (0.69× n, 0.12× b) | 1.71 µs (0.70× n, 0.12× b) | 1.19 µs (0.17× b) | 0.2 µs |
| 1e7 | 2.33 µs (0.50× n, 0.10× b) | 1.65 µs (0.71× n, 0.14× b) | 1.67 µs (0.70× n, 0.14× b) | 1.17 µs (0.20× b) | 0.23 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.832 µs (1.0× n, 0.12× b) | 0.831 µs (1.0× n, 0.12× b) | 0.822 µs (1.0× n, 0.12× b) | 0.821 µs (0.12× b) | 0.1 µs |
| 1e5 | 0.852 µs (1.0× n, 0.11× b) | 0.851 µs (1.0× n, 0.11× b) | 0.831 µs (1.0× n, 0.11× b) | 0.841 µs (0.11× b) | 0.09 µs |
| 1e7 | 0.801 µs (1.0× n, 0.12× b) | 0.821 µs (1.0× n, 0.12× b) | 0.821 µs (1.0× n, 0.12× b) | 0.821 µs (0.12× b) | 0.1 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.822 µs (1.0× n, 0.10× b) | 0.831 µs (1.0× n, 0.10× b) | 0.821 µs (1.0× n, 0.10× b) | 0.812 µs (0.10× b) | 0.0801 µs |
| 1e5 | 0.822 µs (1.0× n, 0.09× b) | 0.821 µs (1.0× n, 0.09× b) | 0.822 µs (1.0× n, 0.09× b) | 0.821 µs (0.09× b) | 0.0701 µs |
| 1e7 | 0.792 µs (1.0× n, 0.11× b) | 0.792 µs (1.0× n, 0.11× b) | 0.801 µs (1.0× n, 0.11× b) | 0.791 µs (0.11× b) | 0.09 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.38 µs (2.3× n, 1.4× b) | 1.4 µs (2.3× n, 1.3× b) | 1.54 µs (2.1× n, 1.2× b) | 3.22 µs (0.59× b) | 1.88 µs |
| 1e5 | 11.7 µs (17.0× n, 14.7× b) | 15.3 µs (12.9× n, 11.1× b) | 24.1 µs (8.2× n, 7.1× b) | 198 µs (0.86× b) | 171 µs |
| 1e7 | 2.13 ms (9.4× n, 8.1× b) | 2.15 ms (9.3× n, 8.0× b) | 2.48 ms (8.1× n, 6.9× b) | 20 ms (0.86× b) | 17.2 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.52 µs (1.5× n, 0.66× b) | 1.53 µs (1.5× n, 0.66× b) | 1.56 µs (1.4× n, 0.65× b) | 2.22 µs (0.46× b) | 1.01 µs |
| 1e5 | 1.93 µs (2.6× n, 42.7× b) | 1.9 µs (2.6× n, 43.3× b) | 2.03 µs (2.4× n, 40.6× b) | 4.96 µs (16.6× b) | 82.5 µs |
| 1e7 | 1.7 µs (2.9× n, 5045.1× b) | 1.7 µs (2.9× n, 5042.4× b) | 1.99 µs (2.5× n, 4308.5× b) | 4.9 µs (1753.6× b) | 8.59 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.37 µs (1.4× n, 0.15× b) | 1.38 µs (1.4× n, 0.15× b) | 1.5 µs (1.3× n, 0.14× b) | 1.95 µs (0.11× b) | 0.21 µs |
| 1e5 | 1.55 µs (2.5× n, 0.12× b) | 1.62 µs (2.4× n, 0.12× b) | 2.05 µs (1.9× n, 0.09× b) | 3.92 µs (0.05× b) | 0.19 µs |
| 1e7 | 1.46 µs (2.7× n, 0.14× b) | 1.52 µs (2.6× n, 0.14× b) | 2.02 µs (1.9× n, 0.10× b) | 3.92 µs (0.05× b) | 0.211 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| fast | 2.21 ms (8.7× n, 8.3× b) | 1.89 ms (10.2× n, 9.7× b) | 2.04 ms (9.4× n, 9.0× b) | 19.2 ms (1.0× b) | 18.3 ms |
| pairwise | 2.07 ms (9.7× n, 8.8× b) | 2.02 ms (10.0× n, 9.0× b) | 2.12 ms (9.5× n, 8.6× b) | 20.1 ms (0.91× b) | 18.3 ms |
| compensated | 2.35 ms (5.0× n, 7.8× b) | 2.75 ms (4.2× n, 6.6× b) | 5.2 ms (2.2× n, 3.5× b) | 11.7 ms (1.6× b) | 18.3 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | avx512 | avx2 | sse2 | none | base | blas |
|---|---|---|---|---|---|---|
| fast | 4.47 ms (3.2× n, 6.4× b) | 4.36 ms (3.3× n, 6.5× b) | 4.43 ms (3.2× n, 6.4× b) | 14.3 ms (2.0× b) | 28.4 ms | 10 ms (2.8× b) |
| pairwise | 4.66 ms (3.1× n, 6.0× b) | 4.49 ms (3.2× n, 6.2× b) | 4.6 ms (3.2× n, 6.1× b) | 14.5 ms (1.9× b) | 28.1 ms | 10 ms (2.8× b) |
| compensated | 4.49 ms (3.1× n, 6.2× b) | 4.29 ms (3.2× n, 6.5× b) | 5.81 ms (2.4× n, 4.8× b) | 13.8 ms (2.0× b) | 27.8 ms | 8.51 ms (3.3× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 27.8 ms (6.3× n, 6.7× b) | 45.2 ms (3.9× n, 4.1× b) | 128 ms (1.4× n, 1.5× b) | 176 ms (1.1× b) | 187 ms |
| fast | 16.6 ms (10.6× n, 11.0× b) | 29.1 ms (6.0× n, 6.3× b) | 62.7 ms (2.8× n, 2.9× b) | 176 ms (1.0× b) | 184 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 14.8 ms (3.1× n, 3.1× b) | 33.3 ms (1.4× n, 1.4× b) | 92 ms (0.49× n, 0.50× b) | 45.2 ms (1.0× b) | 46.2 ms |
| fast | 10.6 ms (5.0× n, 4.3× b) | 21.1 ms (2.5× n, 2.1× b) | 47.8 ms (1.1× n, 0.95× b) | 52.7 ms (0.86× b) | 45.2 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 37.5 ms (3.8× n, 3.9× b) | 58.7 ms (2.5× n, 2.5× b) | 180 ms (0.80× n, 0.81× b) | 144 ms (1.0× b) | 145 ms |
| fast | 17.7 ms (8.1× n, 8.2× b) | 25.3 ms (5.7× n, 5.7× b) | 60.8 ms (2.4× n, 2.4× b) | 144 ms (1.0× b) | 145 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 38.2 ms (6.7× n, 7.0× b) | 63.3 ms (4.1× n, 4.2× b) | 192 ms (1.3× n, 1.4× b) | 257 ms (1.0× b) | 268 ms |
| fast | 19.5 ms (13.1× n, 13.8× b) | 29.3 ms (8.7× n, 9.2× b) | 80.4 ms (3.2× n, 3.3× b) | 255 ms (1.1× b) | 269 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 26.3 ms (8.1× n, 3.5× b) | 41.6 ms (5.1× n, 2.2× b) | 171 ms (1.3× n, 0.54× b) | 214 ms (0.43× b) | 92.3 ms |
| fast | 12.5 ms (17.1× n, 7.4× b) | 13.7 ms (15.6× n, 6.8× b) | 26.5 ms (8.1× n, 3.5× b) | 214 ms (0.43× b) | 92.3 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 60.4 ms (2.3× n, 2.7× b) | 110 ms (1.2× n, 1.5× b) | 342 ms (0.40× n, 0.47× b) | 136 ms (1.2× b) | 160 ms |
| fast | 60.3 ms (2.3× n, 2.7× b) | 110 ms (1.2× n, 1.5× b) | 342 ms (0.40× n, 0.47× b) | 136 ms (1.2× b) | 160 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 65.7 ms (1.3× n, 1.4× b) | 109 ms (0.79× n, 0.82× b) | 312 ms (0.28× n, 0.29× b) | 85.9 ms (1.0× b) | 89.2 ms |
| fast | 65.7 ms (1.3× n, 1.4× b) | 109 ms (0.79× n, 0.82× b) | 312 ms (0.28× n, 0.29× b) | 85.9 ms (1.0× b) | 89.2 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 65.1 ms (4.5× n, 4.4× b) | 72.8 ms (4.0× n, 3.9× b) | 186 ms (1.6× n, 1.5× b) | 294 ms (1.0× b) | 286 ms |
| fast | 46.9 ms (6.3× n, 6.1× b) | 32.6 ms (9.0× n, 8.8× b) | 66.5 ms (4.4× n, 4.3× b) | 294 ms (1.0× b) | 287 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 54.7 ms (5.5× n, 5.4× b) | 60.8 ms (5.0× n, 4.8× b) | 164 ms (1.8× n, 1.8× b) | 303 ms (1.0× b) | 295 ms |
| fast | 46.5 ms (6.5× n, 6.3× b) | 50.1 ms (6.0× n, 5.9× b) | 124 ms (2.4× n, 2.4× b) | 303 ms (1.0× b) | 295 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 129 ms (3.8× n, 3.7× b) | 195 ms (2.5× n, 2.5× b) | 532 ms (0.91× n, 0.90× b) | 485 ms (1.0× b) | 479 ms |
| fast | 82.2 ms (5.9× n, 5.8× b) | 101 ms (4.8× n, 4.8× b) | 239 ms (2.0× n, 2.0× b) | 485 ms (1.0× b) | 479 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 79.3 ms (6.4× n, 6.3× b) | 93.7 ms (5.4× n, 5.3× b) | 265 ms (1.9× n, 1.9× b) | 508 ms (1.0× b) | 501 ms |
| fast | 71.8 ms (7.1× n, 7.0× b) | 82.3 ms (6.2× n, 6.1× b) | 225 ms (2.3× n, 2.2× b) | 509 ms (1.0× b) | 501 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 150 ms (5.4× n, 5.4× b) | 215 ms (3.8× n, 3.8× b) | 552 ms (1.5× n, 1.5× b) | 812 ms (1.0× b) | 811 ms |
| fast | 114 ms (7.2× n, 7.2× b) | 152 ms (5.3× n, 5.3× b) | 362 ms (2.2× n, 2.2× b) | 813 ms (1.0× b) | 813 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 80.5 ms (2.2× n, 2.3× b) | 109 ms (1.6× n, 1.7× b) | 204 ms (0.87× n, 0.92× b) | 179 ms (1.1× b) | 188 ms |
| fast | 68.4 ms (2.6× n, 2.7× b) | 82.1 ms (2.2× n, 2.2× b) | 145 ms (1.2× n, 1.3× b) | 178 ms (1.0× b) | 183 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | avx512 | avx2 | sse2 | none | base |
|---|---|---|---|---|---|
| accurate | 178 ms (4.2× n, 4.1× b) | 283 ms (2.6× n, 2.6× b) | 738 ms (1.0× n, 1.0× b) | 744 ms (1.0× b) | 737 ms |
| fast | 117 ms (6.4× n, 6.3× b) | 162 ms (4.6× n, 4.5× b) | 394 ms (1.9× n, 1.9× b) | 744 ms (1.0× b) | 737 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.17 µs | 1.09 µs (0.16× b) | 0.701 µs (0.24× b) |
| 1e1 | 0.18 µs | 1.1 µs (0.16× b) | 0.711 µs (0.25× b) |
| 1e2 | 0.34 µs | 1.12 µs (0.30× b) | 0.751 µs (0.45× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.0701 µs | 1.31 µs (0.05× b) |
| 1e1 | 0.1 µs | 1.42 µs (0.07× b) |
| 1e2 | 0.18 µs | 1.46 µs (0.12× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.0711 µs | 0.962 µs (0.07× b) |
| 1e1 | 0.12 µs | 1 µs (0.12× b) |
| 1e2 | 0.541 µs | 1.09 µs (0.50× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.0701 µs | 1.3 µs (0.05× b) |
| 1e1 | 0.14 µs | 1.33 µs (0.11× b) |
| 1e2 | 0.25 µs | 1.4 µs (0.18× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.0701 µs | 1.3 µs (0.05× b) |
| 1e1 | 0.0801 µs | 1.31 µs (0.06× b) |
| 1e2 | 0.18 µs | 1.36 µs (0.13× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.25 µs | 1.47 µs (0.17× b) |
| 1e1 | 0.391 µs | 1.48 µs (0.26× b) |
| 1e2 | 0.511 µs | 1.51 µs (0.34× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.1 µs | 1.52 µs (0.07× b) |
| 1e1 | 0.131 µs | 1.57 µs (0.08× b) |
| 1e2 | 0.3 µs | 1.73 µs (0.17× b) |

### add (integer64)

`simd_add(xi64, yi64)` vs base `xi64 + yi64`

| n | base | auto |
|---|---|---|
| 1e0 | 38.1 µs | 3.9 µs (9.8× b) |
| 1e1 | 36.6 µs | 3.79 µs (9.7× b) |
| 1e2 | 37.5 µs | 3.88 µs (9.7× b) |

