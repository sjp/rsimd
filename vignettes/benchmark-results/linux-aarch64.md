# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-07T07:00:10Z |
| rsimd | 0.1.0 (git 8a2d3e839cdb544674eb5fd938b645cac69469fd) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `gcc -std=gnu2x` |
| CFLAGS | `-g -O2` |
| OS | Linux 6.17.0-1022-azure (aarch64) |
| CPU | implementer 0x41 part 0xd49, 4 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16, sve, sve2, sve128 |
| Tiers | sve2, sve, neon, none (auto: sve2) |
| bench | 1.1.4 |
| Run time | 1014 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.48 µs (2.2× n, 8.5× b) | 1.45 µs (2.2× n, 8.7× b) | 1.6 µs (2.0× n, 7.9× b) | 3.22 µs (3.9× b) | 12.6 µs |
| 1e5 | 29.4 µs (6.7× n, 44.4× b) | 30.3 µs (6.5× n, 43.2× b) | 33.3 µs (6.0× n, 39.3× b) | 198 µs (6.6× b) | 1.31 ms |
| 1e7 | 3.85 ms (5.2× n, 33.2× b) | 3.65 ms (5.5× n, 35.0× b) | 3.87 ms (5.2× n, 33.0× b) | 20 ms (6.4× b) | 128 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.44 µs (1.3× n, 0.58× b) | 1.42 µs (1.3× n, 0.58× b) | 1.46 µs (1.3× n, 0.57× b) | 1.86 µs (0.45× b) | 0.832 µs |
| 1e5 | 18.2 µs (3.3× n, 3.2× b) | 20.1 µs (3.0× n, 2.9× b) | 22.3 µs (2.7× n, 2.7× b) | 60.5 µs (1.0× b) | 59.2 µs |
| 1e7 | 1.78 ms (3.6× n, 3.6× b) | 1.99 ms (3.2× n, 3.2× b) | 2.18 ms (3.0× n, 3.0× b) | 6.44 ms (1.0× b) | 6.45 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.51 µs (2.1× n, 27.1× b) | 1.52 µs (2.1× n, 26.9× b) | 1.58 µs (2.0× n, 25.9× b) | 3.23 µs (12.7× b) | 41 µs |
| 1e5 | 29.4 µs (6.7× n, 124.9× b) | 30.5 µs (6.5× n, 120.4× b) | 33.3 µs (5.9× n, 110.2× b) | 198 µs (18.5× b) | 3.67 ms |
| 1e7 | 4.1 ms (4.9× n, 87.4× b) | 3.85 ms (5.2× n, 93.1× b) | 4.03 ms (5.0× n, 89.0× b) | 20 ms (17.9× b) | 358 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | sve2 | sve | neon | none | base | blas |
|---|---|---|---|---|---|---|
| 1e3 | 2.03 µs (1.7× n, 7.3× b) | 2.17 µs (1.6× n, 6.8× b) | 2.22 µs (1.6× n, 6.7× b) | 3.48 µs (4.3× b) | 14.8 µs | 1.28 µs (11.6× b) |
| 1e5 | 47.7 µs (3.9× n, 28.8× b) | 49.1 µs (3.8× n, 28.0× b) | 58.9 µs (3.1× n, 23.3× b) | 185 µs (7.4× b) | 1.37 ms | 84.6 µs (16.3× b) |
| 1e7 | 7.18 ms (2.7× n, 21.4× b) | 7.37 ms (2.6× n, 20.8× b) | 7.21 ms (2.7× n, 21.3× b) | 19.2 ms (8.0× b) | 153 ms | 15.6 ms (9.8× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.5 µs (0.81× n, 0.22× b) | 3.46 µs (0.82× n, 0.22× b) | 2.27 µs (1.3× n, 0.34× b) | 2.84 µs (0.27× b) | 0.76 µs |
| 1e5 | 210 µs (1.3× n, 1.0× b) | 206 µs (1.4× n, 1.0× b) | 214 µs (1.3× n, 0.94× b) | 282 µs (0.71× b) | 200 µs |
| 1e7 | 25.2 ms (1.2× n, 1.0× b) | 25.2 ms (1.2× n, 1.0× b) | 25.9 ms (1.2× n, 1.0× b) | 31 ms (0.80× b) | 24.7 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.4 µs (1.7× n, 0.81× b) | 2.95 µs (1.4× n, 0.66× b) | 3.02 µs (1.3× n, 0.64× b) | 4.06 µs (0.48× b) | 1.94 µs |
| 1e5 | 110 µs (1.9× n, 3.7× b) | 101 µs (2.1× n, 4.0× b) | 122 µs (1.7× n, 3.3× b) | 208 µs (2.0× b) | 407 µs |
| 1e7 | 15.9 ms (1.6× n, 3.0× b) | 16.1 ms (1.6× n, 2.9× b) | 16.6 ms (1.5× n, 2.8× b) | 25.2 ms (1.9× b) | 46.9 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.54 µs (1.4× n, 0.67× b) | 2.5 µs (1.4× n, 0.68× b) | 4.05 µs (0.88× n, 0.42× b) | 3.55 µs (0.48× b) | 1.7 µs |
| 1e5 | 216 µs (1.5× n, 1.1× b) | 218 µs (1.5× n, 1.1× b) | 366 µs (0.87× n, 0.67× b) | 317 µs (0.78× b) | 246 µs |
| 1e7 | 27.9 ms (1.2× n, 1.2× b) | 27.6 ms (1.3× n, 1.2× b) | 37.8 ms (0.92× n, 0.87× b) | 34.9 ms (0.94× b) | 33 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.34 µs (0.88× n, 1.7× b) | 2.98 µs (1.3× n, 2.5× b) | 3.02 µs (1.3× n, 2.4× b) | 3.8 µs (1.9× b) | 7.32 µs |
| 1e5 | 204 µs (1.4× n, 3.0× b) | 208 µs (1.4× n, 2.9× b) | 210 µs (1.4× n, 2.9× b) | 296 µs (2.0× b) | 606 µs |
| 1e7 | 25.1 ms (1.3× n, 2.5× b) | 25.2 ms (1.3× n, 2.5× b) | 25 ms (1.3× n, 2.5× b) | 32 ms (2.0× b) | 62.6 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 5.67 µs (1.1× n, 0.82× b) | 5.84 µs (1.1× n, 0.79× b) | 6.13 µs (1.0× n, 0.76× b) | 6.26 µs (0.74× b) | 4.64 µs |
| 1e5 | 607 µs (1.1× n, 1.0× b) ⚠ | 613 µs (1.1× n, 0.95× b) ⚠ | 652 µs (1.0× n, 0.89× b) ⚠ | 663 µs (0.88× b) | 582 µs |
| 1e7 | 61 ms (1.1× n, 1.0× b) ⚠ | 61.4 ms (1.1× n, 1.0× b) ⚠ | 65.6 ms (1.0× n, 0.91× b) ⚠ | 66.7 ms (0.89× b) | 59.4 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.04 µs (1.1× n, 0.43× b) | 1.06 µs (1.1× n, 0.42× b) | 1.03 µs (1.1× n, 0.44× b) | 1.15 µs (0.39× b) | 0.448 µs |
| 1e5 | 23.3 µs (1.6× n, 1.5× b) | 23.1 µs (1.6× n, 1.5× b) | 28.6 µs (1.3× n, 1.2× b) | 37.1 µs (0.93× b) | 34.5 µs |
| 1e7 | 3.84 ms (1.1× n, 1.1× b) | 3.6 ms (1.2× n, 1.2× b) | 3.86 ms (1.1× n, 1.1× b) | 4.37 ms (1.0× b) | 4.35 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.96 µs (1.2× n, 0.45× b) | 0.984 µs (1.2× n, 0.44× b) | 0.952 µs (1.3× n, 0.45× b) | 1.2 µs (0.36× b) | 0.432 µs |
| 1e5 | 13.7 µs (2.3× n, 2.2× b) | 14.4 µs (2.2× n, 2.1× b) | 9.61 µs (3.3× n, 3.2× b) | 31.6 µs (1.0× b) | 30.6 µs |
| 1e7 | 1.49 ms (2.4× n, 2.5× b) | 1.46 ms (2.4× n, 2.5× b) | 1 ms (3.5× n, 3.7× b) | 3.51 ms (1.1× b) | 3.7 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.58 µs (1.0× n, 0.56× b) | 1.49 µs (1.1× n, 0.59× b) | 1.76 µs (0.93× n, 0.50× b) | 1.65 µs (0.53× b) | 0.881 µs |
| 1e5 | 53.8 µs (2.5× n, 2.4× b) | 109 µs (1.2× n, 1.2× b) | 127 µs (1.0× n, 1.0× b) | 133 µs (1.0× b) | 130 µs |
| 1e7 | 13.9 ms (1.1× n, 1.1× b) | 14.3 ms (1.1× n, 1.1× b) | 15.4 ms (1.0× n, 1.0× b) | 15.8 ms (1.0× b) | 15.7 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.06 µs (1.1× n, 1.2× b) | 1.95 µs (1.2× n, 1.3× b) | 2.01 µs (1.2× n, 1.2× b) | 2.34 µs (1.1× b) | 2.48 µs |
| 1e5 | 107 µs (1.3× n, 2.7× b) | 104 µs (1.4× n, 2.8× b) | 101 µs (1.4× n, 2.8× b) | 143 µs (2.0× b) | 288 µs |
| 1e7 | 14.1 ms (1.2× n, 2.5× b) | 14 ms (1.2× n, 2.5× b) | 13.9 ms (1.2× n, 2.6× b) | 16.9 ms (2.1× b) | 35.8 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.4 µs (1.4× n, 0.58× b) | 1.31 µs (1.5× n, 0.62× b) | 1.3 µs (1.5× n, 0.63× b) | 1.9 µs (0.43× b) | 0.816 µs |
| 1e5 | 23 µs (3.4× n, 2.6× b) | 23.8 µs (3.3× n, 2.5× b) | 21 µs (3.7× n, 2.8× b) | 77.2 µs (0.77× b) | 59.2 µs |
| 1e7 | 2.5 ms (3.3× n, 2.6× b) | 2.57 ms (3.2× n, 2.5× b) | 2.3 ms (3.6× n, 2.8× b) | 8.24 ms (0.78× b) | 6.42 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.75 µs (1.4× n, 0.93× b) | 1.69 µs (1.4× n, 1.0× b) | 2.03 µs (1.2× n, 0.80× b) | 2.44 µs (0.67× b) | 1.63 µs |
| 1e5 | 37.2 µs (2.7× n, 2.7× b) | 37.6 µs (2.6× n, 2.7× b) | 49.3 µs (2.0× n, 2.0× b) | 98.7 µs (1.0× b) | 100 µs |
| 1e7 | 4.21 ms (2.6× n, 6.6× b) | 4.29 ms (2.5× n, 6.5× b) | 5.51 ms (2.0× n, 5.1× b) | 10.7 ms (2.6× b) | 27.9 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.4 µs (1.7× n, 1.5× b) | 1.52 µs (1.5× n, 1.4× b) | 1.65 µs (1.4× n, 1.3× b) | 2.33 µs (0.89× b) | 2.06 µs |
| 1e5 | 36.2 µs (2.7× n, 2.9× b) | 36.2 µs (2.7× n, 2.9× b) | 52.3 µs (1.9× n, 2.0× b) | 98.6 µs (1.1× b) | 105 µs |
| 1e7 | 4.44 ms (2.5× n, 6.7× b) | 4.42 ms (2.5× n, 6.7× b) | 6.04 ms (1.8× n, 4.9× b) | 11.2 ms (2.7× b) | 29.7 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.9 µs (2.4× n, 0.31× b) | 3.78 µs (2.5× n, 0.32× b) | 4.08 µs (2.3× n, 0.30× b) | 9.35 µs (0.13× b) | 1.22 µs |
| 1e5 | 112 µs (7.0× n, 1.5× b) | 115 µs (6.8× n, 1.4× b) | 213 µs (3.7× n, 0.77× b) | 783 µs (0.21× b) | 165 µs |
| 1e7 | 20 ms (3.9× n, 1.0× b) | 20 ms (3.9× n, 1.0× b) | 22.1 ms (3.5× n, 0.87× b) | 78 ms (0.25× b) | 19.2 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.96 µs (1.4× n, 0.85× b) | 2.88 µs (1.4× n, 0.87× b) | 3.03 µs (1.4× n, 0.83× b) | 4.12 µs (0.61× b) | 2.5 µs |
| 1e5 | 59.4 µs (3.1× n, 4.2× b) | 57.2 µs (3.2× n, 4.4× b) | 73.1 µs (2.5× n, 3.4× b) | 185 µs (1.3× b) | 249 µs |
| 1e7 | 7.71 ms (2.4× n, 3.7× b) | 7.53 ms (2.4× n, 3.8× b) | 7.75 ms (2.3× n, 3.7× b) | 18.2 ms (1.6× b) | 28.4 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.6 µs (1.6× n, 1.0× b) | 2.57 µs (1.6× n, 1.0× b) | 3.01 µs (1.4× n, 0.88× b) | 4.06 µs (0.65× b) | 2.66 µs |
| 1e5 | 33.1 µs (5.5× n, 8.2× b) | 32.5 µs (5.6× n, 8.3× b) | 40.6 µs (4.5× n, 6.7× b) | 182 µs (1.5× b) | 270 µs |
| 1e7 | 3.4 ms (5.2× n, 8.2× b) | 3.4 ms (5.2× n, 8.2× b) | 4.13 ms (4.3× n, 6.7× b) | 17.8 ms (1.6× b) | 27.8 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.92 µs (2.9× n, 1.8× b) | 1.88 µs (2.9× n, 1.8× b) | 1.97 µs (2.8× n, 1.8× b) | 5.5 µs (0.63× b) | 3.46 µs |
| 1e5 | 144 µs (3.6× n, 2.0× b) | 138 µs (3.7× n, 2.1× b) | 149 µs (3.5× n, 1.9× b) | 513 µs (0.56× b) | 287 µs |
| 1e7 | 16.4 ms (3.2× n, 3.3× b) | 16.2 ms (3.2× n, 3.4× b) | 17 ms (3.0× n, 3.2× b) | 51.6 ms (1.1× b) | 54.8 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.21 µs (2.2× n, 11.8× b) | 2.18 µs (2.2× n, 12.0× b) | 2.43 µs (2.0× n, 10.8× b) | 4.8 µs (5.5× b) | 26.2 µs |
| 1e5 | 93.4 µs (5.9× n, 25.6× b) | 90.3 µs (6.1× n, 26.5× b) | 117 µs (4.7× n, 20.5× b) | 547 µs (4.4× b) | 2.39 ms |
| 1e7 | 19.4 ms (3.2× n, 17.2× b) | 18.8 ms (3.3× n, 17.8× b) | 20.8 ms (3.0× n, 16.0× b) | 61.7 ms (5.4× b) | 334 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.71 µs (0.76× n, 0.48× b) | 2.37 µs (1.2× n, 0.76× b) | 2.77 µs (1.0× n, 0.65× b) | 2.81 µs (0.64× b) | 1.8 µs |
| 1e5 | 101 µs (1.5× n, 1.5× b) | 102 µs (1.5× n, 1.5× b) | 153 µs (1.0× n, 1.0× b) | 153 µs (1.0× b) | 149 µs |
| 1e7 | 27.3 ms (1.1× n, 1.1× b) | 27.2 ms (1.1× n, 1.1× b) | 31.4 ms (1.0× n, 1.0× b) | 30.6 ms (1.0× b) | 30.5 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.82 µs (1.0× n, 1.7× b) | 4.8 µs (1.0× n, 1.7× b) | 4.82 µs (1.0× n, 1.7× b) | 4.87 µs (1.7× b) | 8.38 µs |
| 1e5 | 376 µs (1.4× n, 2.6× b) | 510 µs (1.0× n, 1.9× b) | 512 µs (1.0× n, 1.9× b) | 510 µs (1.9× b) | 960 µs |
| 1e7 | 51.1 ms (1.0× n, 1.9× b) | 51.2 ms (1.0× n, 1.9× b) | 51.3 ms (1.0× n, 1.9× b) | 51.1 ms (1.9× b) | 96.1 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.38 µs (2.0× n, 3.5× b) | 2.38 µs (2.0× n, 3.5× b) | 3.35 µs (1.4× n, 2.5× b) | 4.82 µs (1.7× b) | 8.36 µs |
| 1e5 | 263 µs (1.9× n, 3.6× b) | 262 µs (1.9× n, 3.7× b) | 358 µs (1.4× n, 2.7× b) | 508 µs (1.9× b) | 958 µs |
| 1e7 | 28 ms (1.8× n, 3.4× b) | 27.8 ms (1.8× n, 3.5× b) | 36.4 ms (1.4× n, 2.6× b) | 51 ms (1.9× b) | 95.8 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 48.3 µs (1.3× n, 0.36× b) | 46.5 µs (1.4× n, 0.38× b) | 49.7 µs (1.3× n, 0.35× b) | 62.8 µs (0.28× b) | 17.5 µs |
| 1e5 | 4.09 ms (1.3× n, 0.46× b) | 4.04 ms (1.4× n, 0.46× b) | 4.33 ms (1.3× n, 0.43× b) | 5.5 ms (0.34× b) | 1.86 ms |
| 1e7 | 406 ms (1.4× n, 0.46× b) | 400 ms (1.4× n, 0.47× b) | 435 ms (1.3× n, 0.43× b) | 563 ms (0.33× b) | 187 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 6.28 µs (1.1× n, 0.30× b) | 6.3 µs (1.1× n, 0.30× b) | 6.54 µs (1.1× n, 0.29× b) | 7.2 µs (0.26× b) | 1.88 µs |
| 1e5 | 472 µs (1.2× n, 1.0× b) | 479 µs (1.2× n, 0.95× b) | 519 µs (1.1× n, 0.88× b) | 564 µs (0.80× b) | 454 µs |
| 1e7 | 52.8 ms (1.2× n, 1.0× b) | 52.8 ms (1.1× n, 1.0× b) | 56.9 ms (1.1× n, 0.94× b) | 60.7 ms (0.88× b) | 53.7 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 10.7 µs (1.1× n, 0.60× b) | 10.8 µs (1.1× n, 0.60× b) | 11.3 µs (1.1× n, 0.57× b) | 12.1 µs (0.53× b) | 6.45 µs |
| 1e5 | 930 µs (1.1× n, 1.0× b) | 942 µs (1.1× n, 1.0× b) | 1 ms (1.0× n, 0.92× b) ⚠ | 1.04 ms (0.89× b) | 929 µs |
| 1e7 | 91.7 ms (1.1× n, 1.0× b) | 91.6 ms (1.1× n, 1.0× b) | 98.7 ms (1.0× n, 0.95× b) ⚠ | 103 ms (0.91× b) | 93.3 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.11 µs (2.5× n, 34.7× b) | 3.18 µs (2.5× n, 34.0× b) | 3.52 µs (2.2× n, 30.7× b) | 7.83 µs (13.8× b) | 108 µs |
| 1e5 | 174 µs (3.8× n, 63.4× b) | 174 µs (3.8× n, 63.3× b) | 211 µs (3.1× n, 52.2× b) | 665 µs (16.6× b) | 11 ms |
| 1e7 | 20 ms (3.3× n, 55.1× b) | 19.7 ms (3.4× n, 55.8× b) | 22 ms (3.0× n, 50.1× b) | 66.1 ms (16.7× b) | 1.1 s |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 11.6 µs (0.76× n, 0.53× b) | 11.6 µs (0.76× n, 0.53× b) | 12.9 µs (0.68× n, 0.48× b) | 8.82 µs (0.70× b) | 6.16 µs |
| 1e5 | 1.21 ms (0.86× n, 0.72× b) ⚠ | 1.21 ms (0.86× n, 0.72× b) ⚠ | 1.34 ms (0.78× n, 0.65× b) ⚠ | 1.04 ms (0.84× b) | 872 µs |
| 1e7 | 120 ms (0.86× n, 0.73× b) ⚠ | 120 ms (0.86× n, 0.73× b) ⚠ | 133 ms (0.78× n, 0.66× b) ⚠ | 104 ms (0.84× b) | 87.4 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 6.97 µs (1.9× n, 0.80× b) | 6.93 µs (1.9× n, 0.81× b) | 8.68 µs (1.5× n, 0.64× b) | 13 µs (0.43× b) | 5.58 µs |
| 1e5 | 644 µs (1.9× n, 1.3× b) | 641 µs (2.0× n, 1.3× b) | 780 µs (1.6× n, 1.1× b) | 1.25 ms (0.66× b) | 822 µs |
| 1e7 | 65.1 ms (1.9× n, 1.3× b) | 64.9 ms (1.9× n, 1.3× b) | 76.4 ms (1.6× n, 1.1× b) | 125 ms (0.66× b) | 82.5 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.58 µs (2.1× n, 8.1× b) | 1.58 µs (2.1× n, 8.1× b) | 1.73 µs (1.9× n, 7.4× b) | 3.33 µs (3.9× b) | 12.8 µs |
| 1e5 | 33.3 µs (5.9× n, 39.4× b) | 33.6 µs (5.9× n, 39.0× b) | 47.2 µs (4.2× n, 27.8× b) | 198 µs (6.6× b) | 1.31 ms |
| 1e7 | 3.83 ms (5.2× n, 33.4× b) | 3.85 ms (5.2× n, 33.2× b) | 5.26 ms (3.8× n, 24.3× b) | 19.9 ms (6.4× b) | 128 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.6 µs (1.4× n, 0.57× b) | 1.55 µs (1.4× n, 0.58× b) | 1.66 µs (1.3× n, 0.54× b) | 2.2 µs (0.41× b) | 0.904 µs |
| 1e5 | 35.2 µs (2.8× n, 1.7× b) | 35.7 µs (2.8× n, 1.7× b) | 41.5 µs (2.4× n, 1.4× b) | 98.6 µs (0.61× b) | 59.8 µs |
| 1e7 | 3.82 ms (2.7× n, 1.8× b) | 3.71 ms (2.8× n, 1.9× b) | 4.66 ms (2.2× n, 1.5× b) | 10.3 ms (0.68× b) | 7.03 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.74 µs (1.3× n, 0.48× b) | 1.38 µs (1.6× n, 0.60× b) | 1.37 µs (1.6× n, 0.61× b) | 2.23 µs (0.37× b) | 0.832 µs |
| 1e5 | 15.4 µs (6.0× n, 3.9× b) | 15.4 µs (6.0× n, 3.9× b) | 17.1 µs (5.4× n, 3.5× b) | 91.9 µs (0.65× b) | 59.9 µs |
| 1e7 | 1.51 ms (6.2× n, 4.3× b) | 1.5 ms (6.2× n, 4.4× b) | 1.68 ms (5.6× n, 3.9× b) | 9.31 ms (0.70× b) | 6.54 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.08 µs (2.0× n, 2.1× b) | 1.08 µs (2.0× n, 2.1× b) | 1.22 µs (1.8× n, 1.8× b) | 2.14 µs (1.0× b) | 2.24 µs |
| 1e5 | 27.2 µs (4.7× n, 6.4× b) | 26.8 µs (4.8× n, 6.5× b) | 28.7 µs (4.5× n, 6.1× b) | 129 µs (1.4× b) | 174 µs |
| 1e7 | 3.08 ms (4.2× n, 8.4× b) | 3.3 ms (4.0× n, 7.9× b) | 3.24 ms (4.0× n, 8.0× b) | 13.1 ms (2.0× b) | 25.9 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.39 µs (1.3× n, 0.70× b) | 2.35 µs (1.3× n, 0.71× b) | 2.26 µs (1.4× n, 0.74× b) | 3.15 µs (0.53× b) | 1.68 µs |
| 1e5 | 17.3 µs (5.5× n, 8.7× b) | 8.67 µs (10.9× n, 17.3× b) | 11.6 µs (8.1× n, 12.9× b) | 94.8 µs (1.6× b) | 150 µs |
| 1e7 | 2.64 ms (4.2× n, 6.5× b) | 2.62 ms (4.2× n, 6.6× b) | 2.61 ms (4.2× n, 6.6× b) | 11 ms (1.6× b) | 17.2 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.62 µs (1.3× n, 8.9× b) | 1.61 µs (1.3× n, 9.0× b) | 1.65 µs (1.3× n, 8.8× b) | 2.08 µs (6.9× b) | 14.4 µs |
| 1e5 | 5.91 µs (10.3× n, 220.1× b) | 5.5 µs (11.0× n, 236.7× b) | 9.09 µs (6.7× n, 143.2× b) | 60.6 µs (21.5× b) | 1.3 ms |
| 1e7 | 424 µs (13.9× n, 490.9× b) | 376 µs (15.7× n, 554.3× b) | 747 µs (7.9× n, 278.6× b) | 5.91 ms (35.2× b) | 208 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.54 µs (2.1× n, 7.2× b) | 1.62 µs (2.0× n, 6.9× b) | 1.62 µs (2.0× n, 6.8× b) | 3.3 µs (3.4× b) | 11.1 µs |
| 1e5 | 2.14 µs (4.4× n, 514.0× b) | 2.14 µs (4.4× n, 514.0× b) | 2.58 µs (3.7× n, 427.2× b) | 9.52 µs (115.8× b) | 1.1 ms |
| 1e7 | 2.21 µs (4.3× n, 49590.7× b) | 2.19 µs (4.4× n, 49952.7× b) | 2.6 µs (3.7× n, 42113.9× b) | 9.57 µs (11444.0× b) | 109 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.44 µs (0.93× n, 0.19× b) | 1.42 µs (0.94× n, 0.19× b) | 1.45 µs (0.92× n, 0.19× b) | 1.34 µs (0.20× b) | 0.272 µs |
| 1e5 | 1.91 µs (0.64× n, 0.13× b) | 2.01 µs (0.61× n, 0.12× b) | 2.1 µs (0.58× n, 0.12× b) | 1.22 µs (0.20× b) | 0.248 µs |
| 1e7 | 1.9 µs (0.67× n, 0.16× b) | 1.97 µs (0.65× n, 0.15× b) | 2.06 µs (0.62× n, 0.14× b) | 1.28 µs (0.23× b) | 0.296 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.856 µs (1.0× n, 0.19× b) | 1.02 µs (0.88× n, 0.16× b) | 0.936 µs (1.0× n, 0.17× b) | 0.896 µs (0.18× b) | 0.16 µs |
| 1e5 | 0.84 µs (1.0× n, 0.17× b) | 0.848 µs (1.0× n, 0.17× b) | 0.848 µs (1.0× n, 0.17× b) | 0.848 µs (0.17× b) | 0.144 µs |
| 1e7 | 0.776 µs (1.0× n, 0.21× b) | 0.8 µs (1.0× n, 0.20× b) | 0.792 µs (1.0× n, 0.20× b) | 0.792 µs (0.20× b) | 0.16 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.848 µs (1.0× n, 0.16× b) | 0.832 µs (1.0× n, 0.16× b) | 0.832 µs (1.0× n, 0.16× b) | 0.856 µs (0.16× b) | 0.136 µs |
| 1e5 | 0.792 µs (1.1× n, 0.15× b) | 0.824 µs (1.0× n, 0.15× b) | 0.84 µs (1.0× n, 0.14× b) | 0.856 µs (0.14× b) | 0.12 µs |
| 1e7 | 0.776 µs (1.0× n, 0.20× b) | 0.776 µs (1.0× n, 0.20× b) | 0.768 µs (1.0× n, 0.20× b) | 0.792 µs (0.19× b) | 0.152 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.71 µs (2.0× n, 7.2× b) | 1.67 µs (2.1× n, 7.4× b) | 1.82 µs (1.9× n, 6.8× b) | 3.45 µs (3.6× b) | 12.3 µs |
| 1e5 | 33.2 µs (6.1× n, 39.0× b) | 33.8 µs (6.1× n, 38.4× b) | 47.4 µs (4.3× n, 27.3× b) | 204 µs (6.3× b) | 1.29 ms |
| 1e7 | 3.57 ms (5.7× n, 35.6× b) | 3.62 ms (5.6× n, 35.1× b) | 5.11 ms (4.0× n, 24.8× b) | 20.2 ms (6.3× b) | 127 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.59 µs (1.4× n, 0.52× b) | 1.83 µs (1.2× n, 0.45× b) | 1.74 µs (1.3× n, 0.47× b) | 2.25 µs (0.37× b) | 0.824 µs |
| 1e5 | 2.5 µs (2.1× n, 25.3× b) | 2.5 µs (2.1× n, 25.3× b) | 2.82 µs (1.9× n, 22.5× b) | 5.22 µs (12.1× b) | 63.5 µs |
| 1e7 | 2.52 µs (2.1× n, 3042.4× b) | 2.53 µs (2.1× n, 3038.2× b) | 2.82 µs (1.8× n, 2719.8× b) | 5.19 µs (1479.3× b) | 7.68 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.37 µs (1.6× n, 0.19× b) | 1.3 µs (1.7× n, 0.20× b) | 1.76 µs (1.2× n, 0.15× b) | 2.2 µs (0.12× b) | 0.256 µs |
| 1e5 | 1.6 µs (3.1× n, 0.15× b) | 1.61 µs (3.1× n, 0.15× b) | 1.86 µs (2.7× n, 0.13× b) | 4.94 µs (0.05× b) | 0.248 µs |
| 1e7 | 1.6 µs (3.0× n, 0.17× b) | 1.61 µs (3.0× n, 0.17× b) | 1.82 µs (2.7× n, 0.15× b) | 4.88 µs (0.06× b) | 0.28 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| fast | 3.01 ms (6.6× n, 42.4× b) | 3.25 ms (6.1× n, 39.3× b) | 3.48 ms (5.7× n, 36.7× b) | 19.8 ms (6.5× b) | 128 ms |
| pairwise | 3.07 ms (6.8× n, 41.6× b) | 3.05 ms (6.9× n, 41.8× b) | 3.78 ms (5.6× n, 33.7× b) | 21 ms (6.1× b) | 128 ms |
| compensated | 7.36 ms (2.4× n, 17.3× b) | 7.7 ms (2.3× n, 16.6× b) | 11.5 ms (1.6× n, 11.1× b) | 18 ms (7.1× b) | 128 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | sve2 | sve | neon | none | base | blas |
|---|---|---|---|---|---|---|
| fast | 6.83 ms (2.8× n, 22.4× b) | 6.82 ms (2.8× n, 22.5× b) | 7.01 ms (2.7× n, 21.9× b) | 19.2 ms (8.0× b) | 153 ms | 14.8 ms (10.4× b) |
| pairwise | 6.79 ms (2.9× n, 22.6× b) | 6.78 ms (2.9× n, 22.6× b) | 7.62 ms (2.6× n, 20.2× b) | 20 ms (7.7× b) | 153 ms | 15.8 ms (9.7× b) |
| compensated | 9.54 ms (2.2× n, 16.1× b) | 9.65 ms (2.2× n, 15.9× b) | 14.5 ms (1.4× n, 10.6× b) | 20.8 ms (7.4× b) | 154 ms | 16 ms (9.6× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 137 ms (1.3× n, 1.3× b) | 136 ms (1.3× n, 1.3× b) | 139 ms (1.3× n, 1.3× b) | 178 ms (1.0× b) | 184 ms |
| fast | 75.5 ms (2.4× n, 2.4× b) | 75.3 ms (2.4× n, 2.4× b) | 94.6 ms (1.9× n, 1.9× b) | 177 ms (1.0× b) | 183 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 64 ms (1.0× n, 1.1× b) | 64 ms (1.0× n, 1.1× b) | 64 ms (1.0× n, 1.1× b) | 63.9 ms (1.1× b) | 67.3 ms |
| fast | 70.9 ms (0.89× n, 0.94× b) | 68.4 ms (0.93× n, 1.0× b) | 75 ms (0.84× n, 0.89× b) | 63.4 ms (1.1× b) | 66.8 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 190 ms (1.0× n, 1.0× b) | 190 ms (1.0× n, 1.0× b) | 169 ms (1.1× n, 1.1× b) | 194 ms (1.0× b) | 190 ms |
| fast | 81.5 ms (2.4× n, 2.3× b) | 81.3 ms (2.4× n, 2.3× b) | 88.2 ms (2.2× n, 2.2× b) | 194 ms (1.0× b) | 191 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 166 ms (1.5× n, 1.6× b) | 166 ms (1.5× n, 1.6× b) | 173 ms (1.5× n, 1.5× b) | 258 ms (1.0× b) | 259 ms |
| fast | 91.1 ms (2.8× n, 2.8× b) | 90.5 ms (2.8× n, 2.9× b) | 102 ms (2.5× n, 2.5× b) | 258 ms (1.0× b) | 259 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 123 ms (1.0× n, 0.94× b) | 123 ms (1.0× n, 0.94× b) | 144 ms (0.82× n, 0.81× b) | 118 ms (1.0× b) | 116 ms |
| fast | 52.7 ms (2.3× n, 2.2× b) | 52.6 ms (2.3× n, 2.2× b) | 55.3 ms (2.2× n, 2.1× b) | 119 ms (1.0× b) | 116 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 220 ms (1.0× n, 0.84× b) | 220 ms (1.0× n, 0.84× b) | 223 ms (1.0× n, 0.83× b) | 220 ms (0.84× b) | 186 ms |
| fast | 220 ms (1.0× n, 0.84× b) | 220 ms (1.0× n, 0.84× b) | 222 ms (1.0× n, 0.84× b) | 221 ms (0.84× b) | 186 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 128 ms (1.0× n, 1.0× b) | 127 ms (1.0× n, 1.0× b) | 128 ms (1.0× n, 1.0× b) | 127 ms (1.0× b) | 126 ms |
| fast | 128 ms (1.0× n, 1.0× b) | 128 ms (1.0× n, 1.0× b) | 128 ms (1.0× n, 1.0× b) | 128 ms (1.0× b) | 126 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 196 ms (1.2× n, 1.2× b) | 196 ms (1.2× n, 1.2× b) | 210 ms (1.1× n, 1.1× b) | 232 ms (1.0× b) | 236 ms |
| fast | 128 ms (1.8× n, 1.8× b) | 128 ms (1.8× n, 1.8× b) | 130 ms (1.8× n, 1.8× b) | 231 ms (1.0× b) | 236 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 180 ms (1.8× n, 1.7× b) | 180 ms (1.8× n, 1.7× b) | 182 ms (1.8× n, 1.7× b) | 320 ms (1.0× b) | 313 ms |
| fast | 150 ms (2.1× n, 2.1× b) | 149 ms (2.2× n, 2.1× b) | 166 ms (1.9× n, 1.9× b) | 320 ms (1.0× b) | 314 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 476 ms (1.0× n, 1.0× b) | 477 ms (1.0× n, 1.0× b) | 498 ms (1.0× n, 0.93× b) | 474 ms (1.0× b) | 463 ms |
| fast | 256 ms (1.9× n, 1.8× b) | 253 ms (1.9× n, 1.8× b) | 288 ms (1.6× n, 1.6× b) | 475 ms (1.0× b) | 462 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 282 ms (2.1× n, 2.0× b) | 283 ms (2.1× n, 2.0× b) | 264 ms (2.2× n, 2.2× b) | 588 ms (1.0× b) | 574 ms |
| fast | 249 ms (2.4× n, 2.3× b) | 249 ms (2.4× n, 2.3× b) | 252 ms (2.3× n, 2.3× b) | 589 ms (1.0× b) | 574 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 514 ms (1.4× n, 1.4× b) | 526 ms (1.4× n, 1.4× b) | 537 ms (1.4× n, 1.3× b) | 726 ms (1.0× b) | 718 ms |
| fast | 376 ms (1.9× n, 1.9× b) | 379 ms (1.9× n, 1.9× b) | 421 ms (1.7× n, 1.7× b) | 726 ms (1.0× b) | 718 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 231 ms (1.0× n, 0.94× b) | 232 ms (1.0× n, 0.94× b) | 244 ms (0.91× n, 0.89× b) | 222 ms (1.0× b) | 217 ms |
| fast | 167 ms (1.3× n, 1.3× b) | 168 ms (1.3× n, 1.3× b) | 194 ms (1.1× n, 1.1× b) | 222 ms (1.0× b) | 217 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 672 ms (1.2× n, 1.1× b) | 666 ms (1.2× n, 1.2× b) | 696 ms (1.1× n, 1.1× b) | 775 ms (1.0× b) | 768 ms |
| fast | 396 ms (2.0× n, 1.9× b) | 396 ms (2.0× n, 1.9× b) | 456 ms (1.7× n, 1.7× b) | 774 ms (1.0× b) | 768 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.24 µs | 1.12 µs (0.21× b) | 0.784 µs (0.31× b) |
| 1e1 | 0.32 µs | 1.12 µs (0.29× b) | 0.776 µs (0.41× b) |
| 1e2 | 1.3 µs | 1.22 µs (1.1× b) | 0.832 µs (1.6× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 1.34 µs (0.09× b) |
| 1e1 | 0.152 µs | 1.24 µs (0.12× b) |
| 1e2 | 0.36 µs | 1.54 µs (0.23× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 0.92 µs (0.13× b) |
| 1e1 | 0.16 µs | 0.96 µs (0.17× b) |
| 1e2 | 0.696 µs | 1.55 µs (0.45× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.112 µs | 1.4 µs (0.08× b) |
| 1e1 | 0.184 µs | 1.26 µs (0.15× b) |
| 1e2 | 0.32 µs | 1.62 µs (0.20× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 1.5 µs (0.08× b) |
| 1e1 | 0.144 µs | 1.44 µs (0.10× b) |
| 1e2 | 0.256 µs | 1.73 µs (0.15× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.328 µs | 1.41 µs (0.23× b) |
| 1e1 | 0.46 µs | 1.42 µs (0.32× b) |
| 1e2 | 1.6 µs | 1.57 µs (1.0× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.144 µs | 3.62 µs (0.04× b) |
| 1e1 | 0.288 µs | 4.72 µs (0.06× b) |
| 1e2 | 0.512 µs | 4.45 µs (0.12× b) |

