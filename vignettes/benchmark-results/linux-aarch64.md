# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-11T22:13:57Z |
| rsimd | 0.1.0 (git d95c7863d200e98e30b6f9679a28d91f9d7d1615) |
| rsimd library | `/home/runner/work/_temp/Library/rsimd` (built R 4.6.1; aarch64-unknown-linux-gnu; 2026-10-11 21:53:59 UTC; unix) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `gcc -std=gnu2x` |
| CFLAGS | `-g -O2` |
| OS | Linux 6.17.0-1022-azure (aarch64) |
| CPU | implementer 0x41 part 0xd49, 4 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16, sve, sve2, sve128 |
| Tiers | sve2, sve, neon, none (auto: sve2) |
| bench | 1.1.4 |
| Run time | 1059 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.43 µs (2.9× n, 8.8× b) | 1.43 µs (2.9× n, 8.8× b) | 1.4 µs (2.9× n, 9.0× b) | 4.11 µs (3.1× b) | 12.6 µs |
| 1e5 | 19.5 µs (14.5× n, 67.1× b) | 19.5 µs (14.5× n, 66.9× b) | 18.4 µs (15.3× n, 70.9× b) | 283 µs (4.6× b) | 1.31 ms |
| 1e7 | 2.45 ms (11.5× n, 52.1× b) | 2.42 ms (11.7× n, 52.7× b) | 2.44 ms (11.6× n, 52.3× b) | 28.3 ms (4.5× b) | 128 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.36 µs (1.3× n, 0.61× b) | 1.37 µs (1.3× n, 0.60× b) | 1.41 µs (1.3× n, 0.59× b) | 1.78 µs (0.46× b) | 0.824 µs |
| 1e5 | 20 µs (3.0× n, 3.0× b) | 21.4 µs (2.8× n, 2.8× b) | 23.9 µs (2.5× n, 2.5× b) | 60.6 µs (1.0× b) | 59.2 µs |
| 1e7 | 1.93 ms (3.2× n, 3.2× b) | 2.07 ms (3.0× n, 3.0× b) | 2.31 ms (2.7× n, 2.7× b) | 6.26 ms (1.0× b) | 6.14 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.5 µs (2.8× n, 27.4× b) | 1.47 µs (2.8× n, 27.9× b) | 1.4 µs (3.0× n, 29.3× b) | 4.18 µs (9.8× b) | 41 µs |
| 1e5 | 17.4 µs (16.2× n, 210.2× b) | 18.4 µs (15.3× n, 198.6× b) | 18.2 µs (15.5× n, 201.2× b) | 283 µs (12.9× b) | 3.66 ms |
| 1e7 | 2.6 ms (10.9× n, 137.6× b) | 2.53 ms (11.2× n, 141.2× b) | 2.36 ms (12.0× n, 151.7× b) | 28.3 ms (12.7× b) | 358 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | sve2 | sve | neon | none | base | blas |
|---|---|---|---|---|---|---|
| 1e3 | 1.82 µs (1.8× n, 8.2× b) | 1.81 µs (1.8× n, 8.3× b) | 1.81 µs (1.8× n, 8.3× b) | 3.34 µs (4.5× b) | 14.9 µs | 1.2 µs (12.5× b) |
| 1e5 | 32.3 µs (5.3× n, 47.1× b) | 32.6 µs (5.3× n, 46.7× b) | 32.2 µs (5.4× n, 47.2× b) | 173 µs (8.8× b) | 1.52 ms | 84.9 µs (17.9× b) |
| 1e7 | 6.02 ms (3.0× n, 25.3× b) | 5.81 ms (3.1× n, 26.3× b) | 5.95 ms (3.0× n, 25.6× b) | 18 ms (8.5× b) | 153 ms | 14 ms (10.9× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.26 µs (1.1× n, 0.31× b) | 2.22 µs (1.1× n, 0.31× b) | 2.33 µs (1.1× n, 0.30× b) | 2.48 µs (0.28× b) | 0.696 µs |
| 1e5 | 222 µs (1.1× n, 1.0× b) | 217 µs (1.1× n, 1.0× b) | 228 µs (1.1× n, 1.0× b) | 244 µs (0.89× b) | 217 µs |
| 1e7 | 24.4 ms (1.0× n, 1.0× b) | 24.5 ms (1.0× n, 1.0× b) | 24.9 ms (1.0× n, 1.0× b) | 25.5 ms (0.95× b) | 24.2 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.33 µs (1.5× n, 0.78× b) | 2.26 µs (1.5× n, 0.80× b) | 2.42 µs (1.4× n, 0.75× b) | 3.38 µs (0.53× b) | 1.81 µs |
| 1e5 | 170 µs (1.5× n, 2.7× b) | 163 µs (1.6× n, 2.8× b) | 167 µs (1.6× n, 2.7× b) | 260 µs (1.8× b) | 456 µs |
| 1e7 | 16.6 ms (1.5× n, 2.7× b) | 16.1 ms (1.5× n, 2.8× b) | 16.5 ms (1.5× n, 2.8× b) | 24.9 ms (1.8× b) | 45.3 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.49 µs (1.2× n, 0.54× b) | 2.42 µs (1.2× n, 0.55× b) | 2.58 µs (1.1× n, 0.52× b) | 2.89 µs (0.46× b) | 1.34 µs |
| 1e5 | 246 µs (1.1× n, 1.1× b) | 239 µs (1.1× n, 1.1× b) | 248 µs (1.1× n, 1.1× b) | 273 µs (1.0× b) | 266 µs |
| 1e7 | 27.5 ms (1.1× n, 1.1× b) | 26.8 ms (1.1× n, 1.1× b) | 27.9 ms (1.0× n, 1.1× b) | 29.1 ms (1.1× b) | 30.6 ms |

### scalar_sub (double)

`simd_sub(1, x)` vs base `1 - x`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.27 µs (1.7× n, 0.31× b) | 2.3 µs (1.7× n, 0.30× b) | 2.33 µs (1.6× n, 0.30× b) | 3.84 µs (0.18× b) | 0.696 µs |
| 1e5 | 220 µs (1.8× n, 0.94× b) | 225 µs (1.7× n, 0.92× b) | 227 µs (1.7× n, 0.91× b) | 388 µs (0.53× b) | 207 µs |
| 1e7 | 21.5 ms (1.7× n, 1.0× b) | 22.3 ms (1.7× n, 0.94× b) | 22.1 ms (1.7× n, 0.95× b) | 37 ms (0.57× b) | 20.9 ms |

### fma_scalar (double)

`simd_fma(x, 2, 1)` vs base `x * 2 + 1`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.4 µs (1.3× n, 0.48× b) | 2.36 µs (1.3× n, 0.49× b) | 2.5 µs (1.3× n, 0.46× b) | 3.15 µs (0.37× b) | 1.16 µs |
| 1e5 | 232 µs (1.4× n, 1.1× b) | 227 µs (1.4× n, 1.1× b) | 241 µs (1.3× n, 1.0× b) | 314 µs (0.79× b) | 249 µs |
| 1e7 | 22.7 ms (1.4× n, 1.1× b) | 22.1 ms (1.4× n, 1.2× b) | 23.6 ms (1.3× n, 1.1× b) | 30.6 ms (0.83× b) | 25.4 ms |

### lt_scalar (double)

`simd_lt(x, 0)` vs base `x < 0`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.19 µs (1.3× n, 0.66× b) | 2.18 µs (1.3× n, 0.67× b) | 2.28 µs (1.2× n, 0.64× b) | 2.83 µs (0.51× b) | 1.46 µs |
| 1e5 | 146 µs (1.5× n, 1.4× b) | 147 µs (1.5× n, 1.4× b) | 159 µs (1.3× n, 1.3× b) | 214 µs (1.0× b) | 205 µs |
| 1e7 | 15.1 ms (1.4× n, 1.4× b) | 15.1 ms (1.4× n, 1.4× b) | 16.5 ms (1.3× n, 1.3× b) | 21.4 ms (1.0× b) | 20.9 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.58 µs (1.2× n, 2.8× b) | 2.56 µs (1.2× n, 2.8× b) | 2.61 µs (1.2× n, 2.8× b) | 3.03 µs (2.4× b) | 7.22 µs |
| 1e5 | 231 µs (1.2× n, 2.8× b) | 230 µs (1.2× n, 2.8× b) | 232 µs (1.2× n, 2.7× b) | 273 µs (2.3× b) | 637 µs |
| 1e7 | 24.1 ms (1.1× n, 2.6× b) | 24.2 ms (1.1× n, 2.5× b) | 24.1 ms (1.1× n, 2.6× b) | 27.3 ms (2.3× b) | 61.7 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 5.62 µs (1.1× n, 0.82× b) | 5.67 µs (1.1× n, 0.82× b) | 6.14 µs (1.0× n, 0.75× b) | 6.27 µs (0.74× b) | 4.62 µs |
| 1e5 | 630 µs (1.1× n, 1.0× b) ⚠ | 640 µs (1.1× n, 1.0× b) ⚠ | 679 µs (1.0× n, 0.90× b) ⚠ | 691 µs (0.88× b) | 609 µs |
| 1e7 | 60.5 ms (1.1× n, 1.0× b) ⚠ | 61 ms (1.1× n, 1.0× b) ⚠ | 65.4 ms (1.0× n, 0.90× b) ⚠ | 66.5 ms (0.89× b) | 59 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.928 µs (1.2× n, 0.47× b) | 0.992 µs (1.1× n, 0.44× b) | 0.992 µs (1.1× n, 0.44× b) | 1.1 µs (0.40× b) | 0.44 µs |
| 1e5 | 23.6 µs (1.6× n, 1.5× b) | 24.9 µs (1.5× n, 1.4× b) | 29.2 µs (1.3× n, 1.2× b) | 37.8 µs (0.95× b) | 35.8 µs |
| 1e7 | 2.78 ms (1.6× n, 1.5× b) | 2.81 ms (1.5× n, 1.5× b) | 3.35 ms (1.3× n, 1.3× b) | 4.32 ms (1.0× b) | 4.24 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.832 µs (1.2× n, 0.51× b) | 0.84 µs (1.2× n, 0.51× b) | 0.8 µs (1.3× n, 0.53× b) | 1.03 µs (0.41× b) | 0.425 µs |
| 1e5 | 14.1 µs (2.3× n, 2.2× b) | 14.3 µs (2.2× n, 2.1× b) | 9.76 µs (3.3× n, 3.1× b) | 32.1 µs (1.0× b) | 30.6 µs |
| 1e7 | 1.47 ms (2.3× n, 2.2× b) | 1.48 ms (2.3× n, 2.2× b) | 945 µs (3.5× n, 3.4× b) | 3.33 ms (1.0× b) | 3.26 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.47 µs (1.1× n, 0.58× b) | 1.44 µs (1.1× n, 0.59× b) | 1.38 µs (1.2× n, 0.62× b) | 1.64 µs (0.52× b) | 0.856 µs |
| 1e5 | 122 µs (1.2× n, 1.2× b) | 120 µs (1.2× n, 1.2× b) | 114 µs (1.3× n, 1.2× b) | 145 µs (1.0× b) | 141 µs |
| 1e7 | 12.8 ms (1.2× n, 1.2× b) | 12.8 ms (1.2× n, 1.2× b) | 12.7 ms (1.2× n, 1.2× b) | 15.2 ms (1.0× b) | 15 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.96 µs (1.1× n, 1.2× b) | 1.95 µs (1.1× n, 1.2× b) | 1.93 µs (1.1× n, 1.2× b) | 2.22 µs (1.1× b) | 2.38 µs |
| 1e5 | 128 µs (1.3× n, 2.7× b) | 128 µs (1.3× n, 2.7× b) | 127 µs (1.3× n, 2.7× b) | 161 µs (2.2× b) | 346 µs |
| 1e7 | 13.5 ms (1.2× n, 2.6× b) | 13.6 ms (1.2× n, 2.6× b) | 13 ms (1.3× n, 2.7× b) | 16.8 ms (2.1× b) | 35.3 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.26 µs (1.5× n, 0.65× b) | 1.26 µs (1.5× n, 0.65× b) | 1.27 µs (1.5× n, 0.64× b) | 1.86 µs (0.44× b) | 0.816 µs |
| 1e5 | 23.4 µs (3.3× n, 2.5× b) | 23.8 µs (3.3× n, 2.5× b) | 21.1 µs (3.7× n, 2.8× b) | 78 µs (0.76× b) | 59.2 µs |
| 1e7 | 2.29 ms (3.5× n, 2.7× b) | 2.31 ms (3.4× n, 2.6× b) | 2 ms (4.0× n, 3.0× b) | 7.94 ms (0.77× b) | 6.08 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.76 µs (1.4× n, 0.86× b) | 1.77 µs (1.4× n, 0.86× b) | 1.81 µs (1.4× n, 0.84× b) | 2.48 µs (0.61× b) | 1.51 µs |
| 1e5 | 33.4 µs (3.0× n, 8.2× b) | 33.9 µs (2.9× n, 8.1× b) | 38.9 µs (2.5× n, 7.0× b) | 98.9 µs (2.8× b) | 274 µs |
| 1e7 | 3.55 ms (3.0× n, 8.0× b) | 3.54 ms (3.0× n, 8.0× b) | 4.21 ms (2.5× n, 6.7× b) | 10.5 ms (2.7× b) | 28.4 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.33 µs (1.7× n, 1.4× b) | 1.33 µs (1.7× n, 1.4× b) | 1.44 µs (1.6× n, 1.3× b) | 2.26 µs (0.80× b) | 1.82 µs |
| 1e5 | 33.5 µs (2.9× n, 8.3× b) | 33.3 µs (3.0× n, 8.4× b) | 41.9 µs (2.4× n, 6.7× b) | 98.6 µs (2.8× b) | 279 µs |
| 1e7 | 4.2 ms (2.7× n, 6.8× b) | 4.12 ms (2.7× n, 7.0× b) | 5.11 ms (2.2× n, 5.6× b) | 11.3 ms (2.5× b) | 28.7 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.74 µs (2.2× n, 0.25× b) | 4.89 µs (2.1× n, 0.24× b) | 5.1 µs (2.0× n, 0.23× b) | 10.5 µs (0.11× b) | 1.18 µs |
| 1e5 | 200 µs (4.0× n, 0.89× b) | 200 µs (4.0× n, 0.89× b) | 228 µs (3.5× n, 0.78× b) | 802 µs (0.22× b) | 178 µs |
| 1e7 | 19.8 ms (3.9× n, 0.94× b) | 19.6 ms (4.0× n, 0.95× b) | 22.3 ms (3.5× n, 0.84× b) | 77.7 ms (0.24× b) | 18.6 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.86 µs (1.4× n, 0.82× b) | 2.87 µs (1.4× n, 0.82× b) | 2.98 µs (1.4× n, 0.79× b) | 4.11 µs (0.57× b) | 2.35 µs |
| 1e5 | 56.9 µs (3.2× n, 4.6× b) | 57.2 µs (3.2× n, 4.6× b) | 70.6 µs (2.6× n, 3.7× b) | 183 µs (1.4× b) | 262 µs |
| 1e7 | 7.05 ms (2.6× n, 4.0× b) | 7.11 ms (2.6× n, 4.0× b) | 7.72 ms (2.4× n, 3.6× b) | 18.3 ms (1.5× b) | 28.1 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.58 µs (1.6× n, 1.0× b) | 2.62 µs (1.6× n, 1.0× b) | 2.66 µs (1.5× n, 1.0× b) | 4.08 µs (0.63× b) | 2.57 µs |
| 1e5 | 33 µs (5.5× n, 8.6× b) | 32.9 µs (5.5× n, 8.6× b) | 40.6 µs (4.5× n, 7.0× b) | 182 µs (1.6× b) | 284 µs |
| 1e7 | 3.38 ms (5.3× n, 8.2× b) | 3.35 ms (5.4× n, 8.2× b) | 4.31 ms (4.2× n, 6.4× b) | 17.9 ms (1.5× b) | 27.5 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.83 µs (2.9× n, 1.9× b) | 1.83 µs (2.9× n, 1.9× b) | 1.74 µs (3.1× n, 2.0× b) | 5.35 µs (0.64× b) | 3.43 µs |
| 1e5 | 150 µs (3.5× n, 3.5× b) | 150 µs (3.5× n, 3.5× b) | 143 µs (3.7× n, 3.7× b) | 527 µs (1.0× b) | 530 µs |
| 1e7 | 15.8 ms (3.3× n, 3.5× b) | 15.9 ms (3.2× n, 3.4× b) | 15.2 ms (3.4× n, 3.6× b) | 51.4 ms (1.1× b) | 54.6 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.11 µs (2.2× n, 12.2× b) | 2.12 µs (2.2× n, 12.1× b) | 2.17 µs (2.1× n, 11.9× b) | 4.63 µs (5.6× b) | 25.7 µs |
| 1e5 | 180 µs (3.5× n, 13.3× b) | 180 µs (3.5× n, 13.3× b) | 187 µs (3.4× n, 12.8× b) | 633 µs (3.8× b) | 2.4 ms |
| 1e7 | 18.9 ms (3.3× n, 17.7× b) | 18.9 ms (3.3× n, 17.7× b) | 19.8 ms (3.1× n, 16.9× b) | 61.8 ms (5.4× b) | 334 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.23 µs (1.2× n, 0.79× b) | 2.22 µs (1.2× n, 0.80× b) | 2.67 µs (1.0× n, 0.66× b) | 2.69 µs (0.66× b) | 1.77 µs |
| 1e5 | 101 µs (3.2× n, 3.1× b) | 275 µs (1.2× n, 1.1× b) | 329 µs (1.0× n, 1.0× b) | 320 µs (1.0× b) | 316 µs |
| 1e7 | 27.4 ms (1.1× n, 1.1× b) | 27.2 ms (1.1× n, 1.1× b) | 32.1 ms (1.0× n, 1.0× b) | 30.9 ms (1.0× b) | 30.8 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 4.74 µs (1.0× n, 1.8× b) | 4.77 µs (1.0× n, 1.8× b) | 4.78 µs (1.0× n, 1.8× b) | 4.78 µs (1.8× b) | 8.38 µs |
| 1e5 | 533 µs (1.0× n, 1.8× b) | 534 µs (1.0× n, 1.8× b) | 537 µs (1.0× n, 1.8× b) | 533 µs (1.8× b) | 980 µs |
| 1e7 | 51.3 ms (1.0× n, 1.9× b) | 51.3 ms (1.0× n, 1.9× b) | 51.5 ms (1.0× n, 1.9× b) | 51 ms (1.9× b) | 96.2 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 2.34 µs (2.0× n, 3.6× b) | 2.36 µs (2.0× n, 3.5× b) | 3.27 µs (1.5× n, 2.6× b) | 4.78 µs (1.8× b) | 8.38 µs |
| 1e5 | 287 µs (1.9× n, 3.4× b) | 288 µs (1.9× n, 3.4× b) | 384 µs (1.4× n, 2.5× b) | 535 µs (1.8× b) | 980 µs |
| 1e7 | 28.3 ms (1.8× n, 3.4× b) | 28.1 ms (1.8× n, 3.4× b) | 36.6 ms (1.4× n, 2.6× b) | 51.2 ms (1.9× b) | 96.1 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 45.7 µs (1.4× n, 0.38× b) | 45.8 µs (1.4× n, 0.38× b) | 49.4 µs (1.3× n, 0.35× b) | 62 µs (0.28× b) | 17.4 µs |
| 1e5 | 4.05 ms (1.4× n, 0.46× b) | 4.04 ms (1.4× n, 0.47× b) | 4.4 ms (1.3× n, 0.43× b) | 5.68 ms (0.33× b) | 1.88 ms |
| 1e7 | 400 ms (1.4× n, 0.46× b) | 399 ms (1.4× n, 0.47× b) | 435 ms (1.3× n, 0.43× b) | 562 ms (0.33× b) | 186 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.33 µs (1.0× n, 0.56× b) | 3.34 µs (1.0× n, 0.56× b) | 3.72 µs (0.89× n, 0.50× b) | 3.3 µs (0.57× b) | 1.87 µs |
| 1e5 | 494 µs (1.0× n, 1.0× b) | 497 µs (1.0× n, 1.0× b) | 544 µs (0.92× n, 0.92× b) | 502 µs (1.0× b) | 500 µs |
| 1e7 | 51.1 ms (1.0× n, 1.0× b) | 51.1 ms (1.0× n, 1.0× b) | 56.3 ms (0.93× n, 0.94× b) | 52.1 ms (1.0× b) | 52.7 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 7.7 µs (1.1× n, 0.84× b) | 7.72 µs (1.1× n, 0.84× b) | 8.43 µs (1.0× n, 0.77× b) | 8.39 µs (0.77× b) | 6.48 µs |
| 1e5 | 964 µs (1.1× n, 1.0× b) ⚠ | 966 µs (1.1× n, 1.0× b) ⚠ | 1.03 ms (1.0× n, 0.94× b) ⚠ | 1.02 ms (0.95× b) | 968 µs |
| 1e7 | 90.7 ms (1.1× n, 1.0× b) ⚠ | 91.1 ms (1.1× n, 1.0× b) ⚠ | 97.8 ms (1.0× n, 0.95× b) ⚠ | 97.3 ms (1.0× b) | 92.9 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 3.14 µs (1.8× n, 34.4× b) | 3.09 µs (1.8× n, 35.0× b) | 3.41 µs (1.7× n, 31.7× b) | 5.66 µs (19.1× b) | 108 µs |
| 1e5 | 175 µs (2.5× n, 63.1× b) | 174 µs (2.5× n, 63.2× b) | 212 µs (2.1× n, 52.1× b) | 443 µs (24.9× b) | 11 ms |
| 1e7 | 20.4 ms (2.2× n, 54.0× b) | 20.5 ms (2.2× n, 53.8× b) | 23.4 ms (1.9× n, 47.1× b) | 45 ms (24.5× b) | 1.1 s |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 11.6 µs (0.75× n, 0.52× b) | 11.6 µs (0.75× n, 0.52× b) | 12.9 µs (0.67× n, 0.46× b) | 8.68 µs (0.69× b) | 5.98 µs |
| 1e5 | 1.24 ms (0.86× n, 0.73× b) ⚠ | 1.23 ms (0.87× n, 0.73× b) ⚠ | 1.37 ms (0.78× n, 0.66× b) ⚠ | 1.07 ms (0.84× b) | 897 µs |
| 1e7 | 121 ms (0.86× n, 0.73× b) ⚠ | 120 ms (0.86× n, 0.73× b) ⚠ | 133 ms (0.78× n, 0.66× b) ⚠ | 104 ms (0.84× b) | 87.6 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 6.54 µs (1.3× n, 0.83× b) | 6.58 µs (1.3× n, 0.83× b) | 7.81 µs (1.1× n, 0.70× b) | 8.7 µs (0.63× b) | 5.45 µs |
| 1e5 | 687 µs (1.3× n, 1.3× b) | 688 µs (1.3× n, 1.2× b) | 806 µs (1.1× n, 1.1× b) | 896 µs (1.0× b) | 859 µs |
| 1e7 | 64.9 ms (1.3× n, 1.3× b) | 65.1 ms (1.3× n, 1.3× b) | 76.4 ms (1.1× n, 1.1× b) | 84.4 ms (1.0× b) | 82.4 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.63 µs (2.6× n, 7.7× b) | 1.5 µs (2.8× n, 8.4× b) | 1.67 µs (2.5× n, 7.6× b) | 4.17 µs (3.0× b) | 12.6 µs |
| 1e5 | 33.4 µs (8.5× n, 39.3× b) | 34.2 µs (8.3× n, 38.4× b) | 47.3 µs (6.0× n, 27.7× b) | 283 µs (4.6× b) | 1.31 ms |
| 1e7 | 3.65 ms (7.7× n, 34.9× b) | 3.62 ms (7.8× n, 35.3× b) | 5.29 ms (5.3× n, 24.1× b) | 28.2 ms (4.5× b) | 128 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.61 µs (1.5× n, 0.52× b) | 1.5 µs (1.6× n, 0.56× b) | 1.55 µs (1.5× n, 0.54× b) | 2.36 µs (0.36× b) | 0.84 µs |
| 1e5 | 28.4 µs (3.5× n, 2.1× b) | 28.6 µs (3.5× n, 2.1× b) | 32.7 µs (3.0× n, 1.8× b) | 98.9 µs (0.61× b) | 59.8 µs |
| 1e7 | 3.02 ms (3.7× n, 2.4× b) | 3.13 ms (3.5× n, 2.4× b) | 3.41 ms (3.2× n, 2.2× b) | 11.1 ms (0.67× b) | 7.38 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.46 µs (1.6× n, 0.56× b) | 1.47 µs (1.5× n, 0.56× b) | 1.48 µs (1.5× n, 0.56× b) | 2.28 µs (0.36× b) | 0.824 µs |
| 1e5 | 15.6 µs (5.9× n, 3.8× b) | 15.8 µs (5.9× n, 3.8× b) | 17.5 µs (5.3× n, 3.4× b) | 92.3 µs (0.65× b) | 59.7 µs |
| 1e7 | 1.55 ms (6.1× n, 4.3× b) | 1.54 ms (6.1× n, 4.4× b) | 1.7 ms (5.5× n, 4.0× b) | 9.42 ms (0.71× b) | 6.71 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.02 µs (2.0× n, 2.1× b) | 1.01 µs (2.1× n, 2.1× b) | 1.04 µs (2.0× n, 2.0× b) | 2.07 µs (1.0× b) | 2.1 µs |
| 1e5 | 26.8 µs (4.8× n, 9.3× b) | 27 µs (4.8× n, 9.2× b) | 29.7 µs (4.3× n, 8.4× b) | 129 µs (1.9× b) | 249 µs |
| 1e7 | 2.98 ms (4.4× n, 8.7× b) | 3.23 ms (4.1× n, 8.0× b) | 3.1 ms (4.2× n, 8.3× b) | 13.1 ms (2.0× b) | 25.8 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.6 µs (1.5× n, 1.0× b) | 1.61 µs (1.5× n, 1.0× b) | 1.58 µs (1.5× n, 1.0× b) | 2.42 µs (0.65× b) | 1.57 µs |
| 1e5 | 28.4 µs (4.0× n, 5.5× b) | 28.4 µs (4.0× n, 5.5× b) | 27.9 µs (4.0× n, 5.6× b) | 113 µs (1.4× b) | 155 µs |
| 1e7 | 2.67 ms (4.1× n, 5.7× b) | 2.67 ms (4.1× n, 5.7× b) | 2.6 ms (4.2× n, 5.9× b) | 11 ms (1.4× b) | 15.3 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.27 µs (1.4× n, 11.2× b) | 1.27 µs (1.4× n, 11.2× b) | 1.3 µs (1.4× n, 11.0× b) | 1.82 µs (7.8× b) | 14.2 µs |
| 1e5 | 5.64 µs (10.7× n, 230.0× b) | 5.2 µs (11.6× n, 249.5× b) | 8.82 µs (6.8× n, 147.0× b) | 60.4 µs (21.5× b) | 1.3 ms |
| 1e7 | 423 µs (14.0× n, 489.4× b) | 376 µs (15.7× n, 551.4× b) | 748 µs (7.9× n, 276.9× b) | 5.91 ms (35.0× b) | 207 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.5 µs (2.8× n, 7.2× b) | 1.5 µs (2.8× n, 7.2× b) | 1.53 µs (2.7× n, 7.1× b) | 4.17 µs (2.6× b) | 10.8 µs |
| 1e5 | 1.9 µs (6.9× n, 579.9× b) | 1.92 µs (6.8× n, 572.6× b) | 1.8 µs (7.2× n, 610.8× b) | 13 µs (84.5× b) | 1.1 ms |
| 1e7 | 1.93 µs (6.7× n, 56875.9× b) | 1.95 µs (6.6× n, 56176.6× b) | 1.86 µs (7.0× n, 58828.7× b) | 13 µs (8450.7× b) | 110 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.46 µs (0.90× n, 0.18× b) | 1.46 µs (0.90× n, 0.18× b) | 1.5 µs (0.87× n, 0.18× b) | 1.31 µs (0.20× b) | 0.264 µs |
| 1e5 | 2.01 µs (0.65× n, 0.12× b) | 2.04 µs (0.64× n, 0.12× b) | 2.24 µs (0.58× n, 0.11× b) | 1.3 µs (0.19× b) | 0.248 µs |
| 1e7 | 2.01 µs (0.67× n, 0.14× b) | 2.08 µs (0.65× n, 0.14× b) | 2.25 µs (0.60× n, 0.13× b) | 1.34 µs (0.21× b) | 0.288 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.744 µs (1.1× n, 0.22× b) | 0.736 µs (1.1× n, 0.22× b) | 0.752 µs (1.1× n, 0.21× b) | 0.792 µs (0.20× b) | 0.16 µs |
| 1e5 | 0.736 µs (1.0× n, 0.20× b) | 0.728 µs (1.0× n, 0.20× b) | 0.728 µs (1.0× n, 0.20× b) | 0.752 µs (0.19× b) | 0.144 µs |
| 1e7 | 0.784 µs (1.0× n, 0.20× b) | 0.768 µs (1.0× n, 0.21× b) | 0.776 µs (1.0× n, 0.21× b) | 0.792 µs (0.20× b) | 0.16 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 0.744 µs (1.0× n, 0.18× b) | 0.728 µs (1.0× n, 0.19× b) | 0.712 µs (1.0× n, 0.19× b) | 0.744 µs (0.18× b) | 0.136 µs |
| 1e5 | 0.728 µs (1.0× n, 0.16× b) | 0.72 µs (1.0× n, 0.17× b) | 0.72 µs (1.0× n, 0.17× b) | 0.736 µs (0.16× b) | 0.12 µs |
| 1e7 | 0.752 µs (1.0× n, 0.19× b) | 0.744 µs (1.0× n, 0.19× b) | 0.744 µs (1.0× n, 0.19× b) | 0.768 µs (0.19× b) | 0.144 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.7 µs (2.5× n, 7.3× b) | 1.7 µs (2.5× n, 7.3× b) | 1.86 µs (2.2× n, 6.7× b) | 4.19 µs (3.0× b) | 12.4 µs |
| 1e5 | 34.3 µs (8.4× n, 37.8× b) | 34.8 µs (8.3× n, 37.2× b) | 47.4 µs (6.1× n, 27.3× b) | 288 µs (4.5× b) | 1.3 ms |
| 1e7 | 3.96 ms (7.2× n, 32.0× b) | 3.87 ms (7.4× n, 32.8× b) | 5.42 ms (5.3× n, 23.4× b) | 28.7 ms (4.4× b) | 127 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.66 µs (1.4× n, 0.50× b) | 1.66 µs (1.4× n, 0.50× b) | 1.71 µs (1.4× n, 0.48× b) | 2.37 µs (0.35× b) | 0.824 µs |
| 1e5 | 2.45 µs (2.2× n, 25.9× b) | 2.44 µs (2.2× n, 25.9× b) | 2.69 µs (2.0× n, 23.5× b) | 5.48 µs (11.5× b) | 63.3 µs |
| 1e7 | 2.46 µs (2.2× n, 3113.6× b) | 2.47 µs (2.2× n, 3103.5× b) | 2.72 µs (2.0× n, 2824.1× b) | 5.49 µs (1397.9× b) | 7.67 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| 1e3 | 1.46 µs (1.5× n, 0.18× b) | 1.46 µs (1.5× n, 0.18× b) | 1.49 µs (1.5× n, 0.18× b) | 2.23 µs (0.12× b) | 0.264 µs |
| 1e5 | 1.88 µs (2.7× n, 0.13× b) | 1.86 µs (2.7× n, 0.13× b) | 2 µs (2.5× n, 0.12× b) | 5.06 µs (0.05× b) | 0.248 µs |
| 1e7 | 1.87 µs (2.7× n, 0.15× b) | 1.93 µs (2.6× n, 0.15× b) | 2.01 µs (2.5× n, 0.14× b) | 5.05 µs (0.06× b) | 0.28 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| fast | 2.64 ms (10.7× n, 48.4× b) | 2.57 ms (11.0× n, 49.6× b) | 2.7 ms (10.5× n, 47.3× b) | 28.3 ms (4.5× b) | 128 ms |
| pairwise | 3.06 ms (9.5× n, 41.7× b) | 3 ms (9.7× n, 42.6× b) | 2.8 ms (10.4× n, 45.7× b) | 29 ms (4.4× b) | 128 ms |
| compensated | 7.29 ms (2.5× n, 17.5× b) | 7.3 ms (2.5× n, 17.5× b) | 9.26 ms (2.0× n, 13.8× b) | 18.2 ms (7.0× b) | 128 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | sve2 | sve | neon | none | base | blas |
|---|---|---|---|---|---|---|
| fast | 6.63 ms (2.7× n, 22.9× b) | 6.68 ms (2.7× n, 22.7× b) | 6.14 ms (2.9× n, 24.7× b) | 18 ms (8.4× b) | 152 ms | 14.7 ms (10.4× b) |
| pairwise | 6.3 ms (3.0× n, 24.1× b) | 6.48 ms (2.9× n, 23.4× b) | 6.13 ms (3.1× n, 24.8× b) | 18.9 ms (8.0× b) | 152 ms | 14.6 ms (10.4× b) |
| compensated | 8.2 ms (2.5× n, 18.6× b) | 8.35 ms (2.4× n, 18.2× b) | 10.8 ms (1.9× n, 14.1× b) | 20.2 ms (7.5× b) | 152 ms | 14.1 ms (10.8× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 137 ms (1.3× n, 1.3× b) | 137 ms (1.3× n, 1.3× b) | 138 ms (1.3× n, 1.3× b) | 177 ms (1.0× b) | 183 ms |
| fast | 75.4 ms (2.4× n, 2.4× b) | 74.6 ms (2.4× n, 2.5× b) | 94.3 ms (1.9× n, 1.9× b) | 177 ms (1.0× b) | 183 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 63.2 ms (1.0× n, 1.1× b) | 63.3 ms (1.0× n, 1.1× b) | 63.2 ms (1.0× n, 1.1× b) | 63.2 ms (1.1× b) | 66.6 ms |
| fast | 71.1 ms (0.89× n, 0.93× b) | 67.4 ms (0.94× n, 1.0× b) | 74.8 ms (0.85× n, 0.89× b) | 63.3 ms (1.0× b) | 66.4 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 190 ms (1.0× n, 1.0× b) | 189 ms (1.0× n, 1.0× b) | 168 ms (1.1× n, 1.1× b) | 193 ms (1.0× b) | 190 ms |
| fast | 81.7 ms (2.4× n, 2.3× b) | 81 ms (2.4× n, 2.3× b) | 87.7 ms (2.2× n, 2.2× b) | 193 ms (1.0× b) | 190 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 166 ms (1.6× n, 1.6× b) | 166 ms (1.6× n, 1.6× b) | 174 ms (1.5× n, 1.5× b) | 259 ms (1.0× b) | 260 ms |
| fast | 90.8 ms (2.8× n, 2.9× b) | 91.2 ms (2.8× n, 2.9× b) | 101 ms (2.5× n, 2.6× b) | 258 ms (1.0× b) | 261 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 123 ms (1.0× n, 0.93× b) | 123 ms (1.0× n, 0.93× b) | 144 ms (0.81× n, 0.80× b) | 117 ms (1.0× b) | 114 ms |
| fast | 52.4 ms (2.2× n, 2.2× b) | 52.5 ms (2.2× n, 2.2× b) | 55 ms (2.1× n, 2.1× b) | 117 ms (1.0× b) | 115 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 174 ms (1.1× n, 1.1× b) | 174 ms (1.1× n, 1.1× b) | 175 ms (1.1× n, 1.1× b) | 186 ms (1.0× b) | 186 ms |
| fast | 174 ms (1.1× n, 1.1× b) | 176 ms (1.1× n, 1.1× b) | 176 ms (1.1× n, 1.1× b) | 186 ms (1.0× b) | 186 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× b) | 125 ms |
| fast | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× n, 1.0× b) | 127 ms (1.0× b) | 125 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 195 ms (1.2× n, 1.2× b) | 193 ms (1.2× n, 1.2× b) | 209 ms (1.1× n, 1.1× b) | 228 ms (1.0× b) | 235 ms |
| fast | 127 ms (1.8× n, 1.8× b) | 127 ms (1.8× n, 1.8× b) | 129 ms (1.8× n, 1.8× b) | 228 ms (1.0× b) | 235 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 178 ms (1.8× n, 1.8× b) | 178 ms (1.8× n, 1.8× b) | 180 ms (1.8× n, 1.7× b) | 319 ms (1.0× b) | 313 ms |
| fast | 149 ms (2.1× n, 2.1× b) | 148 ms (2.2× n, 2.1× b) | 164 ms (1.9× n, 1.9× b) | 319 ms (1.0× b) | 312 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 492 ms (1.0× n, 0.93× b) | 493 ms (1.0× n, 0.93× b) | 514 ms (0.91× n, 0.89× b) | 470 ms (1.0× b) | 459 ms |
| fast | 271 ms (1.7× n, 1.7× b) | 272 ms (1.7× n, 1.7× b) | 308 ms (1.5× n, 1.5× b) | 469 ms (1.0× b) | 459 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 278 ms (2.1× n, 2.1× b) | 281 ms (2.1× n, 2.0× b) | 262 ms (2.2× n, 2.2× b) | 586 ms (1.0× b) | 573 ms |
| fast | 247 ms (2.4× n, 2.3× b) | 248 ms (2.4× n, 2.3× b) | 249 ms (2.4× n, 2.3× b) | 586 ms (1.0× b) | 573 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 517 ms (1.4× n, 1.4× b) | 518 ms (1.4× n, 1.4× b) | 536 ms (1.3× n, 1.3× b) | 721 ms (1.0× b) | 716 ms |
| fast | 376 ms (1.9× n, 1.9× b) | 378 ms (1.9× n, 1.9× b) | 419 ms (1.7× n, 1.7× b) | 721 ms (1.0× b) | 716 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 230 ms (1.0× n, 0.94× b) | 231 ms (1.0× n, 0.94× b) | 243 ms (0.91× n, 0.89× b) | 221 ms (1.0× b) | 216 ms |
| fast | 167 ms (1.3× n, 1.3× b) | 167 ms (1.3× n, 1.3× b) | 192 ms (1.1× n, 1.1× b) | 221 ms (1.0× b) | 216 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | sve2 | sve | neon | none | base |
|---|---|---|---|---|---|
| accurate | 688 ms (1.1× n, 1.1× b) | 691 ms (1.1× n, 1.1× b) | 712 ms (1.1× n, 1.1× b) | 768 ms (1.0× b) | 769 ms |
| fast | 418 ms (1.8× n, 1.8× b) | 417 ms (1.8× n, 1.8× b) | 471 ms (1.6× n, 1.6× b) | 769 ms (1.0× b) | 769 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.232 µs | 1.26 µs (0.18× b) | 0.776 µs (0.30× b) |
| 1e1 | 0.32 µs | 1.3 µs (0.25× b) | 0.8 µs (0.40× b) |
| 1e2 | 1.26 µs | 1.3 µs (1.0× b) | 0.784 µs (1.6× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 1.34 µs (0.09× b) |
| 1e1 | 0.144 µs | 1.47 µs (0.10× b) |
| 1e2 | 0.296 µs | 1.66 µs (0.18× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 0.896 µs (0.13× b) |
| 1e1 | 0.168 µs | 0.96 µs (0.18× b) |
| 1e2 | 0.592 µs | 1.48 µs (0.40× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.104 µs | 1.31 µs (0.08× b) |
| 1e1 | 0.176 µs | 1.34 µs (0.13× b) |
| 1e2 | 0.312 µs | 1.5 µs (0.21× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.12 µs | 1.53 µs (0.08× b) |
| 1e1 | 0.136 µs | 1.55 µs (0.09× b) |
| 1e2 | 0.24 µs | 1.72 µs (0.14× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.328 µs | 1.5 µs (0.22× b) |
| 1e1 | 0.472 µs | 1.58 µs (0.30× b) |
| 1e2 | 1.51 µs | 1.54 µs (1.0× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.152 µs | 1.32 µs (0.12× b) |
| 1e1 | 0.184 µs | 1.45 µs (0.13× b) |
| 1e2 | 0.464 µs | 1.7 µs (0.27× b) |

### add (integer64)

`simd_add(xi64, yi64)` vs base `xi64 + yi64`

| n | base | auto |
|---|---|---|
| 1e0 | 47.2 µs | 4.82 µs (9.8× b) |
| 1e1 | 47.6 µs | 4.94 µs (9.6× b) |
| 1e2 | 48.5 µs | 5.22 µs (9.3× b) |

