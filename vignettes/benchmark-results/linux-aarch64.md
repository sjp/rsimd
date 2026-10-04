# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-04T11:04:16Z |
| rsimd | 0.0.0.9000 (git b364156c6eeb1e54c2f8e636498ec0766063c691-dirty) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `aarch64-linux-gnu-gcc` |
| CFLAGS | `-g -O2 -ffile-prefix-map=/build/reproducible-path/r-base-4.6.1=. -fstack-protector-strong -fstack-clash-protection -Wformat -Werror=format-security -mbranch-protection=standard -Wdate-time -D_FORTIFY_SOURCE=2` |
| OS | Linux 7.0.14-orbstack-00380-ga7e0a2dc9535 (aarch64) |
| CPU | implementer 0x61 part 0x000, 8 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16 |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 58 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.83 µs (1.4× n, 2.3× b) | 5.33 µs (1.7× b) | 8.92 µs |
| 1e5 | 17.2 µs (11.1× n, 70.1× b) | 190 µs (6.3× b) | 1.21 ms |
| 1e7 | 1.45 ms (12.7× n, 83.1× b) | 18.4 ms (6.5× b) | 120 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.71 µs (1.1× n, 0.18× b) | 4.17 µs (0.16× b) | 0.667 µs |
| 1e5 | 21 µs (2.9× n, 2.6× b) | 59.8 µs (0.93× b) | 55.3 µs |
| 1e7 | 1.69 ms (3.3× n, 3.2× b) | 5.53 ms (1.0× b) | 5.37 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.83 µs (1.4× n, 5.6× b) | 6.79 µs (4.0× b) | 26.9 µs |
| 1e5 | 19.4 µs (10.1× n, 162.0× b) | 197 µs (16.0× b) | 3.15 ms |
| 1e7 | 1.47 ms (11.8× n, 202.5× b) | 17.4 ms (17.2× b) | 298 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 6.58 µs (1.2× n, 1.4× b) | 7.75 µs (1.2× b) | 9.54 µs | 1.71 µs (5.6× b) |
| 1e5 | 25 µs (6.3× n, 48.0× b) | 159 µs (7.6× b) | 1.2 ms | 150 µs (8.0× b) |
| 1e7 | 2.41 ms (6.3× n, 49.4× b) | 15.1 ms (7.9× b) | 119 ms | 14.4 ms (8.3× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.58 µs (1.3× n, 0.20× b) | 4.56 µs (0.16× b) | 0.709 µs |
| 1e5 | 24.8 µs (4.2× n, 1.2× b) | 103 µs (0.30× b) | 30.8 µs |
| 1e7 | 3.73 ms (2.9× n, 1.1× b) | 10.8 ms (0.37× b) | 3.96 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.62 µs (1.2× n, 0.37× b) | 4.29 µs (0.31× b) | 1.33 µs |
| 1e5 | 37.1 µs (3.2× n, 7.4× b) | 119 µs (2.3× b) | 276 µs |
| 1e7 | 3.63 ms (3.1× n, 9.9× b) | 11.2 ms (3.2× b) | 36.1 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4 µs (1.3× n, 0.18× b) | 5.04 µs (0.14× b) | 0.709 µs |
| 1e5 | 39 µs (3.8× n, 1.6× b) | 148 µs (0.41× b) | 60.8 µs |
| 1e7 | 4.92 ms (3.3× n, 1.9× b) | 16.2 ms (0.59× b) | 9.51 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.79 µs (1.0× n, 1.2× b) | 3.88 µs (1.2× b) | 4.71 µs |
| 1e5 | 29.7 µs (3.1× n, 14.3× b) | 91.6 µs (4.6× b) | 425 µs |
| 1e7 | 4.78 ms (1.8× n, 10.6× b) | 8.68 ms (5.8× b) | 50.8 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.75 µs (1.3× n, 0.45× b) | 6.17 µs (0.34× b) | 2.13 µs |
| 1e5 | 192 µs (1.7× n, 1.0× b) | 328 µs (0.61× b) | 200 µs |
| 1e7 | 18.1 ms (1.7× n, 1.1× b) | 31.3 ms (0.65× b) | 20.3 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.25 µs (1.1× n, 0.27× b) | 1.38 µs (0.24× b) | 0.333 µs |
| 1e5 | 11 µs (2.7× n, 2.6× b) | 29.8 µs (1.0× b) | 28.3 µs |
| 1e7 | 1.26 ms (2.2× n, 2.1× b) | 2.72 ms (1.0× b) | 2.7 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.29 µs (1.2× n, 0.26× b) | 1.5 µs (0.22× b) | 0.333 µs |
| 1e5 | 6.38 µs (4.6× n, 4.4× b) | 29.6 µs (0.94× b) | 27.8 µs |
| 1e7 | 604 µs (4.5× n, 4.5× b) | 2.74 ms (1.0× b) | 2.7 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.29 µs (1.3× n, 0.22× b) | 3 µs (0.17× b) | 0.501 µs |
| 1e5 | 27.5 µs (3.2× n, 1.0× b) | 88 µs (0.32× b) | 28.6 µs |
| 1e7 | 2.87 ms (3.2× n, 1.2× b) | 9.04 ms (0.40× b) | 3.58 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 7.54 µs (1.3× n, 0.12× b) | 9.71 µs (0.09× b) | 0.876 µs |
| 1e5 | 61 µs (10.0× n, 1.1× b) | 608 µs (0.11× b) | 67.6 µs |
| 1e7 | 5.62 ms (11.1× n, 1.4× b) | 62.4 ms (0.12× b) | 7.64 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.92 µs (1.4× n, 1.9× b) | 5.33 µs (1.4× b) | 7.42 µs |
| 1e5 | 17.4 µs (11.2× n, 56.1× b) | 195 µs (5.0× b) | 975 µs |
| 1e7 | 1.56 ms (11.8× n, 63.9× b) | 18.4 ms (5.4× b) | 99.8 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.77 µs (1.2× n, 0.06× b) | 4.42 µs (0.05× b) | 0.208 µs |
| 1e5 | 22 µs (2.7× n, 0.02× b) | 60.2 µs (<0.01× b) | 0.333 µs |
| 1e7 | 1.77 ms (3.2× n, <0.01× b) | 5.72 ms (<0.01× b) | 0.167 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.13 µs (1.1× n, 0.07× b) | 1.21 µs (0.07× b) | 0.084 µs |
| 1e5 | 1.13 µs (1.0× n, 0.04× b) | 1.12 µs (0.04× b) | 0.042 µs |
| 1e7 | 1.08 µs (1.1× n, 0.08× b) | 1.17 µs (0.07× b) | 0.083 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.2× n, 0.08× b) | 1.33 µs (0.06× b) | 0.083 µs |
| 1e5 | 1.12 µs (1.4× n, 0.15× b) | 1.54 µs (0.11× b) | 0.166 µs |
| 1e7 | 1.12 µs (1.0× n, 0.07× b) | 1.12 µs (0.07× b) | 0.084 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.48 ms (12.5× n, 78.3× b) | 18.5 ms (6.3× b) | 116 ms |
| pairwise | 1.33 ms (14.1× n, 86.5× b) | 18.8 ms (6.1× b) | 115 ms |
| compensated | 4.03 ms (5.3× n, 28.5× b) | 21.4 ms (5.4× b) | 115 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 2.75 ms (5.6× n, 46.9× b) | 15.4 ms (8.4× b) | 129 ms | 14.5 ms (8.9× b) |
| pairwise | 2.44 ms (6.5× n, 48.6× b) | 16 ms (7.4× b) | 119 ms | 14.4 ms (8.3× b) |
| compensated | 5.22 ms (3.1× n, 22.7× b) | 16.1 ms (7.4× b) | 119 ms | 14.4 ms (8.2× b) |

