# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-05T05:26:06Z |
| rsimd | 0.1.0 (git 798a247b52f06ac5c1e52445498c74dda4d88720-dirty) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `aarch64-linux-gnu-gcc` |
| CFLAGS | `-g -O2 -ffile-prefix-map=/build/reproducible-path/r-base-4.6.1=. -fstack-protector-strong -fstack-clash-protection -Wformat -Werror=format-security -mbranch-protection=standard -Wdate-time -D_FORTIFY_SOURCE=2` |
| OS | Linux 7.0.14-orbstack-00380-ga7e0a2dc9535 (aarch64) |
| CPU | implementer 0x61 part 0x000, 8 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16 |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 322 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.42 µs (1.5× n, 3.0× b) | 5.08 µs (2.0× b) | 10.4 µs |
| 1e5 | 17.1 µs (10.7× n, 71.1× b) | 183 µs (6.7× b) | 1.22 ms |
| 1e7 | 1.43 ms (13.1× n, 80.3× b) | 18.8 ms (6.1× b) | 115 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.46 µs (1.1× n, 0.19× b) | 3.79 µs (0.18× b) | 0.667 µs |
| 1e5 | 21.8 µs (2.7× n, 2.5× b) | 58.8 µs (0.92× b) | 54.1 µs |
| 1e7 | 1.68 ms (3.2× n, 3.2× b) | 5.36 ms (1.0× b) | 5.33 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.29 µs (1.4× n, 6.4× b) | 5.96 µs (4.6× b) | 27.6 µs |
| 1e5 | 18.9 µs (9.9× n, 161.6× b) | 187 µs (16.4× b) | 3.06 ms |
| 1e7 | 1.43 ms (12.6× n, 206.7× b) | 17.9 ms (16.4× b) | 295 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 6 µs (1.2× n, 1.8× b) | 7.38 µs (1.4× b) | 10.6 µs | 1.63 µs (6.5× b) |
| 1e5 | 24.7 µs (6.3× n, 47.0× b) | 156 µs (7.4× b) | 1.16 ms | 143 µs (8.1× b) |
| 1e7 | 2.44 ms (6.4× n, 49.5× b) | 15.6 ms (7.7× b) | 121 ms | 14.4 ms (8.4× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.42 µs (1.2× n, 0.16× b) | 4.13 µs (0.13× b) | 0.542 µs |
| 1e5 | 26.2 µs (3.6× n, 1.0× b) | 95.6 µs (0.29× b) | 27.4 µs |
| 1e7 | 3.65 ms (3.0× n, 1.1× b) | 10.8 ms (0.36× b) | 3.94 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.46 µs (1.2× n, 0.37× b) | 4.17 µs (0.31× b) | 1.29 µs |
| 1e5 | 34.8 µs (3.2× n, 7.6× b) | 110 µs (2.4× b) | 263 µs |
| 1e7 | 3.65 ms (3.1× n, 10.9× b) | 11.2 ms (3.6× b) | 39.9 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.83 µs (1.3× n, 0.18× b) | 4.88 µs (0.15× b) | 0.708 µs |
| 1e5 | 33 µs (4.1× n, 1.7× b) | 136 µs (0.42× b) | 56.7 µs |
| 1e7 | 4.79 ms (3.2× n, 1.5× b) | 15.2 ms (0.46× b) | 7.04 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.25 µs (1.1× n, 1.3× b) | 3.67 µs (1.2× b) | 4.31 µs |
| 1e5 | 25.5 µs (2.8× n, 15.5× b) | 71.5 µs (5.5× b) | 395 µs |
| 1e7 | 3.61 ms (2.2× n, 13.5× b) | 7.77 ms (6.3× b) | 48.7 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.04 µs (1.2× n, 0.41× b) | 6.08 µs (0.34× b) | 2.08 µs |
| 1e5 | 171 µs (1.6× n, 1.1× b) | 279 µs (0.67× b) | 187 µs |
| 1e7 | 17.5 ms (1.6× n, 1.2× b) | 28.3 ms (0.72× b) | 20.5 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.17 µs (1.1× n, 0.29× b) | 1.33 µs (0.25× b) | 0.334 µs |
| 1e5 | 10 µs (2.8× n, 2.7× b) | 28.3 µs (0.94× b) | 26.6 µs |
| 1e7 | 1.32 ms (2.0× n, 2.1× b) | 2.7 ms (1.0× b) | 2.74 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.13 µs (1.2× n, 0.30× b) | 1.33 µs (0.25× b) | 0.333 µs |
| 1e5 | 5.46 µs (5.2× n, 4.9× b) | 28.2 µs (0.95× b) | 26.7 µs |
| 1e7 | 634 µs (4.3× n, 4.2× b) | 2.7 ms (1.0× b) | 2.69 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.25 µs (1.4× n, 0.22× b) | 3.08 µs (0.16× b) | 0.5 µs |
| 1e5 | 23.6 µs (4.7× n, 1.1× b) | 110 µs (0.25× b) | 27.1 µs |
| 1e7 | 2.67 ms (4.2× n, 1.2× b) | 11.2 ms (0.29× b) | 3.25 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 6.42 µs (1.4× n, 0.14× b) | 8.87 µs (0.10× b) | 0.875 µs |
| 1e5 | 54.7 µs (10.3× n, 1.2× b) | 566 µs (0.12× b) | 66.7 µs |
| 1e7 | 5.26 ms (11.5× n, 1.4× b) | 60.4 ms (0.12× b) | 7.26 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.63 µs (1.2× n, 0.36× b) | 5.46 µs (0.31× b) | 1.67 µs |
| 1e5 | 31.8 µs (2.2× n, 4.4× b) | 69.5 µs (2.0× b) | 139 µs |
| 1e7 | 2.76 ms (2.5× n, 5.3× b) | 6.9 ms (2.1× b) | 14.7 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.29 µs (1.4× n, 0.41× b) | 6 µs (0.30× b) | 1.77 µs |
| 1e5 | 18.3 µs (3.2× n, 7.6× b) | 58.9 µs (2.4× b) | 139 µs |
| 1e7 | 1.38 ms (4.3× n, 10.6× b) | 5.89 ms (2.5× b) | 14.7 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.13 µs (1.3× n, 0.60× b) | 4 µs (0.47× b) | 1.88 µs |
| 1e5 | 30.4 µs (13.1× n, 5.5× b) | 398 µs (0.42× b) | 166 µs |
| 1e7 | 3.25 ms (13.6× n, 5.6× b) | 44.2 ms (0.41× b) | 18.3 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.29 µs (1.4× n, 4.1× b) | 4.71 µs (2.9× b) | 13.6 µs |
| 1e5 | 50.4 µs (10.7× n, 24.7× b) | 538 µs (2.3× b) | 1.25 ms |
| 1e7 | 5.21 ms (10.9× n, 26.2× b) | 56.9 ms (2.4× b) | 136 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.29 µs (1.1× n, 0.18× b) | 3.54 µs (0.16× b) | 0.584 µs |
| 1e5 | 30.4 µs (1.0× n, 0.89× b) | 31.2 µs (0.86× b) | 27 µs |
| 1e7 | 3.89 ms (1.0× n, 0.95× b) | 3.9 ms (0.95× b) | 3.69 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.42 µs (1.1× n, 1.4× b) | 3.83 µs (1.2× b) | 4.67 µs |
| 1e5 | 47.6 µs (1.8× n, 9.8× b) | 86.1 µs (5.4× b) | 466 µs |
| 1e7 | 5.4 ms (1.7× n, 8.8× b) | 9.31 ms (5.1× b) | 47.8 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.75 µs (1.0× n, 1.2× b) | 3.88 µs (1.2× b) | 4.67 µs |
| 1e5 | 56.7 µs (1.5× n, 8.3× b) | 85.9 µs (5.5× b) | 470 µs |
| 1e7 | 6.54 ms (1.4× n, 7.4× b) | 9.44 ms (5.1× b) | 48.6 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 21.8 µs (1.4× n, 0.30× b) | 31.5 µs (0.21× b) | 6.62 µs |
| 1e5 | 1.68 ms (1.6× n, 0.46× b) | 2.64 ms (0.29× b) | 766 µs |
| 1e7 | 168 ms (1.6× n, 0.46× b) | 265 ms (0.29× b) | 77.8 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.08 µs (1.6× n, 0.30× b) | 6.46 µs (0.19× b) | 1.21 µs |
| 1e5 | 66.8 µs (4.8× n, 1.6× b) | 322 µs (0.33× b) | 105 µs |
| 1e7 | 8.01 ms (4.0× n, 1.5× b) | 32.3 ms (0.37× b) | 11.9 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.63 µs (1.2× n, 0.47× b) | 6.58 µs (0.40× b) | 2.63 µs |
| 1e5 | 221 µs (1.5× n, 1.1× b) | 333 µs (0.72× b) | 240 µs |
| 1e7 | 22.9 ms (1.5× n, 1.1× b) | 33.8 ms (0.76× b) | 25.5 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.46 µs (2.2× n, 12.0× b) | 12 µs (5.5× b) | 65.5 µs |
| 1e5 | 128 µs (6.0× n, 64.5× b) | 765 µs (10.8× b) | 8.24 ms |
| 1e7 | 12.4 ms (6.2× n, 66.7× b) | 76.6 ms (10.8× b) | 826 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.96 µs (0.87× n, 0.41× b) | 5.21 µs (0.47× b) | 2.46 µs |
| 1e5 | 391 µs (1.5× n, 1.3× b) | 593 µs (0.86× b) | 513 µs |
| 1e7 | 44.3 ms (1.4× n, 1.3× b) | 63.5 ms (0.87× b) | 55.5 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 6.42 µs (1.6× n, 0.38× b) | 10.1 µs (0.24× b) | 2.46 µs |
| 1e5 | 185 µs (3.2× n, 1.2× b) | 587 µs (0.38× b) | 224 µs |
| 1e7 | 19.2 ms (3.1× n, 1.2× b) | 58.9 ms (0.41× b) | 23.9 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.67 µs (1.4× n, 2.4× b) | 5.17 µs (1.7× b) | 8.71 µs |
| 1e5 | 17.1 µs (10.7× n, 61.4× b) | 183 µs (5.7× b) | 1.05 ms |
| 1e7 | 1.39 ms (13.0× n, 70.4× b) | 18 ms (5.4× b) | 97.6 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.5 µs (1.1× n, 0.05× b) | 3.96 µs (0.04× b) | 0.167 µs |
| 1e5 | 20.4 µs (2.8× n, <0.01× b) | 57.1 µs (<0.01× b) | 0.167 µs |
| 1e7 | 1.68 ms (3.2× n, <0.01× b) | 5.35 ms (<0.01× b) | 0.125 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.0× n, 0.08× b) | 1.08 µs (0.08× b) | 0.084 µs |
| 1e5 | 1.04 µs (1.0× n, 0.04× b) | 1.04 µs (0.04× b) | 0.043 µs |
| 1e7 | 1.08 µs (1.1× n, 0.08× b) | 1.17 µs (0.07× b) | 0.083 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.0× n, 0.08× b) | 1.08 µs (0.08× b) | 0.083 µs |
| 1e5 | 1.08 µs (1.0× n, 0.08× b) | 1.08 µs (0.08× b) | 0.084 µs |
| 1e7 | 1.13 µs (0.93× n, 0.04× b) | 1.04 µs (0.04× b) | 0.043 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.4 ms (12.8× n, 81.6× b) | 18 ms (6.4× b) | 114 ms |
| pairwise | 1.3 ms (14.5× n, 94.5× b) | 18.8 ms (6.5× b) | 122 ms |
| compensated | 4.03 ms (5.3× n, 28.4× b) | 21.4 ms (5.4× b) | 114 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 2.47 ms (6.1× n, 48.1× b) | 15.1 ms (7.9× b) | 119 ms | 14.4 ms (8.2× b) |
| pairwise | 2.38 ms (6.7× n, 49.9× b) | 16 ms (7.5× b) | 119 ms | 14.4 ms (8.3× b) |
| compensated | 5.21 ms (3.1× n, 22.7× b) | 16.1 ms (7.4× b) | 119 ms | 14.4 ms (8.2× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 47.7 ms (2.6× n, 2.8× b) | 125 ms (1.1× b) | 133 ms |
| fast | 25.9 ms (4.8× n, 4.9× b) | 125 ms (1.0× b) | 128 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 38.9 ms (0.62× n, 0.62× b) | 24.1 ms (1.0× b) | 24 ms |
| fast | 20.9 ms (1.1× n, 1.1× b) | 23.6 ms (1.0× b) | 24.1 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 62.7 ms (1.3× n, 1.3× b) | 81.2 ms (1.0× b) | 79.6 ms |
| fast | 25.7 ms (3.2× n, 3.1× b) | 81.2 ms (1.0× b) | 79.9 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 64.8 ms (2.5× n, 2.6× b) | 164 ms (1.0× b) | 170 ms |
| fast | 34.1 ms (4.8× n, 5.0× b) | 164 ms (1.0× b) | 170 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 43 ms (1.7× n, 1.3× b) | 75.1 ms (0.72× b) | 54.1 ms |
| fast | 11.5 ms (6.5× n, 4.7× b) | 75.1 ms (0.72× b) | 54.1 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 63 ms (2.0× n, 2.0× b) | 127 ms (1.0× b) | 126 ms |
| fast | 24.7 ms (5.0× n, 5.3× b) | 124 ms (1.1× b) | 131 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 55.4 ms (3.2× n, 3.2× b) | 179 ms (1.0× b) | 179 ms |
| fast | 44.7 ms (4.0× n, 4.0× b) | 179 ms (1.0× b) | 179 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 184 ms (1.5× n, 1.6× b) | 285 ms (1.0× b) | 287 ms |
| fast | 88.2 ms (3.2× n, 3.3× b) | 285 ms (1.0× b) | 291 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 92 ms (3.1× n, 3.0× b) | 284 ms (1.0× b) | 276 ms |
| fast | 83.8 ms (3.4× n, 3.3× b) | 284 ms (1.0× b) | 277 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 184 ms (2.0× n, 2.0× b) | 370 ms (1.0× b) | 376 ms |
| fast | 137 ms (2.7× n, 2.7× b) | 370 ms (1.0× b) | 374 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 295 ms (1.4× n, 1.4× b) | 428 ms (1.0× b) | 417 ms |
| fast | 165 ms (2.6× n, 2.6× b) | 432 ms (1.0× b) | 421 ms |

