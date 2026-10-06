# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-06T03:47:17Z |
| rsimd | 0.1.0 (git 0cba3207c1e74b227b406e081ce68b6279a2a82f-dirty) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `aarch64-linux-gnu-gcc` |
| CFLAGS | `-g -O2 -ffile-prefix-map=/build/reproducible-path/r-base-4.6.1=. -fstack-protector-strong -fstack-clash-protection -Wformat -Werror=format-security -mbranch-protection=standard -Wdate-time -D_FORTIFY_SOURCE=2` |
| OS | Linux 7.0.14-orbstack-00380-ga7e0a2dc9535 (aarch64) |
| CPU | implementer 0x61 part 0x000, 8 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16 |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 406 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.46 µs (1.5× n, 3.0× b) | 5.13 µs (2.0× b) | 10.4 µs |
| 1e5 | 17.4 µs (10.7× n, 65.9× b) | 185 µs (6.2× b) | 1.14 ms |
| 1e7 | 1.44 ms (12.4× n, 80.2× b) | 17.9 ms (6.4× b) | 115 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.46 µs (1.1× n, 0.19× b) | 3.83 µs (0.17× b) | 0.667 µs |
| 1e5 | 20.6 µs (2.8× n, 2.6× b) | 57.3 µs (0.93× b) | 53.1 µs |
| 1e7 | 1.68 ms (3.2× n, 3.2× b) | 5.38 ms (1.0× b) | 5.36 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.42 µs (1.4× n, 6.0× b) | 6.17 µs (4.3× b) | 26.5 µs |
| 1e5 | 18.3 µs (10.1× n, 167.3× b) | 184 µs (16.6× b) | 3.05 ms |
| 1e7 | 1.42 ms (12.7× n, 209.0× b) | 18 ms (16.4× b) | 296 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 6.15 µs (1.2× n, 1.7× b) | 7.38 µs (1.4× b) | 10.5 µs | 1.67 µs (6.3× b) |
| 1e5 | 25.2 µs (6.3× n, 46.3× b) | 157 µs (7.4× b) | 1.16 ms | 143 µs (8.1× b) |
| 1e7 | 2.57 ms (5.9× n, 60.2× b) | 15.2 ms (10.2× b) | 155 ms | 14.6 ms (10.6× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.5 µs (1.2× n, 0.14× b) | 4.21 µs (0.12× b) | 0.501 µs |
| 1e5 | 25.8 µs (3.7× n, 1.1× b) | 95.5 µs (0.29× b) | 27.6 µs |
| 1e7 | 4.32 ms (2.9× n, 1.4× b) | 12.6 ms (0.48× b) | 6.07 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.48 µs (1.2× n, 0.35× b) | 4.21 µs (0.29× b) | 1.21 µs |
| 1e5 | 35.3 µs (3.1× n, 7.0× b) | 110 µs (2.2× b) | 248 µs |
| 1e7 | 4.98 ms (2.3× n, 8.6× b) | 11.7 ms (3.7× b) | 43 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.88 µs (1.3× n, 0.18× b) | 4.96 µs (0.14× b) | 0.708 µs |
| 1e5 | 35.9 µs (3.8× n, 1.6× b) | 136 µs (0.43× b) | 58.8 µs |
| 1e7 | 6.32 ms (2.8× n, 1.3× b) | 17.8 ms (0.47× b) | 8.33 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.42 µs (1.1× n, 1.3× b) | 3.88 µs (1.1× b) | 4.33 µs |
| 1e5 | 25.7 µs (3.3× n, 15.2× b) | 84.2 µs (4.6× b) | 391 µs |
| 1e7 | 6.35 ms (1.9× n, 7.8× b) | 12.2 ms (4.1× b) | 49.3 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.17 µs (1.3× n, 0.39× b) | 6.67 µs (0.30× b) | 2 µs |
| 1e5 | 170 µs (1.9× n, 1.1× b) | 318 µs (0.59× b) | 189 µs |
| 1e7 | 20.7 ms (1.5× n, 1.1× b) | 31.4 ms (0.70× b) | 21.9 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.21 µs (1.2× n, 0.28× b) | 1.46 µs (0.23× b) | 0.334 µs |
| 1e5 | 10 µs (2.8× n, 2.7× b) | 28.4 µs (0.94× b) | 26.7 µs |
| 1e7 | 1.4 ms (2.0× n, 1.9× b) | 2.74 ms (1.0× b) | 2.72 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.13 µs (1.2× n, 0.30× b) | 1.38 µs (0.24× b) | 0.334 µs |
| 1e5 | 5.54 µs (5.1× n, 4.8× b) | 28.3 µs (0.94× b) | 26.6 µs |
| 1e7 | 652 µs (4.2× n, 4.2× b) | 2.71 ms (1.0× b) | 2.71 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.29 µs (1.4× n, 0.20× b) | 3.21 µs (0.14× b) | 0.459 µs |
| 1e5 | 23.9 µs (4.6× n, 1.1× b) | 109 µs (0.25× b) | 27.3 µs |
| 1e7 | 2.96 ms (4.2× n, 1.4× b) | 12.6 ms (0.33× b) | 4.12 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 6.67 µs (1.4× n, 0.13× b) | 9.08 µs (0.09× b) | 0.834 µs |
| 1e5 | 54.6 µs (10.4× n, 1.2× b) | 568 µs (0.12× b) | 66.7 µs |
| 1e7 | 7.03 ms (10.9× n, 1.3× b) | 76.5 ms (0.12× b) | 8.97 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.58 µs (1.1× n, 0.35× b) | 5.08 µs (0.32× b) | 1.63 µs |
| 1e5 | 31.8 µs (2.2× n, 4.4× b) | 69.1 µs (2.0× b) | 139 µs |
| 1e7 | 3.29 ms (2.7× n, 5.7× b) | 8.75 ms (2.2× b) | 18.8 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.5 µs (1.3× n, 0.40× b) | 6.04 µs (0.30× b) | 1.79 µs |
| 1e5 | 18.4 µs (10.3× n, 7.5× b) | 191 µs (0.73× b) | 139 µs |
| 1e7 | 1.76 ms (4.2× n, 11.8× b) | 7.48 ms (2.8× b) | 20.7 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.17 µs (1.3× n, 0.61× b) | 4.13 µs (0.47× b) | 1.94 µs |
| 1e5 | 31 µs (11.3× n, 5.3× b) | 350 µs (0.47× b) | 165 µs |
| 1e7 | 6.95 ms (6.9× n, 4.8× b) | 47.8 ms (0.70× b) | 33.3 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.4 µs (1.4× n, 4.0× b) | 4.88 µs (2.8× b) | 13.4 µs |
| 1e5 | 50.6 µs (8.3× n, 24.9× b) | 419 µs (3.0× b) | 1.26 ms |
| 1e7 | 7.49 ms (7.8× n, 26.9× b) | 58.5 ms (3.4× b) | 202 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.46 µs (1.0× n, 0.16× b) | 3.52 µs (0.15× b) | 0.543 µs |
| 1e5 | 33.1 µs (1.0× n, 0.95× b) | 34.1 µs (0.92× b) | 31.3 µs |
| 1e7 | 27.5 ms (0.87× n, 0.20× b) | 24 ms (0.23× b) | 5.48 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.67 µs (1.1× n, 1.3× b) | 4.04 µs (1.2× b) | 4.71 µs |
| 1e5 | 54.5 µs (1.8× n, 9.5× b) | 96.1 µs (5.4× b) | 516 µs |
| 1e7 | 17.3 ms (0.76× n, 4.4× b) | 13.1 ms (5.8× b) | 75.9 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.63 µs (1.1× n, 1.3× b) | 3.83 µs (1.2× b) | 4.67 µs |
| 1e5 | 59.8 µs (1.5× n, 8.6× b) | 91.5 µs (5.6× b) | 514 µs |
| 1e7 | 11.6 ms (3.0× n, 5.4× b) | 34.6 ms (1.8× b) | 63.2 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 21.8 µs (1.5× n, 0.31× b) | 32 µs (0.21× b) | 6.67 µs |
| 1e5 | 1.72 ms (1.6× n, 0.45× b) | 2.72 ms (0.28× b) | 774 µs |
| 1e7 | 197 ms (1.7× n, 0.62× b) | 331 ms (0.37× b) | 122 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.21 µs (1.5× n, 0.30× b) | 6.5 µs (0.19× b) | 1.25 µs |
| 1e5 | 72.5 µs (4.5× n, 1.5× b) | 326 µs (0.32× b) | 105 µs |
| 1e7 | 73.9 ms (0.62× n, 0.32× b) | 45.7 ms (0.52× b) | 23.8 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.58 µs (1.2× n, 0.49× b) | 6.63 µs (0.41× b) | 2.71 µs |
| 1e5 | 229 µs (1.5× n, 1.1× b) | 342 µs (0.72× b) | 247 µs |
| 1e7 | 39.2 ms (1.6× n, 2.1× b) | 64.3 ms (1.3× b) | 82.9 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.58 µs (2.1× n, 11.7× b) | 11.9 µs (5.5× b) | 65.4 µs |
| 1e5 | 130 µs (5.9× n, 63.3× b) | 767 µs (10.7× b) | 8.21 ms |
| 1e7 | 14.4 ms (5.6× n, 58.1× b) | 80.9 ms (10.3× b) | 835 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 5.88 µs (0.89× n, 0.43× b) | 5.21 µs (0.48× b) | 2.5 µs |
| 1e5 | 390 µs (1.6× n, 1.3× b) | 607 µs (0.85× b) | 517 µs |
| 1e7 | 54 ms (1.6× n, 1.3× b) | 86.7 ms (0.83× b) | 72 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 6.17 µs (1.7× n, 0.40× b) | 10.5 µs (0.24× b) | 2.46 µs |
| 1e5 | 142 µs (4.3× n, 1.8× b) | 611 µs (0.41× b) | 249 µs |
| 1e7 | 34.6 ms (3.2× n, 0.80× b) | 109 ms (0.25× b) | 27.6 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.5 µs (1.5× n, 2.5× b) | 5.21 µs (1.7× b) | 8.75 µs |
| 1e5 | 17.6 µs (10.9× n, 54.1× b) | 192 µs (5.0× b) | 954 µs |
| 1e7 | 1.73 ms (11.1× n, 56.7× b) | 19.3 ms (5.1× b) | 98.3 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.58 µs (1.1× n, 0.05× b) | 3.83 µs (0.04× b) | 0.167 µs |
| 1e5 | 20.5 µs (2.8× n, <0.01× b) | 58.5 µs (<0.01× b) | 0.166 µs |
| 1e7 | 1.8 ms (3.1× n, <0.01× b) | 5.57 ms (<0.01× b) | 0.125 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.0× n, 0.08× b) | 1.13 µs (0.07× b) | 0.084 µs |
| 1e5 | 1.08 µs (1.0× n, <0.01× b) | 1.13 µs (<0.01× b) | 0.001 µs |
| 1e7 | 1.29 µs (1.1× n, 0.10× b) | 1.42 µs (0.09× b) | 0.126 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.25 µs (0.90× n, 0.07× b) | 1.13 µs (0.07× b) | 0.084 µs |
| 1e5 | 1.13 µs (1.0× n, 0.07× b) | 1.13 µs (0.07× b) | 0.084 µs |
| 1e7 | 1.25 µs (1.0× n, 0.07× b) | 1.25 µs (0.07× b) | 0.084 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 2.23 ms (10.0× n, 57.9× b) | 22.3 ms (5.8× b) | 129 ms |
| pairwise | 1.54 ms (12.5× n, 74.7× b) | 19.1 ms (6.0× b) | 115 ms |
| compensated | 4.04 ms (5.3× n, 28.5× b) | 21.4 ms (5.4× b) | 115 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 2.42 ms (6.4× n, 54.1× b) | 15.6 ms (8.4× b) | 131 ms | 15.2 ms (8.6× b) |
| pairwise | 2.83 ms (6.5× n, 47.5× b) | 18.4 ms (7.3× b) | 135 ms | 14.7 ms (9.2× b) |
| compensated | 5.25 ms (3.1× n, 24.1× b) | 16.4 ms (7.7× b) | 126 ms | 14.5 ms (8.7× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 52.3 ms (2.7× n, 2.6× b) | 140 ms (1.0× b) | 136 ms |
| fast | 28.8 ms (5.0× n, 4.5× b) | 144 ms (0.90× b) | 129 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 28.6 ms (1.0× n, 1.5× b) | 28.9 ms (1.5× b) | 43 ms |
| fast | 27.4 ms (1.3× n, 1.2× b) | 36.8 ms (0.87× b) | 32.1 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 69.5 ms (1.3× n, 1.2× b) | 91.6 ms (0.94× b) | 86.1 ms |
| fast | 32.3 ms (3.2× n, 3.5× b) | 104 ms (1.1× b) | 114 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 116 ms (1.5× n, 1.5× b) | 175 ms (1.0× b) | 172 ms |
| fast | 34.6 ms (4.9× n, 5.0× b) | 169 ms (1.0× b) | 175 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 43.6 ms (1.7× n, 1.3× b) | 72.7 ms (0.78× b) | 56.8 ms |
| fast | 11.6 ms (7.0× n, 5.7× b) | 81.6 ms (0.81× b) | 66.2 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 112 ms (0.82× n, 0.82× b) | 92 ms (1.0× b) | 91.1 ms |
| fast | 96.5 ms (1.0× n, 0.88× b) | 99.8 ms (0.85× b) | 85 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 58.5 ms (0.92× n, 0.95× b) | 54 ms (1.0× b) | 55.5 ms |
| fast | 56.4 ms (1.0× n, 0.92× b) | 58 ms (0.89× b) | 51.7 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 63.5 ms (2.2× n, 2.1× b) | 140 ms (1.0× b) | 134 ms |
| fast | 28.1 ms (4.8× n, 5.5× b) | 135 ms (1.1× b) | 154 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 70.6 ms (2.9× n, 2.9× b) | 208 ms (1.0× b) | 203 ms |
| fast | 75 ms (2.7× n, 2.7× b) | 204 ms (1.0× b) | 204 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 193 ms (1.6× n, 1.5× b) | 306 ms (1.0× b) | 298 ms |
| fast | 96.4 ms (3.2× n, 3.4× b) | 310 ms (1.1× b) | 332 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 109 ms (2.9× n, 2.9× b) | 320 ms (1.0× b) | 318 ms |
| fast | 102 ms (3.3× n, 3.0× b) | 334 ms (0.93× b) | 310 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 219 ms (1.8× n, 1.9× b) | 389 ms (1.1× b) | 415 ms |
| fast | 145 ms (2.6× n, 2.7× b) | 371 ms (1.1× b) | 392 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 107 ms (1.2× n, 1.3× b) | 127 ms (1.1× b) | 139 ms |
| fast | 82.1 ms (1.6× n, 1.7× b) | 128 ms (1.1× b) | 137 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 351 ms (1.3× n, 1.3× b) | 456 ms (1.0× b) | 470 ms |
| fast | 198 ms (2.2× n, 2.1× b) | 437 ms (1.0× b) | 420 ms |

