# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-06T23:05:04Z |
| rsimd | 0.1.0 (git 25e78a59ab068c6d022b8c72626c3943b6786302-dirty) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `aarch64-linux-gnu-gcc` |
| CFLAGS | `-g -O2 -ffile-prefix-map=/build/reproducible-path/r-base-4.6.1=. -fstack-protector-strong -fstack-clash-protection -Wformat -Werror=format-security -mbranch-protection=standard -Wdate-time -D_FORTIFY_SOURCE=2` |
| OS | Linux 7.0.14-orbstack-00380-ga7e0a2dc9535 (aarch64) |
| CPU | implementer 0x61 part 0x000, 8 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16 |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 362 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.751 µs (3.3× n, 13.8× b) | 2.5 µs (4.1× b) | 10.4 µs |
| 1e5 | 14 µs (12.8× n, 86.9× b) | 179 µs (6.8× b) | 1.22 ms |
| 1e7 | 1.37 ms (13.5× n, 83.7× b) | 18.5 ms (6.2× b) | 115 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.833 µs (1.4× n, 0.80× b) | 1.17 µs (0.57× b) | 0.667 µs |
| 1e5 | 17.3 µs (3.1× n, 3.1× b) | 54 µs (1.0× b) | 53.2 µs |
| 1e7 | 1.67 ms (3.2× n, 3.2× b) | 5.35 ms (1.0× b) | 5.33 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.791 µs (3.2× n, 32.9× b) | 2.5 µs (10.4× b) | 26 µs |
| 1e5 | 14.1 µs (12.7× n, 217.7× b) | 180 µs (17.1× b) | 3.08 ms |
| 1e7 | 1.36 ms (12.9× n, 214.8× b) | 17.6 ms (16.7× b) | 293 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 1 µs (2.3× n, 10.5× b) | 2.33 µs (4.5× b) | 10.5 µs | 1.67 µs (6.3× b) |
| 1e5 | 19.6 µs (7.7× n, 59.1× b) | 150 µs (7.7× b) | 1.16 ms | 143 µs (8.1× b) |
| 1e7 | 2.37 ms (6.4× n, 50.5× b) | 15.1 ms (7.9× b) | 120 ms | 14.4 ms (8.3× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1 µs (1.7× n, 0.54× b) | 1.75 µs (0.31× b) | 0.542 µs |
| 1e5 | 22.5 µs (4.1× n, 1.2× b) | 91.9 µs (0.30× b) | 27.2 µs |
| 1e7 | 3.57 ms (3.0× n, 1.1× b) | 10.6 ms (0.36× b) | 3.85 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.17 µs (1.6× n, 1.1× b) | 1.88 µs (0.67× b) | 1.25 µs |
| 1e5 | 31.8 µs (3.4× n, 7.7× b) | 107 µs (2.3× b) | 246 µs |
| 1e7 | 3.65 ms (3.1× n, 10.9× b) | 11.2 ms (3.6× b) | 39.8 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.33 µs (1.7× n, 0.53× b) | 2.25 µs (0.32× b) | 0.709 µs |
| 1e5 | 41 µs (3.2× n, 1.4× b) | 133 µs (0.42× b) | 56.4 µs |
| 1e7 | 5.19 ms (2.9× n, 1.3× b) | 15.2 ms (0.46× b) | 6.96 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.29 µs (1.4× n, 3.4× b) | 1.83 µs (2.4× b) | 4.33 µs |
| 1e5 | 23.1 µs (3.2× n, 17.0× b) | 73.1 µs (5.4× b) | 393 µs |
| 1e7 | 3.56 ms (2.2× n, 13.5× b) | 7.91 ms (6.1× b) | 48.3 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.25 µs (1.6× n, 0.89× b) | 3.54 µs (0.56× b) | 2 µs |
| 1e5 | 166 µs (1.9× n, 1.1× b) | 309 µs (0.60× b) | 186 µs |
| 1e7 | 17.7 ms (1.8× n, 1.1× b) | 31.2 ms (0.64× b) | 19.8 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.458 µs (1.5× n, 0.73× b) | 0.667 µs (0.50× b) | 0.334 µs |
| 1e5 | 8.92 µs (3.0× n, 3.0× b) | 27.2 µs (1.0× b) | 26.6 µs |
| 1e7 | 1.24 ms (2.2× n, 2.2× b) | 2.69 ms (1.0× b) | 2.69 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.416 µs (1.7× n, 0.80× b) | 0.708 µs (0.47× b) | 0.334 µs |
| 1e5 | 4.42 µs (6.1× n, 6.0× b) | 27.1 µs (1.0× b) | 26.6 µs |
| 1e7 | 603 µs (4.4× n, 4.4× b) | 2.68 ms (1.0× b) | 2.68 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.834 µs (1.9× n, 0.60× b) | 1.58 µs (0.32× b) | 0.5 µs |
| 1e5 | 23 µs (4.6× n, 1.2× b) | 107 µs (0.25× b) | 27 µs |
| 1e7 | 2.8 ms (4.0× n, 1.2× b) | 11.2 ms (0.29× b) | 3.25 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.2× n, 1.6× b) | 1.29 µs (1.3× b) | 1.69 µs |
| 1e5 | 54 µs (1.5× n, 4.7× b) | 78.4 µs (3.2× b) | 254 µs |
| 1e7 | 5.87 ms (1.4× n, 5.3× b) | 8.02 ms (3.9× b) | 31.3 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.543 µs (1.9× n, 1.2× b) | 1.04 µs (0.64× b) | 0.667 µs |
| 1e5 | 8.17 µs (6.6× n, 6.5× b) | 54 µs (1.0× b) | 53 µs |
| 1e7 | 1.22 ms (4.4× n, 4.4× b) | 5.38 ms (1.0× b) | 5.33 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.917 µs (1.3× n, 1.6× b) | 1.21 µs (1.2× b) | 1.5 µs |
| 1e5 | 30 µs (1.8× n, 4.5× b) | 54.2 µs (2.5× b) | 133 µs |
| 1e7 | 2.95 ms (1.8× n, 6.6× b) | 5.39 ms (3.6× b) | 19.6 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.876 µs (1.3× n, 1.2× b) | 1.17 µs (0.93× b) | 1.08 µs |
| 1e5 | 29.9 µs (1.8× n, 2.7× b) | 54.1 µs (1.5× b) | 80.9 µs |
| 1e7 | 3.1 ms (1.9× n, 2.9× b) | 5.77 ms (1.6× b) | 9.15 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.75 µs (2.3× n, 0.50× b) | 4 µs (0.22× b) | 0.875 µs |
| 1e5 | 48.5 µs (11.6× n, 1.4× b) | 561 µs (0.12× b) | 66.5 µs |
| 1e7 | 5.25 ms (11.5× n, 1.4× b) | 60.3 ms (0.12× b) | 7.31 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.29 µs (1.5× n, 1.3× b) | 1.88 µs (0.89× b) | 1.67 µs |
| 1e5 | 28.1 µs (2.3× n, 4.9× b) | 65.5 µs (2.1× b) | 138 µs |
| 1e7 | 2.72 ms (2.4× n, 5.4× b) | 6.46 ms (2.3× b) | 14.7 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.13 µs (1.4× n, 1.6× b) | 1.58 µs (1.1× b) | 1.79 µs |
| 1e5 | 14.7 µs (3.7× n, 9.4× b) | 54.9 µs (2.5× b) | 138 µs |
| 1e7 | 1.37 ms (4.0× n, 10.6× b) | 5.47 ms (2.7× b) | 14.6 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.876 µs (2.1× n, 2.2× b) | 1.83 µs (1.1× b) | 1.94 µs |
| 1e5 | 27.6 µs (12.4× n, 5.9× b) | 343 µs (0.48× b) | 164 µs |
| 1e7 | 3.24 ms (13.6× n, 5.7× b) | 44.1 ms (0.42× b) | 18.4 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.04 µs (2.4× n, 13.0× b) | 2.54 µs (5.3× b) | 13.5 µs |
| 1e5 | 47.6 µs (10.4× n, 25.8× b) | 495 µs (2.5× b) | 1.23 ms |
| 1e7 | 5.19 ms (11.0× n, 25.9× b) | 57 ms (2.4× b) | 134 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.918 µs (1.0× n, 0.59× b) | 0.918 µs (0.59× b) | 0.542 µs |
| 1e5 | 27.8 µs (1.0× n, 1.0× b) | 28.6 µs (0.94× b) | 27 µs |
| 1e7 | 3.7 ms (1.0× n, 1.0× b) | 3.76 ms (1.0× b) | 3.68 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.04 µs (1.4× n, 4.5× b) | 1.46 µs (3.2× b) | 4.67 µs |
| 1e5 | 45.2 µs (1.8× n, 10.4× b) | 83.3 µs (5.6× b) | 469 µs |
| 1e7 | 5.34 ms (1.7× n, 10.1× b) | 9.26 ms (5.8× b) | 53.9 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.17 µs (1.3× n, 4.0× b) | 1.5 µs (3.1× b) | 4.67 µs |
| 1e5 | 54.1 µs (1.5× n, 8.7× b) | 83.3 µs (5.7× b) | 471 µs |
| 1e7 | 6.51 ms (1.4× n, 7.3× b) | 9.17 ms (5.2× b) | 47.4 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 20.1 µs (1.5× n, 0.33× b) | 30.5 µs (0.22× b) | 6.62 µs |
| 1e5 | 1.68 ms (1.6× n, 0.46× b) | 2.67 ms (0.29× b) | 766 µs |
| 1e7 | 168 ms (1.6× n, 0.46× b) | 268 ms (0.29× b) | 77.7 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.92 µs (1.8× n, 0.41× b) | 5.38 µs (0.22× b) | 1.21 µs |
| 1e5 | 65.7 µs (5.2× n, 1.6× b) | 340 µs (0.30× b) | 102 µs |
| 1e7 | 7.93 ms (4.3× n, 1.6× b) | 34.1 ms (0.36× b) | 12.3 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.25 µs (1.3× n, 0.62× b) | 5.5 µs (0.48× b) | 2.63 µs |
| 1e5 | 223 µs (1.5× n, 1.1× b) | 338 µs (0.74× b) | 250 µs |
| 1e7 | 23.8 ms (1.5× n, 1.1× b) | 34.5 ms (0.74× b) | 25.7 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.92 µs (4.3× n, 34.0× b) | 8.21 µs (8.0× b) | 65.3 µs |
| 1e5 | 124 µs (6.2× n, 68.7× b) | 760 µs (11.2× b) | 8.48 ms |
| 1e7 | 12.4 ms (6.2× n, 66.7× b) | 76.1 ms (10.8× b) | 825 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.42 µs (0.83× n, 0.56× b) | 3.67 µs (0.67× b) | 2.46 µs |
| 1e5 | 396 µs (1.5× n, 1.3× b) | 610 µs (0.87× b) | 533 µs |
| 1e7 | 44.4 ms (1.5× n, 1.3× b) | 66.7 ms (0.84× b) | 55.9 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.04 µs (2.4× n, 0.81× b) | 7.25 µs (0.34× b) | 2.46 µs |
| 1e5 | 151 µs (3.9× n, 1.5× b) | 591 µs (0.38× b) | 227 µs |
| 1e7 | 15.4 ms (3.8× n, 1.6× b) | 58.7 ms (0.41× b) | 23.9 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.791 µs (3.1× n, 11.0× b) | 2.46 µs (3.5× b) | 8.71 µs |
| 1e5 | 14.3 µs (12.7× n, 64.8× b) | 182 µs (5.1× b) | 926 µs |
| 1e7 | 1.4 ms (13.0× n, 69.5× b) | 18.3 ms (5.3× b) | 97.4 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.834 µs (1.4× n, 0.20× b) | 1.17 µs (0.14× b) | 0.167 µs |
| 1e5 | 17.3 µs (3.1× n, <0.01× b) | 53.7 µs (<0.01× b) | 0.168 µs |
| 1e7 | 1.73 ms (3.2× n, <0.01× b) | 5.64 ms (<0.01× b) | 0.125 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.376 µs (1.1× n, 0.22× b) | 0.418 µs (0.20× b) | 0.084 µs |
| 1e5 | 0.375 µs (1.0× n, 0.11× b) | 0.376 µs (0.11× b) | 0.042 µs |
| 1e7 | 0.416 µs (1.1× n, 0.30× b) | 0.459 µs (0.27× b) | 0.124 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.376 µs (1.1× n, 0.22× b) | 0.417 µs (0.20× b) | 0.084 µs |
| 1e5 | 0.375 µs (1.1× n, 0.22× b) | 0.417 µs (0.20× b) | 0.084 µs |
| 1e7 | 0.376 µs (1.1× n, 0.11× b) | 0.417 µs (0.10× b) | 0.042 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.46 ms (13.0× n, 81.3× b) | 18.9 ms (6.3× b) | 118 ms |
| pairwise | 1.31 ms (14.4× n, 88.5× b) | 18.8 ms (6.2× b) | 116 ms |
| compensated | 4.17 ms (5.3× n, 27.8× b) | 22.2 ms (5.2× b) | 116 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 2.58 ms (5.9× n, 49.7× b) | 15.3 ms (8.4× b) | 128 ms | 14.4 ms (8.9× b) |
| pairwise | 2.44 ms (6.6× n, 48.9× b) | 16.1 ms (7.4× b) | 119 ms | 14.6 ms (8.2× b) |
| compensated | 5.26 ms (3.1× n, 22.6× b) | 16.2 ms (7.3× b) | 119 ms | 14.7 ms (8.1× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 49.1 ms (2.6× n, 2.6× b) | 130 ms (1.0× b) | 130 ms |
| fast | 32.1 ms (4.1× n, 4.0× b) | 132 ms (1.0× b) | 128 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 26.2 ms (0.92× n, 0.92× b) | 24.2 ms (1.0× b) | 24.1 ms |
| fast | 21.7 ms (1.2× n, 1.1× b) | 25 ms (0.93× b) | 23.4 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 64.9 ms (1.3× n, 1.3× b) | 85.4 ms (1.0× b) | 81.7 ms |
| fast | 26.9 ms (3.1× n, 3.0× b) | 84.7 ms (1.0× b) | 80.5 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 65.3 ms (2.5× n, 2.6× b) | 166 ms (1.0× b) | 173 ms |
| fast | 33.9 ms (4.9× n, 5.1× b) | 167 ms (1.0× b) | 173 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 43.2 ms (1.7× n, 1.3× b) | 73.5 ms (0.77× b) | 56.8 ms |
| fast | 11.8 ms (6.3× n, 4.8× b) | 73.9 ms (0.77× b) | 56.7 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 104 ms (0.93× n, 0.76× b) | 97.2 ms (0.82× b) | 79.6 ms |
| fast | 101 ms (1.0× n, 0.81× b) | 98.2 ms (0.83× b) | 81.6 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 54.8 ms (1.0× n, 0.89× b) | 52.7 ms (0.92× b) | 48.7 ms |
| fast | 52 ms (1.0× n, 0.93× b) | 51.7 ms (0.94× b) | 48.3 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 63.7 ms (2.0× n, 2.0× b) | 125 ms (1.0× b) | 128 ms |
| fast | 24.6 ms (5.1× n, 5.2× b) | 125 ms (1.0× b) | 127 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 56.1 ms (3.2× n, 3.2× b) | 182 ms (1.0× b) | 181 ms |
| fast | 46.3 ms (3.9× n, 3.9× b) | 180 ms (1.0× b) | 180 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 185 ms (1.6× n, 1.6× b) | 295 ms (1.0× b) | 293 ms |
| fast | 89.5 ms (3.2× n, 3.3× b) | 289 ms (1.0× b) | 291 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 94.8 ms (3.0× n, 3.0× b) | 289 ms (1.0× b) | 281 ms |
| fast | 86.1 ms (3.4× n, 3.2× b) | 289 ms (1.0× b) | 278 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 184 ms (2.0× n, 2.0× b) | 369 ms (1.0× b) | 375 ms |
| fast | 137 ms (2.7× n, 2.7× b) | 368 ms (1.0× b) | 375 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 86.8 ms (1.3× n, 1.3× b) | 109 ms (1.1× b) | 115 ms |
| fast | 67.3 ms (1.6× n, 1.7× b) | 110 ms (1.0× b) | 115 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 288 ms (1.5× n, 1.5× b) | 426 ms (1.0× b) | 421 ms |
| fast | 167 ms (2.6× n, 2.5× b) | 430 ms (1.0× b) | 418 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.125 µs | 0.625 µs (0.20× b) | 0.458 µs (0.27× b) |
| 1e1 | 0.209 µs | 0.625 µs (0.33× b) | 0.458 µs (0.46× b) |
| 1e2 | 1.13 µs | 0.626 µs (1.8× b) | 0.459 µs (2.5× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.083 µs | 0.751 µs (0.11× b) |
| 1e1 | 0.084 µs | 0.709 µs (0.12× b) |
| 1e2 | 0.126 µs | 0.751 µs (0.17× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.084 µs | 0.5 µs (0.17× b) |
| 1e1 | 0.084 µs | 0.5 µs (0.17× b) |
| 1e2 | 0.25 µs | 0.667 µs (0.37× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.043 µs | 0.667 µs (0.06× b) |
| 1e1 | 0.125 µs | 0.667 µs (0.19× b) |
| 1e2 | 0.126 µs | 0.708 µs (0.18× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.083 µs | 0.75 µs (0.11× b) |
| 1e1 | 0.084 µs | 0.792 µs (0.11× b) |
| 1e2 | 0.125 µs | 0.793 µs (0.16× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.125 µs | 0.833 µs (0.15× b) |
| 1e1 | 0.292 µs | 0.834 µs (0.35× b) |
| 1e2 | 1.25 µs | 0.875 µs (1.4× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.084 µs | 2.08 µs (0.04× b) |
| 1e1 | 0.084 µs | 2.12 µs (0.04× b) |
| 1e2 | 0.208 µs | 2.21 µs (0.09× b) |

