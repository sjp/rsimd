# rsimd benchmark results

| | |
|---|---|
| Timestamp (UTC) | 2026-10-07T00:53:27Z |
| rsimd | 0.1.0 (git ce52d945d5ac8f136f10153f5ec32175fef5639d-dirty) |
| R | R version 4.6.1 (2026-06-24) |
| Compiler | `aarch64-linux-gnu-gcc` |
| CFLAGS | `-g -O2 -ffile-prefix-map=/build/reproducible-path/r-base-4.6.1=. -fstack-protector-strong -fstack-clash-protection -Wformat -Werror=format-security -mbranch-protection=standard -Wdate-time -D_FORTIFY_SOURCE=2` |
| OS | Linux 7.0.14-orbstack-00380-ga7e0a2dc9535 (aarch64) |
| CPU | implementer 0x61 part 0x000, 8 cores |
| CPU features | neon, fp16, dotprod, i8mm, bf16 |
| Tiers | neon, none (auto: neon) |
| bench | 1.1.4 |
| Run time | 374 s |

Cells show the median time per call, then in parentheses the speedup versus the `none` tier (`n`) and versus base R (`b`); above 1× is faster. ⚠ marks a SIMD tier less than 1.1× faster than `none` on a compute-bound op at n ≥ 1e5, which suggests it is running scalar code. Timings are informational: they vary between machines and runs, and at n = 1e3 the fixed per-call overhead dominates.

## Main table (precision "fast", no NAs)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.792 µs (3.1× n, 13.1× b) | 2.42 µs (4.3× b) | 10.4 µs |
| 1e5 | 14.1 µs (12.8× n, 86.6× b) | 180 µs (6.8× b) | 1.22 ms |
| 1e7 | 1.36 ms (13.7× n, 84.2× b) | 18.7 ms (6.1× b) | 115 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.708 µs (1.6× n, 0.88× b) | 1.17 µs (0.54× b) | 0.626 µs |
| 1e5 | 10.6 µs (5.1× n, 5.0× b) | 54 µs (1.0× b) | 53.1 µs |
| 1e7 | 979 µs (5.5× n, 5.4× b) | 5.36 ms (1.0× b) | 5.33 ms |

### mean (double)

`simd_mean(x)` vs base `mean(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.792 µs (3.1× n, 35.0× b) | 2.46 µs (11.3× b) | 27.7 µs |
| 1e5 | 14 µs (12.8× n, 219.5× b) | 180 µs (17.1× b) | 3.07 ms |
| 1e7 | 1.35 ms (13.3× n, 216.7× b) | 18 ms (16.3× b) | 294 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| n | neon | none | base | blas |
|---|---|---|---|---|
| 1e3 | 1 µs (2.3× n, 10.6× b) | 2.33 µs (4.5× b) | 10.6 µs | 1.67 µs (6.4× b) |
| 1e5 | 19.4 µs (7.8× n, 64.2× b) | 151 µs (8.2× b) | 1.24 ms | 143 µs (8.7× b) |
| 1e7 | 2.37 ms (6.4× n, 50.9× b) | 15.1 ms (8.0× b) | 121 ms | 14.4 ms (8.4× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.04 µs (1.7× n, 0.48× b) | 1.79 µs (0.28× b) | 0.5 µs |
| 1e5 | 22.4 µs (4.1× n, 1.2× b) | 92.5 µs (0.30× b) | 27.3 µs |
| 1e7 | 3.54 ms (3.0× n, 1.1× b) | 10.8 ms (0.36× b) | 3.84 ms |

### add (integer)

`simd_add(xi, yi)` vs base `xi + yi`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.17 µs (1.6× n, 1.0× b) | 1.88 µs (0.64× b) | 1.21 µs |
| 1e5 | 31.8 µs (3.4× n, 7.8× b) | 107 µs (2.3× b) | 249 µs |
| 1e7 | 3.66 ms (3.4× n, 11.1× b) | 12.3 ms (3.3× b) | 40.7 ms |

### fma (double)

`simd_fma(x, y, z)` vs base `x * y + z`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.25 µs (1.8× n, 0.57× b) | 2.29 µs (0.31× b) | 0.709 µs |
| 1e5 | 29.6 µs (4.4× n, 1.9× b) | 131 µs (0.43× b) | 56.1 µs |
| 1e7 | 4.71 ms (3.2× n, 1.5× b) | 15.2 ms (0.45× b) | 6.91 ms |

### pmax (double)

`simd_pmax(x, y)` vs base `pmax(x, y)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.4 µs (1.3× n, 3.1× b) | 1.83 µs (2.3× b) | 4.29 µs |
| 1e5 | 22.9 µs (3.5× n, 17.0× b) | 81.2 µs (4.8× b) | 388 µs |
| 1e7 | 3.51 ms (2.2× n, 13.8× b) | 7.79 ms (6.2× b) | 48.5 ms |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.21 µs (1.6× n, 0.89× b) | 3.5 µs (0.56× b) | 1.96 µs |
| 1e5 | 168 µs (1.8× n, 1.1× b) | 299 µs (0.62× b) | 185 µs |
| 1e7 | 17.7 ms (1.7× n, 1.1× b) | 30.2 ms (0.65× b) | 19.8 ms |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.459 µs (1.5× n, 0.73× b) | 0.667 µs (0.50× b) | 0.333 µs |
| 1e5 | 8.88 µs (3.1× n, 3.0× b) | 27.2 µs (1.0× b) | 26.6 µs |
| 1e7 | 1.23 ms (2.2× n, 2.2× b) | 2.69 ms (1.0× b) | 2.69 ms |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.416 µs (1.6× n, 0.80× b) | 0.667 µs (0.50× b) | 0.334 µs |
| 1e5 | 4.37 µs (6.2× n, 6.1× b) | 27.1 µs (1.0× b) | 26.6 µs |
| 1e7 | 611 µs (4.4× n, 4.4× b) | 2.68 ms (1.0× b) | 2.68 ms |

### is_na (double)

`simd_is_na(x)` vs base `is.na(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.834 µs (1.9× n, 0.60× b) | 1.58 µs (0.32× b) | 0.5 µs |
| 1e5 | 27.5 µs (3.9× n, 1.0× b) | 107 µs (0.25× b) | 27 µs |
| 1e7 | 3.28 ms (3.4× n, 1.0× b) | 11.1 ms (0.29× b) | 3.21 ms |

### which (logical)

`simd_which(xl)` vs base `which(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.08 µs (1.2× n, 1.5× b) | 1.29 µs (1.3× b) | 1.67 µs |
| 1e5 | 53.8 µs (1.5× n, 4.7× b) | 78.3 µs (3.2× b) | 253 µs |
| 1e7 | 5.99 ms (1.3× n, 5.2× b) | 8.01 ms (3.9× b) | 31.2 ms |

### count (logical)

`simd_count(xl)` vs base `sum(xl)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.542 µs (2.0× n, 1.2× b) | 1.08 µs (0.60× b) | 0.646 µs |
| 1e5 | 8.17 µs (6.6× n, 6.5× b) | 54.3 µs (1.0× b) | 53.1 µs |
| 1e7 | 1.22 ms (4.4× n, 4.4× b) | 5.37 ms (1.0× b) | 5.34 ms |

### max_abs (double)

`simd_max_abs(x)` vs base `max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.876 µs (1.3× n, 1.7× b) | 1.17 µs (1.3× b) | 1.5 µs |
| 1e5 | 20.7 µs (2.6× n, 6.4× b) | 54.2 µs (2.5× b) | 133 µs |
| 1e7 | 2.01 ms (2.7× n, 7.2× b) | 5.39 ms (2.7× b) | 14.5 ms |

### which_max_abs (double)

`simd_which_max_abs(x)` vs base `which.max(abs(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.708 µs (1.6× n, 1.5× b) | 1.13 µs (1.0× b) | 1.08 µs |
| 1e5 | 20.7 µs (2.6× n, 3.9× b) | 54.2 µs (1.5× b) | 81.3 µs |
| 1e7 | 2.16 ms (2.7× n, 4.2× b) | 5.75 ms (1.6× b) | 9.11 ms |

### as_integer (double)

`simd_as_integer(xc, mode = "truncating")` vs base `as.integer(xc)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.79 µs (11.3× n, 0.49× b) | 20.3 µs (0.04× b) | 0.874 µs |
| 1e5 | 48.7 µs (41.0× n, 1.4× b) | 2 ms (0.03× b) | 66.6 µs |
| 1e7 | 5.22 ms (38.1× n, 1.4× b) | 199 ms (0.04× b) | 7.23 ms |

### hamming (double)

`simd_hamming(x, xh)` vs base `sum(x != xh)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.33 µs (1.3× n, 1.2× b) | 1.75 µs (1.0× b) | 1.67 µs |
| 1e5 | 28.3 µs (2.3× n, 4.9× b) | 65.4 µs (2.1× b) | 138 µs |
| 1e7 | 2.72 ms (2.4× n, 5.4× b) | 6.47 ms (2.3× b) | 14.6 ms |

### hamming (integer)

`simd_hamming(xi, xih)` vs base `sum(xi != xih)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.17 µs (2.5× n, 1.5× b) | 2.88 µs (0.61× b) | 1.75 µs |
| 1e5 | 14.7 µs (12.7× n, 9.4× b) | 187 µs (0.74× b) | 138 µs |
| 1e7 | 1.37 ms (4.0× n, 10.7× b) | 5.47 ms (2.7× b) | 14.6 ms |

### is_whole (double)

`simd_is_whole(xw)` vs base `xw == trunc(xw)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.917 µs (2.0× n, 2.0× b) | 1.79 µs (1.0× b) | 1.88 µs |
| 1e5 | 30.5 µs (11.8× n, 5.4× b) | 362 µs (0.46× b) | 165 µs |
| 1e7 | 3.49 ms (12.6× n, 5.2× b) | 44.1 ms (0.41× b) | 18.2 ms |

### is_pow2 (double)

`simd_is_pow2(xw)` vs base `xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.04 µs (2.4× n, 12.9× b) | 2.54 µs (5.3× b) | 13.4 µs |
| 1e5 | 47.6 µs (10.3× n, 26.0× b) | 491 µs (2.5× b) | 1.24 ms |
| 1e7 | 5.17 ms (10.9× n, 25.8× b) | 56.4 ms (2.4× b) | 133 ms |

### recip_approx (double)

`simd_recip_approx(xp)` vs base `1/xp`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.917 µs (1.0× n, 0.50× b) | 0.959 µs (0.48× b) | 0.459 µs |
| 1e5 | 28.1 µs (1.0× n, 1.0× b) | 28.8 µs (0.94× b) | 27 µs |
| 1e7 | 3.67 ms (1.0× n, 1.0× b) | 3.74 ms (1.0× b) | 3.7 ms |

### rsqrt (double)

`simd_rsqrt(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.04 µs (1.4× n, 4.7× b) | 1.5 µs (3.3× b) | 4.92 µs |
| 1e5 | 45.3 µs (1.8× n, 10.4× b) | 83.6 µs (5.6× b) | 469 µs |
| 1e7 | 5.37 ms (1.7× n, 8.8× b) | 9.23 ms (5.1× b) | 47.4 ms |

### rsqrt_approx (double)

`simd_rsqrt_approx(xp)` vs base `1/sqrt(xp)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.12 µs (1.3× n, 4.4× b) | 1.5 µs (3.3× b) | 4.92 µs |
| 1e5 | 45.4 µs (1.8× n, 10.3× b) | 83.6 µs (5.6× b) | 466 µs |
| 1e7 | 5.42 ms (1.7× n, 9.8× b) | 9.2 ms (5.8× b) | 53.1 ms |

### rootn (double)

`simd_rootn(xp, 3L)` vs base `xp^(1/3)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 20.1 µs (1.5× n, 0.35× b) | 30 µs (0.23× b) | 7.04 µs |
| 1e5 | 1.68 ms (1.6× n, 0.45× b) | 2.69 ms (0.28× b) | 765 µs |
| 1e7 | 168 ms (1.6× n, 0.46× b) | 268 ms (0.29× b) | 77.8 ms |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2.88 µs (1.8× n, 0.41× b) | 5.21 µs (0.22× b) | 1.17 µs |
| 1e5 | 79.3 µs (4.0× n, 1.3× b) | 320 µs (0.32× b) | 102 µs |
| 1e7 | 7.91 ms (4.1× n, 1.5× b) | 32.5 ms (0.36× b) | 11.8 ms |

### div (complex)

`simd_div(cx, cy)` vs base `cx/cy`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.58 µs (1.2× n, 0.58× b) | 5.38 µs (0.50× b) | 2.67 µs |
| 1e5 | 222 µs (1.5× n, 1.1× b) | 330 µs (0.72× b) | 239 µs |
| 1e7 | 23 ms (1.5× n, 1.1× b) | 33.9 ms (0.75× b) | 25.5 ms |

### prod (complex)

`simd_prod(cu)` vs base `prod(cu)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 2 µs (4.1× n, 33.3× b) | 8.29 µs (8.0× b) | 66.6 µs |
| 1e5 | 124 µs (6.2× n, 66.5× b) | 770 µs (10.7× b) | 8.21 ms |
| 1e7 | 12.4 ms (6.2× n, 66.7× b) | 76.1 ms (10.8× b) | 825 ms |

### abs (complex)

`simd_abs(cx)` vs base `Mod(cx)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 4.33 µs (0.88× n, 0.58× b) | 3.81 µs (0.66× b) | 2.5 µs |
| 1e5 | 386 µs (1.5× n, 1.4× b) | 588 µs (0.89× b) | 523 µs |
| 1e7 | 44 ms (1.4× n, 1.3× b) | 63.6 ms (0.88× b) | 55.7 ms |

### pow_int (complex)

`simd_pow(cs, 3)` vs base `cs^3`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 3.08 µs (2.4× n, 0.78× b) | 7.29 µs (0.33× b) | 2.42 µs |
| 1e5 | 141 µs (4.1× n, 1.6× b) | 572 µs (0.39× b) | 226 µs |
| 1e7 | 15.1 ms (3.8× n, 1.6× b) | 58.2 ms (0.41× b) | 23.9 ms |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.876 µs (2.8× n, 11.9× b) | 2.46 µs (4.2× b) | 10.4 µs |
| 1e5 | 20.7 µs (8.7× n, 58.9× b) | 180 µs (6.8× b) | 1.22 ms |
| 1e7 | 2.01 ms (9.3× n, 56.8× b) | 18.7 ms (6.1× b) | 115 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.791 µs (1.5× n, 1.5× b) | 1.21 µs (1.0× b) | 1.17 µs |
| 1e5 | 19.4 µs (2.8× n, 5.5× b) | 54.2 µs (2.0× b) | 106 µs |
| 1e7 | 1.86 ms (2.9× n, 5.7× b) | 5.4 ms (2.0× b) | 10.7 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.709 µs (1.8× n, 0.94× b) | 1.25 µs (0.53× b) | 0.666 µs |
| 1e5 | 7.46 µs (8.7× n, 7.4× b) | 65.2 µs (0.84× b) | 55 µs |
| 1e7 | 690 µs (9.5× n, 8.0× b) | 6.55 ms (0.84× b) | 5.49 ms |

### is_finite_all (double)

`simd_is_finite_all(x)` vs base `all(is.finite(x))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.5 µs (3.9× n, 2.6× b) | 1.96 µs (0.66× b) | 1.29 µs |
| 1e5 | 11.4 µs (14.0× n, 8.8× b) | 159 µs (0.63× b) | 100 µs |
| 1e7 | 1.28 ms (12.6× n, 10.3× b) | 16 ms (0.82× b) | 13.2 ms |

### bit_and (raw)

`simd_bit_and(rx, ry)` vs base `rx & ry`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 1.13 µs (1.6× n, 0.85× b) | 1.83 µs (0.52× b) | 0.958 µs |
| 1e5 | 4.08 µs (19.8× n, 20.4× b) | 80.7 µs (1.0× b) | 83.3 µs |
| 1e7 | 400 µs (20.3× n, 21.3× b) | 8.1 ms (1.1× b) | 8.53 ms |

### popcount_total (raw)

`simd_popcount_total(rx)` vs base `sum(as.integer(rawToBits(rx)))`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.792 µs (1.4× n, 13.7× b) | 1.13 µs (9.7× b) | 10.9 µs |
| 1e5 | 5.79 µs (5.3× n, 174.0× b) | 30.8 µs (32.8× b) | 1.01 ms |
| 1e7 | 500 µs (5.9× n, 209.4× b) | 2.94 ms (35.6× b) | 105 ms |

## Inputs with 1% NA (precision "fast")

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.791 µs (3.2× n, 9.7× b) | 2.54 µs (3.0× b) | 7.67 µs |
| 1e5 | 1.21 µs (6.8× n, 866.2× b) | 8.25 µs (126.9× b) | 1.05 ms |
| 1e7 | 1.21 µs (6.8× n, 80059.4× b) | 8.25 µs (11730.9× b) | 96.8 ms |

### sum (integer)

`simd_sum(xi)` vs base `sum(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.75 µs (0.89× n, 0.22× b) | 0.667 µs (0.25× b) | 0.166 µs |
| 1e5 | 1.04 µs (0.64× n, 0.12× b) | 0.666 µs (0.19× b) | 0.125 µs |
| 1e7 | 1.04 µs (0.64× n, 0.16× b) | 0.666 µs (0.25× b) | 0.166 µs |

### any_na (double)

`simd_any_na(x)` vs base `anyNA(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.416 µs (1.1× n, 0.20× b) | 0.459 µs (0.18× b) | 0.084 µs |
| 1e5 | 0.416 µs (1.0× n, 0.20× b) | 0.417 µs (0.20× b) | 0.083 µs |
| 1e7 | 0.416 µs (1.1× n, 0.20× b) | 0.459 µs (0.18× b) | 0.084 µs |

### any_na (integer)

`simd_any_na(xi)` vs base `anyNA(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.376 µs (1.1× n, 0.22× b) | 0.416 µs (0.20× b) | 0.084 µs |
| 1e5 | 0.416 µs (1.0× n, 0.20× b) | 0.417 µs (0.20× b) | 0.084 µs |
| 1e7 | 0.375 µs (1.1× n, 0.22× b) | 0.417 µs (0.20× b) | 0.084 µs |

### sum_narm (double)

`simd_sum(x, na.rm = TRUE)` vs base `sum(x, na.rm = TRUE)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.916 µs (2.8× n, 9.5× b) | 2.58 µs (3.4× b) | 8.71 µs |
| 1e5 | 20.7 µs (9.0× n, 58.2× b) | 186 µs (6.5× b) | 1.2 ms |
| 1e7 | 2.02 ms (9.2× n, 56.4× b) | 18.6 ms (6.1× b) | 114 ms |

### min (double)

`simd_min(x)` vs base `min(x)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.833 µs (1.5× n, 1.4× b) | 1.21 µs (1.0× b) | 1.17 µs |
| 1e5 | 1.37 µs (2.2× n, 76.6× b) | 3.08 µs (34.1× b) | 105 µs |
| 1e7 | 1.38 µs (2.2× n, 7674.6× b) | 3.04 µs (3470.1× b) | 10.6 ms |

### min (integer)

`simd_min(xi)` vs base `min(xi)`

| n | neon | none | base |
|---|---|---|---|
| 1e3 | 0.667 µs (1.9× n, 0.19× b) | 1.29 µs (0.10× b) | 0.126 µs |
| 1e5 | 0.875 µs (3.9× n, 0.14× b) | 3.42 µs (0.04× b) | 0.126 µs |
| 1e7 | 1.04 µs (3.3× n, 0.16× b) | 3.42 µs (0.05× b) | 0.166 µs |

## Precision modes (n = 1e7)

### sum (double)

`simd_sum(x)` vs base `sum(x)`

| mode | neon | none | base |
|---|---|---|---|
| fast | 1.35 ms (13.8× n, 84.4× b) | 18.7 ms (6.1× b) | 114 ms |
| pairwise | 1.29 ms (14.6× n, 88.5× b) | 18.9 ms (6.1× b) | 114 ms |
| compensated | 4.02 ms (5.3× n, 28.4× b) | 21.3 ms (5.4× b) | 114 ms |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`; blas: `crossprod(x, y)[1L]`

| mode | neon | none | base | blas |
|---|---|---|---|---|
| fast | 2.4 ms (6.3× n, 49.7× b) | 15.1 ms (7.9× b) | 119 ms | 14.4 ms (8.3× b) |
| pairwise | 2.44 ms (6.6× n, 48.6× b) | 16.1 ms (7.4× b) | 119 ms | 14.4 ms (8.3× b) |
| compensated | 5.2 ms (3.3× n, 22.8× b) | 17.4 ms (6.8× b) | 118 ms | 14.4 ms (8.2× b) |

## Math accuracy modes (n = 1e7)

### sin (double)

`simd_sin(x)` vs base `sin(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 48.5 ms (2.6× n, 2.6× b) | 127 ms (1.0× b) | 127 ms |
| fast | 27.4 ms (4.6× n, 4.7× b) | 127 ms (1.0× b) | 128 ms |

### log (double)

`simd_log(xp)` vs base `log(xp)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 24.8 ms (1.0× n, 0.93× b) | 25.3 ms (0.91× b) | 23.1 ms |
| fast | 20.6 ms (1.2× n, 1.1× b) | 25.3 ms (0.93× b) | 23.5 ms |

### tanh (double)

`simd_tanh(xt)` vs base `tanh(xt)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 69.9 ms (1.2× n, 1.1× b) | 83.4 ms (1.0× b) | 79.9 ms |
| fast | 26.6 ms (3.1× n, 3.0× b) | 83.2 ms (1.0× b) | 80 ms |

### atan2 (double)

`simd_atan2(y, x)` vs base `atan2(y, x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 64.1 ms (2.6× n, 2.7× b) | 164 ms (1.0× b) | 170 ms |
| fast | 35.5 ms (4.6× n, 4.8× b) | 164 ms (1.0× b) | 170 ms |

### hypot (double)

`simd_hypot(x, y)` vs base `sqrt(x * x + y * y)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 42.9 ms (1.7× n, 1.3× b) | 72.4 ms (0.74× b) | 53.8 ms |
| fast | 11.4 ms (6.3× n, 4.7× b) | 72.3 ms (0.74× b) | 53.8 ms |

### pow (double)

`simd_pow(xp, xt)` vs base `xp^xt`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 95.1 ms (1.0× n, 0.83× b) | 96.5 ms (0.82× b) | 79.2 ms |
| fast | 95.2 ms (1.0× n, 0.83× b) | 95.3 ms (0.83× b) | 79.3 ms |

### asinh (double)

`simd_asinh(x)` vs base `asinh(x)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 49.2 ms (1.0× n, 1.0× b) | 49.2 ms (1.0× b) | 47.6 ms |
| fast | 49.2 ms (1.0× n, 1.0× b) | 49.1 ms (1.0× b) | 47.5 ms |

### sqrt (complex)

`simd_sqrt(cs)` vs base `sqrt(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 63.6 ms (1.9× n, 2.0× b) | 124 ms (1.0× b) | 125 ms |
| fast | 24.5 ms (5.0× n, 5.1× b) | 123 ms (1.0× b) | 125 ms |

### exp (complex)

`simd_exp(cs)` vs base `exp(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 56.5 ms (3.1× n, 3.2× b) | 177 ms (1.0× b) | 180 ms |
| fast | 45.9 ms (3.9× n, 3.9× b) | 178 ms (1.0× b) | 179 ms |

### log (complex)

`simd_log(cs)` vs base `log(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 185 ms (1.5× n, 1.5× b) | 284 ms (1.0× b) | 286 ms |
| fast | 88.7 ms (3.2× n, 3.2× b) | 284 ms (1.0× b) | 286 ms |

### sin (complex)

`simd_sin(cs)` vs base `sin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 94.1 ms (3.0× n, 2.9× b) | 283 ms (1.0× b) | 275 ms |
| fast | 85.3 ms (3.3× n, 3.3× b) | 283 ms (1.0× b) | 278 ms |

### asin (complex)

`simd_asin(cs)` vs base `asin(cs)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 183 ms (2.0× n, 2.0× b) | 368 ms (1.0× b) | 374 ms |
| fast | 137 ms (2.7× n, 2.7× b) | 368 ms (1.0× b) | 374 ms |

### asin_cut (complex)

`simd_asin(ccut)` vs base `asin(ccut)`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 86.8 ms (1.3× n, 1.3× b) | 110 ms (1.0× b) | 115 ms |
| fast | 67.2 ms (1.6× n, 1.7× b) | 110 ms (1.0× b) | 114 ms |

### pow (complex)

`simd_pow(cs, 0.5 + (0+0.5i))` vs base `cs^(0.5 + (0+0.5i))`

| mode | neon | none | base |
|---|---|---|---|
| accurate | 289 ms (1.5× n, 1.4× b) | 424 ms (1.0× b) | 416 ms |
| fast | 165 ms (2.6× n, 2.5× b) | 424 ms (1.0× b) | 415 ms |

## Per-call overhead (auto tier)

### sum (double)

`simd_sum(x)` vs base `sum(x)`; bare: `.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)`

| n | base | auto | bare |
|---|---|---|---|
| 1e0 | 0.126 µs | 0.625 µs (0.20× b) | 0.458 µs (0.28× b) |
| 1e1 | 0.209 µs | 0.667 µs (0.31× b) | 0.459 µs (0.46× b) |
| 1e2 | 1.13 µs | 0.625 µs (1.8× b) | 0.459 µs (2.5× b) |

### add (double)

`simd_add(x, y)` vs base `x + y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.084 µs | 0.709 µs (0.12× b) |
| 1e1 | 0.084 µs | 0.708 µs (0.12× b) |
| 1e2 | 0.125 µs | 0.75 µs (0.17× b) |

### exp (double)

`simd_exp(xe)` vs base `exp(xe)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.043 µs | 0.5 µs (0.09× b) |
| 1e1 | 0.042 µs | 0.501 µs (0.08× b) |
| 1e2 | 0.251 µs | 0.667 µs (0.38× b) |

### eq (double)

`simd_eq(x, y)` vs base `x == y`

| n | base | auto |
|---|---|---|
| 1e0 | 0.042 µs | 0.667 µs (0.06× b) |
| 1e1 | 0.083 µs | 0.668 µs (0.12× b) |
| 1e2 | 0.168 µs | 0.75 µs (0.22× b) |

### as_integer (double)

`simd_as_integer(xc)` vs base `as.integer(xc)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.001 µs | 0.751 µs (<0.01× b) |
| 1e1 | 0.084 µs | 0.792 µs (0.11× b) |
| 1e2 | 0.126 µs | 0.793 µs (0.16× b) |

### dot (double)

`simd_dot(x, y)` vs base `sum(x * y)`

| n | base | auto |
|---|---|---|
| 1e0 | 0.167 µs | 0.834 µs (0.20× b) |
| 1e1 | 0.292 µs | 0.834 µs (0.35× b) |
| 1e2 | 1.21 µs | 0.834 µs (1.5× b) |

### mul (complex)

`simd_mul(cx, cy)` vs base `cx * cy`

| n | base | auto |
|---|---|---|
| 1e0 | 0.084 µs | 1.96 µs (0.04× b) |
| 1e1 | 0.084 µs | 2.04 µs (0.04× b) |
| 1e2 | 0.167 µs | 2.12 µs (0.08× b) |

