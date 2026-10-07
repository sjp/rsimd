# rsimd

[![R-CMD-check](https://github.com/sjp/rsimd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjp/rsimd/actions/workflows/R-CMD-check.yaml)

`rsimd` provides SIMD-accelerated primitives for R atomic vectors: reductions,
elementwise arithmetic, predicates, bitwise operations, conversions and
transcendental functions. The implementations are written once in C on top of the portable
[SIMDe](https://github.com/simd-everywhere/simde) intrinsics library and use
[SLEEF](https://sleef.org/) for vectorised elementary functions. They compile to
native SSE2, AVX2 and AVX-512 code on x86-64 and to NEON and SVE code on arm64.

At run time the package detects which instruction sets the CPU supports and picks
the best available implementation. You can also choose one by name, including a
plain scalar implementation that serves as the reference.


## Installation

```r
install.packages("rsimd")
```

or the development version from GitHub:

```sh
git clone https://github.com/sjp/rsimd.git
R CMD INSTALL rsimd
```

Installing from source needs a C11 compiler. `configure` builds every
implementation tier your compiler supports and prints which ones it built.

## Example

```r
library(rsimd)
simd_available()                     # e.g. "avx2" "sse2" "none"
x <- runif(1e6)
simd_sum(x)                          # like sum(x)
simd_exp(c(0, 1, NA))                # 1.000000 2.718282 NA
simd_add(x, 1)                       # equal lengths or a length-1 operand
simd_with_impl("none", simd_sum(x))  # the scalar reference implementation
simd_precision("compensated")        # how floating-point sums accumulate
v <- simd_vec(c(1.5, NA, 3))         # operators and Summary dispatch to rsimd
v * 2 + 1                            # 4 NA 7
sum(v, na.rm = TRUE)                 # 4.5
```

## Functions

| Family | Examples | Help |
|--------|----------|------|
| Reductions | `simd_sum()`, `simd_mean()`, `simd_min()`, `simd_range()`, `simd_which_max()`, `simd_any()`, `simd_count_na()`, `simd_max_abs()` | `?simd_sum`, `?simd_min`, `?simd_which_min` |
| Linear algebra and statistics | `simd_dot()`, `simd_norm()`, `simd_dist()`, `simd_cosine()`, `simd_var()`, `simd_sd()` | `?simd_dot`, `?simd_var` |
| Scans | `simd_cumsum()`, `simd_cumprod()`, `simd_cummin()`, `simd_cummax()` | `?simd_cumsum` |
| Elementwise arithmetic | `simd_add()`, `simd_div()`, `simd_idiv()`, `simd_mod()`, `simd_fma()`, `simd_pmin()`, `simd_clamp()`, `simd_round()`, `simd_add_wrap()` | `?simd_add`, `?simd_fma`, `?simd_round` |
| Predicates and comparisons | `simd_is_na()`, `simd_is_finite()`, `simd_is_whole()` (with `_any` and `_all` forms), `simd_eq()`, `simd_which()`, `simd_count()`, `simd_hamming()` | `?simd_is_na`, `?simd_eq`, `?simd_hamming` |
| Logic and bits | `simd_and()`, `simd_not()`, `simd_bit_and()`, `simd_shl()`, `simd_rotl()`, `simd_popcount()` | `?simd_and`, `?simd_bit_and`, `?simd_popcount` |
| Conversions | `simd_as_integer()`, `simd_as_double()`, `simd_as_raw()`, `simd_as_integer64()` | `?simd_as_integer` |
| Elementary functions | `simd_exp()`, `simd_log1p()`, `simd_pow()`, `simd_sin()`, `simd_sinpi()`, `simd_tanh()`, `simd_exp2m1()` | `?simd_exp`, `?simd_sin`, `?simd_math_accuracy` |
| Floating-point extras | `simd_ilogb()`, `simd_nextafter()`, `simd_remainder()`, `simd_rootn()`, `simd_rsqrt_approx()` | `?simd_ilogb`, `?simd_recip_approx` |
| Machine learning | `simd_sigmoid()`, `simd_softmax()`, `simd_log_softmax()` | `?simd_softmax` |
| Complex vectors | arithmetic, sums, scans and elementary functions, plus `simd_conj()`, `simd_re()`, `simd_arg()` | `?simd_conj` |
| 64-bit integers | `bit64::integer64` input to sums, extremes, arithmetic, comparisons and bit operations | `?rsimd-integer64` |
| `simd_vec` class | operators, `Math` and `Summary` call the `simd_*()` functions; `simd_unwrap()`, `simd_na_free()` | `?simd_vec` |
| Implementation selection | `simd_available()`, `simd_use()`, `simd_with_impl()`, `simd_precision()`, `simd_cpu_features()` | `?simd_use`, `?simd_precision` |

## Implementation tiers

| Tier | Instruction sets | Platforms |
|------|------------------|-----------|
| `none` | scalar C (the reference) | all |
| `sse2` | SSE2 | x86-64, i686 |
| `avx2` | AVX2 and FMA | x86-64 |
| `avx512` | AVX-512 F, BW, DQ and VL | x86-64 |
| `neon` | Advanced SIMD | arm64 (and a reduced version on 32-bit arm) |
| `sve`, `sve2` | SVE, SVE2 (any vector length) | arm64 Linux and FreeBSD |

`simd_available()` lists the tiers that are both built and supported by the
CPU, best first. `simd_use()`, the `rsimd.impl` option and the `RSIMD_IMPL`
environment variable choose one; the default `"auto"` uses the best.

## Performance

There are several benchmarks for this package some of which are documented in a vignette:

```r
vignette("benchmarks", package = "rsimd")
```

As a general rule, you should expect improvements to be more clearly observable on vectors of at least 1,000 elements.
Longer vectors (e.g. 10,000 or more) can observe significant performance improvements, often over 5 times faster.

When you have shorter vectors (smaller than 1,000 elements), base R is simpler and just as fast.

## Documentation

Four vignettes come with the package:

- `vignette("rsimd", package = "rsimd")`: getting started;
- `vignette("implementation-selection", package = "rsimd")`: CPU detection
  and choosing an implementation;
- `vignette("numerical-semantics", package = "rsimd")`: missing values,
  precision modes and every difference from base R;
- `vignette("benchmarks", package = "rsimd")`: timings against base R.

Results can differ from base R in the last bits of floating-point results, and in a few documented edge cases. `vignette("numerical-semantics")` lists every difference, and `?simd_precision` and `?simd_math_accuracy` describe the accuracy modes. Some modes, such as compensated summation, are more accurate than base R.

## Licence

MIT. The bundled SIMDe (MIT) and SLEEF (Boost Software License 1.0) code is
listed with its copyright holders in `inst/COPYRIGHTS`.
