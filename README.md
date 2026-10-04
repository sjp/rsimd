# rsimd

[![R-CMD-check](https://github.com/sjp/rsimd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjp/rsimd/actions/workflows/R-CMD-check.yaml)

`rsimd` provides SIMD-accelerated primitives for R atomic vectors: reductions,
elementwise arithmetic, predicates, bitwise operations, conversions and
transcendental functions. Kernels are written once in C on top of the portable
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
v <- simd_vec(c(1.5, NA, 3))         # operators and Summary use SIMD kernels
v * 2 + 1                            # 4 NA 7
sum(v, na.rm = TRUE)                 # 4.5
```

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

## Documentation

Four vignettes come with the package:

- `vignette("rsimd", package = "rsimd")`: getting started;
- `vignette("implementation-selection", package = "rsimd")`: CPU detection
  and choosing an implementation;
- `vignette("numerical-semantics", package = "rsimd")`: missing values,
  precision modes and every difference from base R;
- `vignette("benchmarks", package = "rsimd")`: timings against base R.

## Licence

MIT. The bundled SIMDe (MIT) and SLEEF (Boost Software License 1.0) code is
listed with its copyright holders in `inst/COPYRIGHTS`.
