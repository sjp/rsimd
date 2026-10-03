# rsimd

Status: pre-alpha. The API is incomplete and will change.

`rsimd` provides SIMD-accelerated primitives for R atomic vectors: reductions,
elementwise arithmetic, predicates, bitwise operations, conversions and
transcendental functions. Kernels are written once in C on top of the portable
[SIMDe](https://github.com/simd-everywhere/simde) intrinsics library and use
[SLEEF](https://sleef.org/) for vectorised elementary functions. They compile to
native SSE2, AVX2 and AVX-512 code on x86-64 and to NEON and SVE code on arm64.

At run time the package detects which instruction sets the CPU supports and picks
the best available implementation. You can also choose one by name, including a
plain scalar implementation that serves as the reference, through the
`rsimd.impl` option or the `RSIMD_IMPL` environment variable.

## Installation

`rsimd` is not on CRAN yet. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("sjp/rsimd")
```

A C11 compiler is required.
