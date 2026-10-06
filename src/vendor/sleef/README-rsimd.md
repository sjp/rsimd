# SLEEF inline headers vendored in rsimd

This directory holds inline headers generated from [SLEEF](https://sleef.org/) (SIMD
Library for Evaluating Elementary Functions), which rsimd uses for vectorised
elementary functions (`exp`, `log`, `pow`, trigonometric and hyperbolic functions).
SLEEF is distributed under the Boost Software License 1.0 (`LICENSE.txt`);
`inst/COPYRIGHTS` lists the copyright holders.

SLEEF does not ship these headers: its CMake build generates them
(`-DSLEEF_BUILD_INLINE_HEADERS=TRUE`). `tools/vendor_sleef.sh` runs that build at the
tag pinned in the script and writes the files here; `VERSION` records the tag, commit,
generation date and compilers. The headers are used only to compile the package and
are not installed.

## What is vendored

One header per SIMD tier, each self-contained (`static inline` functions, with the
Payne-Hanek reduction tables as `static const` arrays; nothing is linked):

| Header | rsimd tier | Function names |
|--------|-----------|----------------|
| `sleefinline_sse2.h` | `sse2` (x86-64 and i686) | `Sleef_<f>d2_<acc>sse2` |
| `sleefinline_avx2.h` | `avx2` | `Sleef_<f>d4_<acc>avx2` |
| `sleefinline_avx512f.h` | `avx512` | `Sleef_<f>d8_<acc>avx512f` |
| `sleefinline_advsimd.h` | `neon` (arm64 only) | `Sleef_<f>d2_<acc>advsimd` |
| `sleefinline_sve.h` | `sve`, `sve2` | `Sleef_<f>dx_<acc>sve` |

The `none` tier uses C99 libm and 32-bit ARM has no double-precision SIMD, so neither
has a header. `tools/tiers.txt` maps tiers to headers; configure test-compiles each
header with its tier's flags and builds the tier without SLEEF if that fails. Kernels
do not call these functions directly: `src/kernels/common.inc.h` wraps them as
`rsimd_sleef_<f>()` on `rsimd_vf64` and states the accuracy policy.

## Changes from SLEEF's output

`tools/vendor_sleef.sh` edits the generated headers so that they hold double precision
only, do not depend on the machine that generated them, compile cleanly with
`-Wall -pedantic` and have no undefined behaviour that sanitizers report:

- the single-precision (`float`) functions are not generated (about a third of the
  size);
- `SLEEF_FLOAT128_IS_IEEEQP` and `SLEEF_LONGDOUBLE_IS_IEEEQP`, set from a C++ compile
  test on the generating machine and used only by the quad-precision type, are left
  undefined;
- `#pragma STDC FP_CONTRACT OFF` (ignored by GCC with a warning) is removed; the tier
  files are compiled with `-ffp-contract=off`, which SLEEF requires;
- helpers declared plain `static` are made `static inline`;
- the advsimd and sve headers get back the ARM Ltd. copyright notice of the SLEEF
  source files they are generated from;
- in the advsimd header, the 64-bit mask add, subtract and negate helpers use
  unsigned lanes (`vaddq_u64`, `vsubq_u64`) instead of signed ones. `nextafter`
  relies on their wraparound, which is undefined for the signed intrinsics because
  GCC implements them as C arithmetic (UBSan reports it); the instructions are the
  same.

Both the native and the cross-compiled generation routes give byte-identical headers.

## Known limits of SLEEF's functions

Measured with `tools/check_vector_layer.sh` against long double libm, the wrappers are
within SLEEF's documented bounds (1.0 ULP for `_u10`, 0.5 ULP for `_u05`, and
`max(bound, DBL_MIN)` for `atan2`, `sinpi` and `cospi`) on every tier. Outside these
ranges the kernels must handle values themselves:

- `sinh` and `cosh` keep their bound on [-709, 709] only, and `asinh` and `acosh` on
  [-1.34e154, 1.34e154]; beyond that `asinh` and `acosh` return infinities;
- `sinpi`, `cospi` and `sincospi` return 0 and 1 for |x| > 2.5e8 (SLEEF's
  documentation says 1e9), and `sinpi(n)` is -0 for odd n > 0;
- NaN payloads, so R's `NA_real_`, are not preserved.

## Updating

1. Change `SLEEF_TAG` and `SLEEF_COMMIT` in `tools/vendor_sleef.sh`.
2. Run `sh tools/vendor_sleef.sh` on Linux (x86-64 or arm64; the other architecture
   is cross-compiled, see the script) or run the `vendor-sleef` GitHub Actions
   workflow and take its artifact. The script rewrites the headers, `VERSION` and
   `LICENSE.txt` (keeping this README) and refreshes the SLEEF section of
   `inst/COPYRIGHTS`. It fails if the headers exceed their 2.5 MB budget.
3. Rebuild the package and run `sh tools/check_vector_layer.sh` and
   `sh tools/check_size.sh`.
