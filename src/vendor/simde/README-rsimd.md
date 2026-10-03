# SIMDe subset vendored in rsimd

This directory holds a pruned copy of [SIMDe](https://github.com/simd-everywhere/simde)
(SIMD Everywhere), a header-only C library that implements x86 SIMD intrinsics
portably. rsimd's kernels are written once against SIMDe's `simde_`-prefixed x86
intrinsic names; SIMDe compiles them to native SSE/AVX instructions on x86-64 and to
native NEON instructions on arm64.

`VERSION` records the upstream repository, commit and `SIMDE_VERSION`. `COPYING` is
SIMDe's MIT licence; `inst/COPYRIGHTS` lists the per-file copyright holders and the
files that are CC0-1.0. Apart from this README and `VERSION`, every file is
byte-identical to upstream.

## What is vendored

Only the include closure of `tools/simde_probe.c`, which lists the headers rsimd's
kernels use:

- the x86 base set `x86/sse2.h`, `x86/sse4.1.h`, `x86/sse4.2.h`, `x86/avx.h`,
  `x86/avx2.h` and `x86/fma.h`, with what they include (`sse.h`, `sse3.h`, `ssse3.h`,
  `mmx.h` and SIMDe's top-level support headers);
- a curated subset of `x86/avx512/*.h` covering what kernels need for double, 32-bit
  and 64-bit integer lanes and masks.

The full SIMDe tree is about 11 MB; this subset is under 3 MB. The headers are used
only to compile the package and are not installed.

### Why there is no `arm/` tree

SIMDe's x86 headers implement their NEON code paths with the compiler's own
`<arm_neon.h>`, so `arm/neon.h` (a 5.5 MB include closure) is not needed for native
NEON code generation. SIMDe does not map x86 intrinsics onto SVE, so rsimd's SVE tiers
are written directly against the compiler's `<arm_sve.h>` and `arm/sve.h` is not
needed either. The `wasm/` and `mips/` trees and x86 extensions rsimd does not use
(SVML, XOP, GFNI, CLMUL, AES, BMI, F16C) are also left out.

## AVX-512 coverage

SIMDe implements only part of AVX-512 (about 36% of the intrinsics at the vendored
version; see the
[implementation status](https://github.com/simd-everywhere/implementation-status/blob/main/avx512.md)).
Kernels for the `avx512` tier must use intrinsics these headers provide, or guard with
`#if` and fall back to two 256-bit operations.

## Usage rules

- Always call intrinsics with the `simde_` prefix. `SIMDE_ENABLE_NATIVE_ALIASES` is
  never defined, so SIMDe names cannot collide with `<immintrin.h>` or `<arm_neon.h>`.
- `SIMDE_NO_NATIVE` and `SIMDE_ENABLE_OPENMP` are never defined.
- Each tier's translation unit includes only the headers that tier needs (see
  `tools/tiers.txt`). The scalar `none` tier does not include SIMDe.

## Diagnostic pragmas

A few upstream headers (`hedley.h`, `simde-common.h`, `x86/mmx.h`, `x86/sse.h`,
`x86/sse2.h`) contain `#pragma GCC diagnostic ignored` or `#pragma clang diagnostic
ignored` lines. They are scoped by push/pop pragmas, limited to SIMDe's own code, and
mostly silence clang-only or pedantic diagnostics (`-Wvariadic-macros`,
`-Wvector-conversion`, `-Wc11-extensions` and similar). `R CMD check` reports them in
a NOTE. They are kept so that the vendored files stay identical to upstream.

## Updating

1. Change `SIMDE_COMMIT` in `tools/vendor_simde.sh` (and, to use another header, add
   its `#include` to `tools/simde_probe.c`).
2. Run `sh tools/vendor_simde.sh`. It fetches the pinned commit, replaces this
   directory with the new closure (keeping this README), rewrites `VERSION` and
   refreshes the SIMDe section of `inst/COPYRIGHTS`. It fails if the subset exceeds
   its 3.5 MB budget.
3. Run `sh tools/check_simde_subset.sh` and `sh tools/check_size.sh`.
