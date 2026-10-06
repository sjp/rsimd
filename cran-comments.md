## Submission

This is a new submission.

rsimd provides SIMD-accelerated primitives for R atomic vectors. It compiles
its kernels once per instruction-set tier (SSE2, AVX2, AVX-512 on x86-64;
NEON, SVE, SVE2 on arm64; plus a scalar tier) and chooses a tier at run time
after detecting what the CPU and operating system support. The sections
below explain the NOTEs this arrangement causes and answer the questions we
expect it to raise.

## Test environments

* Local: Debian 13 (trixie) aarch64, R 4.6.1, GCC 16.2 and GCC 14,
  Debian clang 21.1 (tiers none, neon, sve, sve2 built; none and neon run
  natively, sve and sve2 under QEMU).
* GitHub Actions, R CMD check --as-cran: Ubuntu x86-64 (R-release, R-devel,
  R-oldrel), Ubuntu arm64 (R-release), macOS arm64 and x86-64 (R-release),
  Windows x86-64 (R-release).
  TODO: confirm all seven jobs pass on the release commit.
* GitHub Actions, scheduled: ASan/UBSan (GCC and clang), valgrind, an R built
  without long double (noLD), Intel SDE (AVX-512 tier) and QEMU (SVE and SVE2
  tiers).
  TODO: confirm the latest scheduled runs pass on the release commit.
* R-hub: gcc-asan, clang-asan, valgrind, nold, rchk, macos-arm64, windows.
  TODO: fill in from the R-hub results.
* win-builder: R-devel, R-release, R-oldrel.
  TODO: fill in from the win-builder results.

## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new submission.
* Files which contain pragma(s) suppressing diagnostics (in the vendored
  SIMDe headers): see "Pragmas suppressing diagnostics" below.
* Compilation used the following non-portable flag(s) (the instruction-set
  flags of the tier files): see "Non-portable compilation flags" below.

The installed size is reported as INFO on Linux where the shared library
keeps its debugging information (30.5 MB on aarch64, of which libs 29.6 MB).
Without debugging information the library is 2.0 MB on aarch64 and the
installed package 3.0 MB. On x86-64, where four tiers (none, SSE2, AVX2,
AVX-512) are built, the library without debugging information is 3.0 MB and
the installed package 3.9 MB (measured with a GCC cross-compiler).
TODO: confirm the x86-64 figures with the size-check workflow.

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Anticipated questions

| Question | Answer |
|----------|--------|
| Why do the compiler flags include `-mavx2`, `-mavx512f` or `-march=armv8-a+sve`? Writing R Extensions advises against `-march`. | Only the files `src/tier_<tier>.c`, which contain the SIMD kernels, get instruction-set flags, through explicit rules that `configure` writes into `src/Makevars`. Their code is reached only through a dispatch table, and only after a run-time check of the CPU and operating system (cpuid and xgetbv on x86, getauxval/elf_aux_info on Linux and FreeBSD, sysctl on macOS, IsProcessorFeaturePresent on Windows). Every other file is compiled with R's flags alone and never calls tier code directly. `-march=native` is never used. The CRAN package RsimdDispatch uses the same arrangement. |
| What if the compiler cannot build a tier? | `configure` test-compiles each tier with R's compiler and flags and leaves out any that fail. The install log lists the tiers built and, for each tier not built, the reason (for example `rsimd: tiers not built: avx512 (compiler rejects -mavx512f -mavx512bw -mavx512dq -mavx512vl -mfma)`). The scalar tier and the architecture's baseline tier (SSE2 on x86-64, NEON on arm64) are always built. |
| What if the CPU does not support a tier that was built? | Run-time detection leaves it out of `simd_available()`, and the default `"auto"` picks the best tier the CPU supports. The tests also run under Intel SDE emulating a CPU with only SSE-era extensions (Merom) and under QEMU emulating an arm64 CPU without SVE (Cortex-A72), where only the lower tiers are available. |
| Bundled code and size? | `src/vendor/simde` is a pruned subset of SIMDe (MIT, 3.0 MB of headers) and `src/vendor/sleef` the inline headers generated from SLEEF 3.9.0 (Boost Software License 1.0, 1.9 MB). Both are header-only, used only at compile time and not installed. Their authors are copyright holders in `Authors@R`, and their licences and copyright notices are in `inst/COPYRIGHTS`. The source tarball is 1.3 MB. Installed size: see above. |
| Do results differ from base R? | In the last bits of floating-point sums (the accumulation order depends on the vector width; `simd_precision()` offers pairwise and compensated modes) and of elementary functions (within SLEEF's 1-ULP bound, or 3.5 ULP if the user opts into `simd_math_accuracy("fast")`), and in a few documented edge cases. The vignette "Numerical semantics and precision" lists every difference. The tests compare every tier with the scalar tier and with base R, with tolerances where results may differ. |
| Do the tests depend on the hardware? | The tests loop over the tiers available on the machine and always include the scalar tier, so they pass whatever the CPU supports. On CRAN they test the best tier and the scalar tier only, which keeps the run time independent of the number of tiers (about 80 s on the aarch64 test machine). |
| Sanitizers, valgrind, noLD? | Checked in continuous integration for every tier (see Test environments) and on R-hub. |
| Threads? | None are used. |
| Alternative BLAS? | Not relevant: the package does not use BLAS or LAPACK. |

## Pragmas suppressing diagnostics

`R CMD check` notes that some files under `src/vendor/simde` contain pragmas
suppressing diagnostics. These are unmodified headers of the vendored, MIT-licensed
SIMDe library (`x86/mmx.h`, `x86/sse2.h`). The pragmas are scoped by push/pop
pragmas to SIMDe's own code and silence clang-only or pedantic diagnostics
(`-Wvariadic-macros`, `-Wvector-conversion`, `-Wc11-extensions` and similar);
none concerns uninitialised values, formats or array bounds. rsimd's own sources
contain no such pragmas. We keep the vendored headers identical to upstream so
that they can be updated mechanically.

## Non-portable compilation flags

`R CMD check` notes instruction-set flags such as `-mavx2 -mfma`,
`-mavx512f -mavx512bw -mavx512dq -mavx512vl` and `-march=armv8-a+sve`. Each of them
is applied to a single translation unit, `src/tier_<tier>.c`, through an explicit
rule that `configure` writes into `src/Makevars` only if a test compile with R's
compiler and flags succeeds. Code in those files runs only after a run-time check
that the CPU and operating system support the instruction set
(`src/cpu_features.c`); every other file is compiled with R's flags alone, so the
package runs on any CPU of the target architecture, and `-march=native` is never
used. Per-file flags are needed because the vendored SIMDe headers select native
instructions at include time from the compiler's predefined macros. The CRAN
package RsimdDispatch uses the same arrangement. Any other flag in that note
(for example `-mbranch-protection=standard` on Debian) comes from R's own
configuration.

## Bundled third-party code

`src/vendor/simde` (MIT, a subset of SIMDe) and `src/vendor/sleef` (Boost Software
License 1.0, inline headers generated from SLEEF 3.9.0) are bundled; their authors
are listed as copyright holders in `Authors@R` and their licences and copyright
notices are in `inst/COPYRIGHTS`. Both are header-only, used only at compile time
and not installed.

## Resubmission history

None yet.
