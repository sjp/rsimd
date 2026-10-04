## Submission notes

### Pragmas suppressing diagnostics

`R CMD check` notes that some files under `src/vendor/simde` contain pragmas
suppressing diagnostics. These are unmodified headers of the vendored, MIT-licensed
SIMDe library (`hedley.h`, `simde-common.h`, `x86/mmx.h`, `x86/sse.h`, `x86/sse2.h`).
The pragmas are scoped by push/pop pragmas to SIMDe's own code and silence clang-only
or pedantic diagnostics (`-Wvariadic-macros`, `-Wvector-conversion`,
`-Wc11-extensions` and similar); none concerns uninitialised values, formats or array
bounds. rsimd's own sources contain no such pragmas. We keep the vendored headers
identical to upstream so that they can be updated mechanically.

### Non-portable compilation flags

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

### Bundled third-party code

`src/vendor/simde` (MIT, a subset of SIMDe) and `src/vendor/sleef` (Boost Software
License 1.0, inline headers generated from SLEEF 3.9.0) are bundled; their authors
are listed as copyright holders in `Authors@R` and their licences and copyright
notices are in `inst/COPYRIGHTS`. Both are header-only, used only at compile time
and not installed.
