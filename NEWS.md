# rsimd 0.0.0.9000

* Initial package skeleton with native routine registration and `simd_version()`.
* Vendored a pruned subset of the SIMDe headers (x86 base set and a curated AVX-512
  subset) in `src/vendor/simde`, with maintainer scripts to regenerate and check it.
* `simd_cpu_features()` reports the SIMD extensions the CPU and operating system support,
  detected at load time on x86, x86-64, arm64 and 32-bit arm (Linux, macOS, Windows and
  FreeBSD). The `RSIMD_CPU_FEATURES_MASK` environment variable switches features off for
  testing fallback paths.
* A `configure` script (also used on Windows) probes which instruction-set tiers
  the compiler can build (`sse2`, `avx2`, `avx512` on x86; `neon`, `sve`, `sve2`
  on arm) and compiles each tier's kernels in a separate file with that tier's
  flags. `RSIMD_DISABLE_TIERS` skips tiers at install time.
* Implementation selection: `simd_tiers()` lists the tier ids,
  `simd_available()` the tiers that are compiled in and supported by the CPU
  (best first), `simd_current()` the one in use, and `simd_use()` selects one
  by name or `"auto"`. `simd_with_impl()` evaluates an expression with a tier
  selected temporarily. The selection can also be set with the `rsimd.impl`
  option or, at load time, the `RSIMD_IMPL` environment variable. Operations a
  tier lacks use the next lower tier's version.
