# rsimd 0.0.0.9000

* Initial package skeleton with native routine registration and `simd_version()`.
* Vendored a pruned subset of the SIMDe headers (x86 base set and a curated AVX-512
  subset) in `src/vendor/simde`, with maintainer scripts to regenerate and check it.
* `simd_cpu_features()` reports the SIMD extensions the CPU and operating system support,
  detected at load time on x86, x86-64, arm64 and 32-bit arm (Linux, macOS, Windows and
  FreeBSD). The `RSIMD_CPU_FEATURES_MASK` environment variable switches features off for
  testing fallback paths.
