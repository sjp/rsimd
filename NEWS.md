# rsimd 0.0.0.9000

* Initial package skeleton with native routine registration and `simd_version()`.
* Vendored a pruned subset of the SIMDe headers (x86 base set and a curated AVX-512
  subset) in `src/vendor/simde`, with maintainer scripts to regenerate and check it.
