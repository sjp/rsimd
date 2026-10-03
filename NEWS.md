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
* Internal access layer for R vectors: compact sequences such as `1:n` are read
  in chunks without being expanded, long vectors are supported, long
  computations can be interrupted, and binary operations accept equal lengths
  or a length-1 operand only.
* `simd_precision()` selects how floating-point sums are accumulated:
  `"fast"` (the default), `"pairwise"` or `"compensated"` (Neumaier, which
  agrees with base R's long double `sum()` without needing long double).
  The options and environment variables are documented in `?rsimd_options`.
* Internal rules for missing values and integer overflow that all functions
  follow: `NA` and `NaN` stay distinct, a reduction over both returns `NA`
  whatever their order, `na.rm = TRUE` removes both, and checked integer
  arithmetic returns `NA` with base R's overflow warning.
* Added `simd_sum()`, the sum of a double, integer or logical vector, with
  `na.rm` and `na_check` arguments and the accumulation mode set by
  `simd_precision()`. Integer sums are exact and become a double when they
  leave the integer range, as in base R.
* Added the core reductions `simd_prod()`, `simd_mean()`, `simd_min()`,
  `simd_max()`, `simd_range()`, `simd_which_min()`, `simd_which_max()`,
  `simd_any()`, `simd_all()`, `simd_any_na()`, `simd_count_na()` and
  `simd_which_na()`, with base R's result types, warnings and missing-value
  rules (`NA` wins over `NaN` whatever the order). Minimum and maximum keep
  the sign of the first zero, as in base R. `simd_mean()` refines its result
  with a second pass in the `"pairwise"` and `"compensated"` modes.
  `simd_count_na()` returns a double.
* Added `simd_sum_sq()`, `simd_sum_abs()`, `simd_norm()`, `simd_dot()`,
  `simd_dist()`, `simd_cosine()`, `simd_var()` and `simd_sd()`, one-pass
  (two for `var`/`sd`) reductions that read integer and logical elements as
  doubles without converting the vector. `dot`, `dist` and `cosine` require
  equal lengths, and their `na.rm = TRUE` drops pairs with a missing element.
  `simd_var()` and `simd_sd()` return `NA` for any missing value, as base R
  does. Sums of products use fused multiply-adds in `"fast"` mode where the
  CPU has them.
* Added the scans `simd_cumsum()`, `simd_cumprod()`, `simd_cummin()` and
  `simd_cummax()` with base R's result types, missing-value placement and
  integer overflow warning. Integer results and running minima/maxima are
  identical to base R's; `simd_cumsum()` is a sequential compensated sum in
  the `"compensated"` mode.
