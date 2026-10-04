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
* Added elementwise arithmetic: `simd_add()`, `simd_sub()`, `simd_mul()`,
  `simd_div()`, `simd_idiv()` (`%/%`), `simd_mod()` (`%%`), `simd_neg()`,
  `simd_abs()`, `simd_sign()`, `simd_copysign()`, `simd_recip()`,
  `simd_sqrt()`, `simd_fma()` (fused), `simd_mul_add()`, `simd_add_mul()`,
  `simd_lerp()`, `simd_pmin()`, `simd_pmax()`, `simd_pmin_num()`,
  `simd_pmax_num()`, `simd_clamp()`, `simd_floor()`, `simd_ceiling()`,
  `simd_trunc()` and `simd_round()`, and the wrapping integer variants
  `simd_add_wrap()`, `simd_sub_wrap()`, `simd_mul_wrap()`, `simd_neg_wrap()`
  and `simd_abs_wrap()`. Operands must have equal lengths or length one.
  Integer results, `pmin`/`pmax` and rounding are identical to base R's,
  including `NA` positions, the sign of zero and the single overflow warning;
  double `%/%` and `%%` are exact. `simd_fma()` and `simd_lerp()` round once
  on every implementation, and `simd_lerp()` is exact at `t = 0` and `t = 1`.
* Added the predicates `simd_is_na()`, `simd_is_nan()`, `simd_is_finite()`,
  `simd_is_infinite()`, `simd_is_negative()` and `simd_is_zero()`, each with
  `_any` and `_all` forms that stop at the first deciding element. They never
  return `NA`, tell `NA` from `NaN` for every payload, and follow base R for
  every type (raw is not finite, as in `is.finite()`); `simd_is_negative()`
  tests the sign bit, so `-0` is negative.
* Added the comparisons `simd_eq()`, `simd_ne()`, `simd_lt()`, `simd_le()`,
  `simd_gt()` and `simd_ge()`, identical to base R's operators for every pair
  of double, integer, logical and raw operands, `NA` positions included.
* Added three-valued logic `simd_and()`, `simd_or()`, `simd_xor()` and
  `simd_not()` (bytewise on raw vectors), and the bitwise operations
  `simd_bit_and()`, `simd_bit_or()`, `simd_bit_xor()`, `simd_bit_not()`,
  `simd_shl()`, `simd_shr()`, `simd_sar()` (arithmetic shift),
  `simd_rotl()`, `simd_rotr()` and the bit counts `simd_popcount()`,
  `simd_popcount_total()`, `simd_lzcnt()` and `simd_tzcnt()` on integer,
  logical and raw vectors. Integer results follow base R's `bitwAnd()`
  family, so a result with the bit pattern of `NA_integer_` is `NA`; raw
  shifts follow `rawShift()`.
* Added the conversions `simd_as_integer()`, `simd_as_double()`,
  `simd_as_logical()` and `simd_as_raw()`. The default `"checked"` mode
  matches `as.integer()` and `as.raw()` including their warnings;
  `"saturating"` clamps and `"truncating"` wraps around, without warnings.
* Bundled inline headers generated from SLEEF 3.9.0 (double precision only) in
  `src/vendor/sleef`, with a maintainer script that regenerates them
  reproducibly. `configure` gives each SIMD tier SLEEF's vectorised elementary
  functions when its header compiles (1.0-ULP variants, 0.5-ULP for `sinpi`,
  `cospi` and `hypot`); `RSIMD_DISABLE_SLEEF=1` turns them off. Every tier is
  compiled with `-ffp-contract=off`.
* Added the elementary functions `simd_exp()`, `simd_exp2()`, `simd_exp10()`,
  `simd_expm1()`, `simd_log()` (with a single `base`), `simd_log2()`,
  `simd_log10()`, `simd_log1p()`, `simd_pow()`, `simd_cbrt()`, `simd_hypot()`,
  `simd_sin()`, `simd_cos()`, `simd_tan()`, `simd_sincos()`, `simd_asin()`,
  `simd_acos()`, `simd_atan()`, `simd_atan2()`, `simd_sinpi()`, `simd_cospi()`,
  `simd_tanpi()`, `simd_sinh()`, `simd_cosh()`, `simd_tanh()`, `simd_asinh()`,
  `simd_acosh()` and `simd_atanh()`, using SLEEF on the SIMD tiers and the C
  math library on `"none"`. They follow base R's rules for missing values and
  its single "NaNs produced" warning; `simd_pow()` reproduces base R's `^`
  exactly at its special cases, and the pi functions are exact at integers and
  half-integers. Odd functions keep `-0` (base R 4.6 returns `0` for some).
* Added `simd_sigmoid()`, the logistic function in a form that is accurate in
  both tails (it keeps the subnormal results that `plogis()` rounds to 0), and
  `simd_softmax()` and `simd_log_softmax()` over a whole vector, computed
  stably from `exp(x - max(x))` with the denominator summed in the
  `simd_precision()` mode and no memory beyond the result. Missing and
  infinite values give what base R's formulas give, preferring `NA` over
  `NaN`, without warnings.
* Complex numbers: `simd_add()`, `simd_sub()` and `simd_neg()` take complex
  operands (a double, integer or logical operand is converted to complex),
  `simd_sum()` sums complex vectors part by part in the `simd_precision()`
  mode, and `simd_is_na()`, `simd_is_nan()`, `simd_is_finite()` and
  `simd_is_infinite()` (with their `_any` and `_all` forms) follow base R's
  rules for complex values. New `simd_conj()`, `simd_re()` and `simd_im()`.
  Missing values propagate part by part, as in base R. Other functions reject
  complex input. `simd_any_na()`, `simd_count_na()` and `simd_which_na()` now
  use SIMD kernels for complex vectors too.
* A getting-started vignette, so far covering complex numbers. Vignettes are
  built with R's own Sweave engine and need no extra packages.
* 64-bit integers: `bit64::integer64` vectors are accepted by the sums,
  extremes, scans, missing-value functions, arithmetic (checked and `_wrap`),
  `pmin()`/`pmax()`/`clamp()`, predicates, comparisons, bit operations and
  conversions, recognised by their class, so bit64 is needed only to create
  and print them. Results are exact 64-bit integers with bit64's warnings for
  overflow and division by zero; `simd_sum()` is `NA` only when the exact
  total does not fit, whatever the order of the elements. Integer and logical
  operands combine with integer64 exactly; a double operand converts the
  integer64 operands to double with a warning. New `simd_as_integer64()`
  converts doubles (with the `"checked"`, `"saturating"` and `"truncating"`
  modes), integers, logicals and raw vectors. See `?rsimd-integer64`.
* Faster `simd_idiv()`, `simd_mod()` and `simd_count_na()`. Integer `%/%`
  and `%%` by a scalar divisor no longer divide (they multiply by a
  precomputed constant, 4 to 16 lanes at a time); with a vector divisor the
  double-lane path does less work per element. `simd_count_na()` keeps
  per-lane counts instead of counting every vector's mask. On arm64 the
  neon tier uses NEON instructions directly for blends, mask tests and
  counts, the high half of 32-bit products, shifts and double-to-integer
  stores, which SIMDe's SSE2 translation took several instructions for.
  The elementwise functions also skip their R-level type checks when no
  operand has a class or is complex, which cuts the fixed cost of a call
  from about 9 to about 3 microseconds.
* `simd_vec()` wraps a vector in the optional `simd_vec` class, on which the
  arithmetic, comparison and logical operators, the `Math` and `Summary` group
  generics, `mean()` and `anyNA()` run the SIMD kernels. An object can be
  pinned to one implementation (`simd_impl(x) <- "neon"`), which every
  function given it then uses whatever `simd_use()` selected; operands pinned
  to different implementations are an error. Objects also carry a known
  NA-free flag (`simd_vec(x, check_na = TRUE)`, `simd_na_free()`) that lets
  reductions and elementwise operations skip their missing-value checks and
  is kept only by operations that cannot create a missing value. Any `simd_`
  function given a `simd_vec` returns a `simd_vec` for vector-valued results;
  reductions, predicates and comparisons return plain values. `integer64`
  data keeps class `c("simd_vec", "integer64")`. See `?simd_vec`.
* Continuous integration checks the package on Linux (x86-64 and arm64),
  macOS (arm64 and x86-64) and Windows, runs the tests once per
  implementation tier, and on a schedule under the address and undefined
  behaviour sanitizers, valgrind, an R without long double, Intel SDE
  (AVX-512) and QEMU (SVE and SVE2).
