# rsimd 0.1.0

First release.

* **Implementation selection.** CPU features are detected at load time
  (`simd_cpu_features()`) on x86-64, arm64 and 32-bit x86 and arm, under
  Linux, macOS, Windows and FreeBSD. `configure` builds every instruction-set
  tier the compiler supports (`sse2`, `avx2` and `avx512` on x86; `neon`,
  `sve` and `sve2` on arm), each in its own file with its own flags, and
  always the scalar `none` tier. `simd_available()` lists the tiers this
  machine can run, `simd_use()` selects one by name or `"auto"`,
  `simd_with_impl()` selects one temporarily, and the `rsimd.impl` option and
  `RSIMD_IMPL` environment variable set it too (`?rsimd_options`).
* **Reductions:** `simd_sum()`, `simd_prod()`, `simd_mean()`, `simd_min()`,
  `simd_max()`, `simd_range()`, `simd_which_min()`, `simd_which_max()`,
  `simd_any()`, `simd_all()`, `simd_any_na()`, `simd_count_na()`,
  `simd_which_na()`, `simd_sum_sq()`, `simd_sum_abs()`, `simd_norm()`,
  `simd_dot()`, `simd_dist()`, `simd_cosine()`, `simd_var()`, `simd_sd()`,
  and the scans `simd_cumsum()`, `simd_cumprod()`, `simd_cummin()` and
  `simd_cummax()`. `simd_precision()` chooses how floating-point sums are
  accumulated: `"fast"`, `"pairwise"` or `"compensated"` (Neumaier).
  `simd_which()` and `simd_count()` are `which()` and `sum()` of a logical
  vector; `simd_max_abs()`, `simd_min_abs()`, `simd_which_max_abs()` and
  `simd_which_min_abs()` are `max(abs(x))` and its kin without the vector
  `abs(x)`, also for complex and integer64 input; `simd_prod_sums()` and
  `simd_prod_diffs()` are `prod(x + y)` and `prod(x - y)` without the
  vector of sums.
* **Elementwise arithmetic:** `simd_add()`, `simd_sub()`, `simd_mul()`,
  `simd_div()`, `simd_idiv()`, `simd_mod()`, `simd_neg()`, `simd_abs()`,
  `simd_sign()`, `simd_copysign()`, `simd_recip()`, `simd_sqrt()`,
  `simd_fma()`, `simd_mul_add()`, `simd_add_mul()`, `simd_lerp()`,
  `simd_pmin()`, `simd_pmax()`, `simd_pmin_num()`, `simd_pmax_num()`,
  `simd_clamp()`, `simd_floor()`, `simd_ceiling()`, `simd_trunc()`,
  `simd_round()`, and wrapping integer variants (`simd_add_wrap()` and
  friends). Operands must have equal lengths or length one.
* **Predicates, comparisons, logic, bits and conversions:** `simd_is_na()`,
  `simd_is_nan()`, `simd_is_finite()`, `simd_is_infinite()`,
  `simd_is_negative()`, `simd_is_zero()` and the number classes
  `simd_is_normal()`, `simd_is_subnormal()`, `simd_is_whole()`,
  `simd_is_even()`, `simd_is_odd()` and `simd_is_pow2()`, with `_any` and
  `_all` forms; `simd_eq()` to `simd_ge()`; the Hamming distances
  `simd_hamming()` (`sum(x != y)` without a logical vector) and
  `simd_hamming_bits()`; `simd_and()`, `simd_or()`, `simd_xor()`,
  `simd_not()`; `simd_bit_and()` and the other bitwise operations, shifts,
  rotations and bit counts; `simd_as_integer()`, `simd_as_double()`,
  `simd_as_logical()` and `simd_as_raw()` with `"checked"`, `"saturating"`
  and `"truncating"` modes.
* **Elementary functions:** exponentials, logarithms, powers, roots,
  trigonometric, pi-scaled and hyperbolic functions (`simd_exp()` to
  `simd_atanh()`), vectorised with the bundled SLEEF 3.9.0 headers on the
  SIMD tiers and the C math library on `"none"`. `simd_math_accuracy("fast")`
  (option `rsimd.math_accuracy`) switches the SIMD tiers to SLEEF's faster
  3.5-ULP variants where they exist, with the same missing values, special
  values and warnings. On arm64 (`"neon"`, `"sve"`, `"sve2"`), `simd_log()`,
  `simd_log2()`, `simd_cosh()`, `simd_asinh()`, `simd_acosh()` and
  `simd_pow()` call the C math library in accurate mode, which is faster
  there than SLEEF's 128-bit vector code, so they are identical to base R
  (`?simd_math_accuracy`).
  `simd_exp2m1()`, `simd_exp10m1()`, `simd_log2p1()` and `simd_log10p1()`
  (C23's `exp2m1()` and the others) compute `2^x - 1`, `10^x - 1`,
  `log2(1 + x)` and `log10(1 + x)` without losing accuracy near 0, within
  2 ULP, and `simd_sincospi()` gives `simd_sinpi()` and `simd_cospi()` in one
  pass. Machine-learning helpers `simd_sigmoid()`,
  `simd_softmax()` and `simd_log_softmax()`.
* **Floating-point extras:** `simd_ilogb()`, `simd_scaleb()`,
  `simd_nextafter()`, `simd_next_up()`, `simd_next_down()`,
  `simd_remainder()` (IEEE remainder) and `simd_rsqrt()`, identical to C and
  base R on every tier; `simd_rootn()`, the real n-th root, within 1 ULP and
  exact for exact roots. Approximations `simd_recip_approx()` and
  `simd_rsqrt_approx()` (within a relative error of 2^-22, from hardware
  estimates and Newton steps) and `simd_mul_add_approx()` (fused where the
  CPU has a fused multiply-add).
* **Complex vectors:** `simd_add()`, `simd_sub()`, `simd_mul()`,
  `simd_div()`, `simd_neg()`, `simd_abs()` (the modulus), `simd_sum()`,
  `simd_prod()`, `simd_mean()`, `simd_cumsum()`, `simd_cumprod()`,
  `simd_eq()`, `simd_ne()` and the predicates take complex input;
  `simd_conj()`, `simd_re()`, `simd_im()` and `simd_arg()`. Multiplication,
  division, `cumsum()` and `cumprod()` are identical to base R's results,
  special values included: base R's rounding depends on the compiler that
  built R, so the variant is chosen when the package loads by comparing
  with base R (`attr(simd_current(), "complex")`). The elementary functions
  `simd_sqrt()`, `simd_exp()`, `simd_log()` (also with a base),
  `simd_log2()`, `simd_log10()`, `simd_pow()`, `simd_sin()` to
  `simd_atan()`, `simd_sinh()` to `simd_atanh()` and `simd_atan2()` take
  complex input too: the `none` implementation is base R exactly, the others
  are within 2 to 5 ULP (stated per function in `?simd_exp`) with base R's
  special values and branch cuts (vectorised: base R's own branch-cut
  formulas, within 2 ULP of base R), and whole-number powers are identical
  to base R everywhere. `simd_acosh()` returns the principal value (C's
  `cacosh`) where base R's `acosh()` has a negative real part. `simd_vec`
  objects use the kernels for complex `*`, `/`, `^` and the Math group.
* **64-bit integers:** `bit64::integer64` vectors are supported by the sums,
  extremes, scans, arithmetic, predicates, comparisons, bit operations and
  conversions, with exact 64-bit results and bit64's warnings;
  `simd_as_integer64()` converts to them (`?rsimd-integer64`). Making an
  integer64 result loads bit64, so that its methods sort, print and compare
  the result correctly, and is an error, with an install hint, if bit64 is
  not installed.
* **`simd_vec` class:** `simd_vec()` wraps a vector so that the operators,
  the `Math` and `Summary` group generics, `mean()` and `anyNA()` use the
  SIMD kernels. An object can be pinned to one implementation
  (`simd_impl<-`) and can carry a known NA-free flag that lets functions
  skip their missing-value checks. The flag is tied to the object it was set
  on, so base functions that copy attributes onto new data (`pmin()`,
  `storage.mode<-`, `unclass()`, `readRDS()`, ...) can never carry a stale
  one; `Re()`, `Im()`, `Mod()`, `Arg()`, `Conj()`, `is.na()`, `[[<-`,
  `all.equal()` and `range(finite = TRUE)` have methods (`?simd_vec`).
  `simd_unwrap()` and `bit64::as.integer64()` return the plain data, an
  integer64 one keeping its class. A `simd_vec` has no names or
  dimensions: setting them returns the plain data with them set, so
  `quantile()`, `summary()` and `setNames()` give base R's results.
* Classed data other than integer64 and `simd_vec` (dates, times,
  `difftime`, factors, ...) is an error in every function and operator,
  rather than being treated as its bare numbers.
* Edge cases made consistent across implementations before release:
  `simd_idiv()` on doubles is the exact floor of the quotient, rounded once,
  for quotients of 2^52 and more too (it could be a few units off there);
  `simd_var()` and `simd_sd()` give `Inf` on every implementation when
  finite values overflow the mean (some lengths gave `NaN`);
  `simd_next_up()` and `simd_next_down()` return every NaN input as it is
  (some NaNs became infinities); and `simd_sin()`, `simd_tan()`,
  `simd_asin()`, `simd_atan()`, `simd_sinh()`, `simd_tanh()`,
  `simd_asinh()`, `simd_atanh()`, `simd_expm1()` and `simd_log1p()` return
  the smallest subnormal numbers themselves, sign included, in both
  accuracy modes (some gave zero); the complex `simd_tan()` and
  `simd_tanh()` of numbers with such a part keep its sign in the result
  (SVE gave a zero of the other sign). A function called from a warning handler
  in the middle of another rsimd call no longer changes that call's
  implementation or `simd_vec` result.
* Results follow base R's types, warnings and missing-value rules; `NA` and
  `NaN` stay distinct. Where a result can differ from base R's (the last
  bits of floating-point sums and elementary functions, and a few
  documented edge cases) the vignette "Numerical semantics and precision"
  lists it. The vignettes "Getting started with rsimd", "Choosing an
  implementation" and "Benchmarks" cover the rest.
