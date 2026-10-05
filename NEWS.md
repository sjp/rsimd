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
  `simd_is_negative()` and `simd_is_zero()` with `_any` and `_all` forms;
  `simd_eq()` to `simd_ge()`; `simd_and()`, `simd_or()`, `simd_xor()`,
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
  values and warnings. Machine-learning helpers `simd_sigmoid()`,
  `simd_softmax()` and `simd_log_softmax()`.
* **Complex vectors:** `simd_add()`, `simd_sub()`, `simd_neg()`,
  `simd_sum()` and the predicates take complex input; `simd_conj()`,
  `simd_re()` and `simd_im()`.
* **64-bit integers:** `bit64::integer64` vectors are supported by the sums,
  extremes, scans, arithmetic, predicates, comparisons, bit operations and
  conversions, with exact 64-bit results and bit64's warnings;
  `simd_as_integer64()` converts to them (`?rsimd-integer64`).
* **`simd_vec` class:** `simd_vec()` wraps a vector so that the operators,
  the `Math` and `Summary` group generics, `mean()` and `anyNA()` use the
  SIMD kernels. An object can be pinned to one implementation
  (`simd_impl<-`) and can carry a known NA-free flag that lets functions
  skip their missing-value checks (`?simd_vec`).
* Results follow base R's types, warnings and missing-value rules; `NA` and
  `NaN` stay distinct. Where a result can differ from base R's (the last
  bits of floating-point sums and elementary functions, and a few
  documented edge cases) the vignette "Numerical semantics and precision"
  lists it. The vignettes "Getting started with rsimd", "Choosing an
  implementation" and "Benchmarks" cover the rest.
