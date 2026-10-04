# Running the rsimd tests

The tests use testthat (3rd edition). Every kernel test runs once for each
implementation tier that this machine can run (`simd_available()`). Each
tier's result is compared with the scalar `none` tier, which serves as the
oracle, and with base R.

Run the tests against the installed package. `testthat::test_local()` loads
the source tree with pkgload instead, and that breaks the tests that check
NAMESPACE and Rd consistency.

```sh
R CMD INSTALL .
cd tests/testthat
NOT_CRAN=true Rscript -e 'testthat::test_dir(".", load_package = "installed", package = "rsimd")'
```

## Test helpers

The `helper-*.R` files are loaded before the tests:

- `helper-tiers.R`: `tiers_to_test()`, `with_each_tier(fn)`,
  `for_each_tier(f)` and `skip_if_no_tier(tier)`.
- `helper-expect.R`: `expect_simd_equal(f, x, ...)` compares `f(x, ...)` on
  every tier with the `none` tier. `expect_simd_matches_base(simd_fn,
  base_fn, x, ...)` compares every tier with base R. Values that are not
  doubles must be identical. Doubles must match within a bound that depends
  on the precision mode (see the comments in the file), and their missing
  and infinite values must match exactly. `expect_simd_identical(f, x, ...)`
  requires every tier to be bit-identical to `none`, including the sign of
  zero, for results that do not depend on the order of operations (min,
  max, which, any, all, missing-value counts).
- `helper-cases.R`: `edge_lengths()`, `edge_doubles()`, `edge_ints()`,
  `edge_lgl()`, `altrep_inputs(n)`, `rand_vec(type, n, ..., seed)` and
  `with_seed(seed, code)`.
- `helper-extended.R`: `skip_unless_extended()`.
- `helper-subset.R`: `test_subset()`, `reduced_lengths()` and
  `skip_if_no_subprocess()` for the reduced runs selected by
  `RSIMD_TEST_SUBSET` (see below). Tests that start a subprocess call
  `skip_if_no_subprocess()`, and kernel tests loop over `tiers_to_test()`
  (which honours `RSIMD_TEST_TIERS`), not `simd_available()`.
- `helper-int64.R`: builds `integer64` vectors from their bit patterns
  without bit64 (`i64()`, `i64_dec()` from decimal strings, `i64_halves()`,
  `i64_from_bits()`), reads them back (`i64_str()`, `i64_parts()`,
  `i64_bits()`, `i64_value()`), and `rand_i64()`, `expect_i64()` (bit for
  bit) and `expect_tiers_i64()` (value and warnings on every tier). The
  integer64 tests that compare with bit64 itself are in
  `test-integer64-bit64.R` and are skipped when bit64 is not installed.

## Environment variables

| Variable | Effect |
|----------|--------|
| `NOT_CRAN=true` | Runs the tests marked `skip_on_cran()`, such as subprocess tests and the C lint. |
| `RSIMD_EXTENDED_TESTS=true` | Also runs the long and randomised tests. They take a few seconds more. |
| `RSIMD_TEST_SUBSET=quick` | Caps input lengths at 1000 for slow environments (valgrind, emulators); the cases around every vector width are kept. |
| `RSIMD_TEST_SUBSET=tier_emulation` | As `quick`, and skips tests that start a subprocess (an emulator does not follow it). |
| `RSIMD_TEST_TIERS=avx512,none` | The kernel tests loop over these tiers only (`none` is always added as the oracle); a tier that is not available is an error. |
| `RSIMD_IMPL`, `RSIMD_CPU_FEATURES_MASK`, `RSIMD_DISABLE_TIERS`, `RSIMD_DEBUG_STRIDE` | Change the implementation in use, the CPU features detected, the tiers built and the chunk size (see `?rsimd_options`). |

The C lint test runs `tools/lint_c.sh`, so it needs a source checkout. It is
skipped under `R CMD check`.

## Sanitizers (ASan and UBSan)

Build the package with sanitizer flags into a scratch library. Put the
flags in a throwaway Makevars file selected with `R_MAKEVARS_USER`, so that
your own `~/.R/Makevars` is left alone. The package's tier flags are still
added by `configure`. Never pass `-march` flags yourself.

```sh
cat > /tmp/Makevars-san <<'EOF'
CC = gcc -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer
EOF
mkdir -p /tmp/sanlib
R_MAKEVARS_USER=/tmp/Makevars-san R CMD INSTALL --no-test-load -l /tmp/sanlib .
```

`--no-test-load` is needed because R itself is not built with the
sanitizers, so the runtime libraries have to be preloaded:

```sh
cd tests/testthat
LD_PRELOAD="$(gcc -print-file-name=libasan.so) $(gcc -print-file-name=libubsan.so)" \
  ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=print_stacktrace=1 \
  NOT_CRAN=true RSIMD_EXTENDED_TESTS=true R_LIBS=/tmp/sanlib \
  Rscript -e 'testthat::test_dir(".", load_package = "installed", package = "rsimd")'
```

Any `runtime error:` or `AddressSanitizer` line in the output is a failure.
`detect_leaks=0` turns off reports of R's own allocations, which are not
leaks.

## valgrind

```sh
R CMD INSTALL .
cd tests/testthat
NOT_CRAN=true R -d "valgrind --track-origins=yes" \
  -e 'testthat::test_dir(".", load_package = "installed", package = "rsimd")'
```

Check the `ERROR SUMMARY` at the end. Subprocess tests (callr) run their
child R outside valgrind.

## Continuous integration

`tools/ci/run_tests.sh` and `tools/ci/test_tiers.sh` wrap the commands above
for the CI workflows, including runs under Intel SDE and QEMU; see
`CONTRIBUTING.md`.

## Reference environment

CRAN's sanitizer checks run on an R that is itself built with ASan and
UBSan. The rocker `r-devel-san` image provides one
(<https://rocker-project.org/images/base/r-devel.html>; see also "Checking
memory access" in Writing R Extensions). With that R, there is no need to
preload the libraries.

## Benchmarks

`bench/sum.R` (not part of the built package) times `simd_sum()` on every
available tier against base `sum()` for 1e7 elements, and
`bench/reductions.R` does the same for `simd_min()` and `simd_any_na()`
(and the integer `simd_min()` with `na_check = FALSE`), and
`bench/linalg.R` for `simd_dot()` and `simd_cumsum()`. They are quick
sanity checks, not a benchmark suite.
