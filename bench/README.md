# Benchmarks

This directory holds the rsimd benchmark suite and a few smaller smoke scripts. Nothing
here is part of the built package (`bench/` is in `.Rbuildignore`), nothing runs during
`R CMD check`, and no test asserts a timing. The numbers are for people: to see what
each tier gains and to spot regressions and misconfiguration by eye.

## Running the suite

Install the package and `bench`, then from the package root:

```sh
R CMD INSTALL .
Rscript bench/run.R --quick           # sizes 1e3 and 1e5, a minute or two at most
Rscript bench/run.R                   # sizes 1e3, 1e5 and 1e7, a few minutes
Rscript bench/run.R --ops sum,exp --sizes 1e4,1e6 --out /tmp/bench
```

| Option | Meaning |
|--------|---------|
| `--quick` | Sizes 1e3 and 1e5 and half the timing budget (at least 10, at most 100 iterations). |
| `--ops a,b` | Only these operations: `sum`, `mean`, `dot`, `add`, `fma`, `pmax`, `exp`, `any_na`, `is_na`, `as_integer`. |
| `--sizes a,b` | Input lengths (default `1e3,1e5,1e7`: L1-resident, cache-resident, DRAM-bound). |
| `--out dir` | Output directory (default `bench/results/`, which git ignores). |

Every operation runs once per tier in `simd_available()` (through `simd_with_impl()`)
and once as the base R equivalent:

| Op | rsimd | Base R | Types |
|----|-------|--------|-------|
| sum | `simd_sum(x)` | `sum(x)` | double, integer |
| mean | `simd_mean(x)` | `mean(x)` | double |
| dot | `simd_dot(x, y)` | `sum(x * y)`, plus `crossprod(x, y)[1]` as a BLAS reference | double |
| add | `simd_add(x, y)` | `x + y` | double, integer |
| fma | `simd_fma(x, y, z)` | `x * y + z` | double |
| pmax | `simd_pmax(x, y)` | `pmax(x, y)` | double |
| exp | `simd_exp(x)` | `exp(x)` | double |
| any_na | `simd_any_na(x)` | `anyNA(x)` | double, integer |
| is_na | `simd_is_na(x)` | `is.na(x)` | double |
| as_integer | `simd_as_integer(x, mode = "truncating")` | `as.integer(x)` | double |

There are four tables:

- **Main**: every op at every size, precision mode `"fast"`, no missing values.
- **1% NA**: `sum` and `any_na` with 1% of the elements set to `NA`, which exercises the
  NA-aware paths. Both base R and rsimd stop `anyNA()` at the first `NA` (base R's
  integer `sum()` does too), so those rows measure call overhead, not throughput.
- **Precision modes**: `sum` and `dot` under `"fast"`, `"pairwise"` and `"compensated"`
  at the largest size, showing what the more accurate modes cost.
- **Math accuracy modes**: `sin`, `log`, `tanh`, `atan2` and `hypot` under
  `simd_math_accuracy()` `"accurate"` and `"fast"` at the largest size, showing what
  SLEEF's 3.5-ULP variants gain. The `none` tier and base R have no modes, so their rows
  differ only by noise. Base R's `hypot` is `sqrt(x * x + y * y)`, which can overflow.

Inputs are generated once per size with `set.seed(20261003)`: doubles from
`runif(n, -100, 100)` (`x`, `y`, `z` independent), `exp` on `runif(n, -50, 50)`,
`as_integer` on `runif(n, -1e6, 1e6)`, integers from `sample.int(2e6, n, TRUE) - 1e6`,
`log` on `10^runif(n, -3, 3)` and `tanh` on `runif(n, -5, 5)`.
Each timing is `bench::mark(check = FALSE, filter_gc = TRUE, memory = FALSE)` with at least
20 and at most 200 iterations and a minimum time of 0.1 s, 0.2 s or 0.5 s for
n < 1e5, < 1e7 and ≥ 1e7. `check = FALSE` because results legitimately differ in the last
bits between precision modes. Correctness is the test suite's job.

## Output

Each run writes three files named `<UTC timestamp>-<os>-<arch>-<git sha>`, and copies the
markdown to `latest.md`:

- `.rds`: a list with `results` (the flat table), `marks` (every `bench_mark` object,
  named `table/op/type/mode/n/impl`) and `metadata`: rsimd, R and bench versions, git SHA
  (with `-dirty` for uncommitted changes), `R CMD config CC` and `CFLAGS`, OS, CPU model,
  core count, `simd_cpu_features()`, `simd_available()`, the tier `"auto"` selects, the
  machine key used by `compare.R`, the sizes and the run time.
- `.csv`: the flat table: `table, op, type, mode, n, impl, median_us, itr_sec, n_itr,
  speedup_none, speedup_base, flag, rsimd_version, timestamp, machine`.
- `.md`: a metadata header and one table per op, rows for sizes (or modes),
  columns for tiers then base R.

## Reading `latest.md`

A cell such as `412 µs (3.1× n, 5.6× b)` is the median time per call, then the speedup
against the `none` tier (`n`) and against base R (`b`). Above 1× means faster. `none`
cells show only the speedup against base R, and base R cells show only the time.

A SIMD tier that is less than 1.1× faster than `none` on a compute-bound op (`sum`, `mean`,
`dot`, `exp`, `pmax`, `any_na`) at n ≥ 1e5 in the main table is marked ⚠, and `run.R`
prints a warning. That pattern usually means the tier is running scalar code: a kernel
missing from the tier, a dispatch fallback, or a build that did not compile the tier with
its flags. It is a hint to investigate, never a failure.

Things that are expected and not a problem:

- **At n = 1e3 the speedups shrink or invert.** Every rsimd call pays a fixed few
  microseconds of R-side work (argument and type checks, the tier synchronisation and
  `.Call`). Base R primitives such as `x + y` and `anyNA(x)` pay well under a microsecond.
- **Memory-bound ops are near parity with base R at large n.** `add`, `fma`, `is_na`
  and `as_integer` at 1e7 are limited by memory bandwidth and by allocating the result,
  which both sides do. The win over `none` stays visible at 1e5.
- **`sum` against base R depends on the platform's `long double`.** Base R accumulates
  double sums in `long double`. On x86-64 that is the 80-bit x87 format, which is slow
  but in hardware. On arm64 Linux it is 128-bit quad precision done in software, so base
  `sum()` and `mean()` are tens of times slower than on x86 (84 ms against 1.4 ms for
  `simd_sum` on 1e7 doubles on an Apple M-series core running Linux). On macOS arm64,
  `long double` is plain double. Ratios against base R for these ops are therefore not
  comparable across platforms.
- **`exp` on arm64 is close to base R.** glibc's aarch64 `exp` is already fast, and SLEEF's
  `u10` `exp` does not beat it by much. On x86-64, SLEEF against a scalar libm loop shows
  the largest gains in the suite.

### Indicative targets

Rough expectations on a modern x86-64 with AVX2 at n = 1e7 (double unless stated). They
are for spotting gross misconfiguration, not pass/fail thresholds.

| Op | avx2 vs base R | avx2 vs `none` | Comment |
|----|----------------|----------------|---------|
| sum (double) | 4–10× | 3–6× | base uses x87 long double; `fast` mode uses 4+ accumulators to break the dependency chain |
| sum (integer) | 2–4× | 2–4× | base accumulates in 64-bit scalar; the NA sentinel check costs a compare |
| mean | 3–8× | 2–4× | two passes in base |
| dot | 3–6× | 3–5× | vs `sum(x*y)`, which allocates a temporary; vs BLAS `crossprod` roughly parity |
| add | 1.2–2× | ~1× | memory-bound; both sides allocate the result |
| fma | 1.5–3× | 1.5–2× | base allocates two temporaries |
| pmax | 3–10× | 2–4× | base `pmax` is slow (attribute handling, recycling) |
| exp | 3–8× | 3–8× | SLEEF `u10` vs libm scalar |
| any_na | 2–6× | 2–4× | early exit, so small-n results depend on NA position |
| is_na | 2–4× | 2–3× | payload-aware mask |
| as_integer | 2–4× | 2–3× | base checks range and warns |

On arm64 (`neon`, 128-bit) expect roughly half these ratios for compute-bound ops and parity
for memory-bound ones. `avx512` over `avx2` typically gives 1.2–1.8× on compute-bound ops
and about 1× on memory-bound ones.

### Noise

Timings move between runs by 5–15% (more at n = 1e3) because of CPU frequency scaling and
turbo, thermal state, other processes, and where the input happens to sit in memory. For
numbers you want to keep, close other programs and run on mains power. Prefer the full run
to `--quick`, and repeat a run before acting on a difference. GitHub's hosted runners share
hosts with other jobs and can land on different CPU models from one week to the next. Treat
their results as rough, and compare only runs whose metadata shows the same CPU.

## Comparing two runs

```sh
Rscript bench/compare.R bench/results/OLD.rds bench/results/NEW.rds
Rscript bench/compare.R OLD.rds NEW.rds --threshold 10 --md   # markdown, 10% threshold
```

For every op, type, size and implementation present in both runs, `compare.R` prints the
old and new medians and `ratio = old / new` (above 1 means the new run is faster). It marks
changes beyond the threshold (default 15%) with ▲ faster or ▼ slower, and ends with a
count. It refuses to compare runs whose machine key differs (OS, architecture, CPU model
and detected CPU features) unless you pass `--force`. Whatever the timings, it exits with
status 0.

## Continuous integration

The `benchmarks` workflow (`.github/workflows/benchmarks.yaml`) runs the suite every
Saturday, and on demand from "Run workflow" with optional `quick`, `ops` and `sizes`
inputs. It runs on `ubuntu-latest` (x86-64), `ubuntu-24.04-arm` and `macos-14`. Each job
uploads `bench/results/` as an artifact kept for 90 days and appends `latest.md` to the job
summary. The workflow never commits results.

## Refreshing the benchmarks vignette

The benchmarks vignette does not run `bench`. It renders snapshots committed under
`vignettes/benchmark-results/`. To refresh them:

1. Run the full suite (`Rscript bench/run.R`) on a quiet machine, or download a
   `benchmarks` workflow artifact.
2. Copy the run's `.md` and `.csv` into `vignettes/benchmark-results/`, named after the
   platform (for example `linux-x86_64.md` and `linux-x86_64.csv`), replacing that
   platform's previous snapshot.
3. Rebuild the vignettes and check the tables and the metadata (date, rsimd version, CPU).

## Smoke scripts

The other scripts are quicker, narrower checks from when each function family was written.
Each prints a small table, per tier and against base R:

| Script | Covers |
|--------|--------|
| `sum.R` | `simd_sum` on 1e7 doubles, integers and integer64 (against bit64 when installed), in every precision mode |
| `reductions.R` | `simd_min`, `simd_any_na`, `simd_count_na`, and `na_check = FALSE` against the default |
| `linalg.R` | `simd_dot` and `simd_cumsum` |
| `arith.R` | `simd_add`, `simd_fma`, `simd_pmax`, and `simd_idiv` with vector and scalar divisors |
| `math.R` | `simd_exp`, `simd_log`, `simd_pow`, `simd_sin`, `simd_tanh` |
| `ml.R` | `simd_sigmoid`, `simd_softmax`, `simd_log_softmax` at three sizes |
| `complex.R` | the complex functions |

Run them with `Rscript bench/<script>`.
