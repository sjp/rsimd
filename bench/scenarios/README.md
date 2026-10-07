# Statistical scenarios

`bench/run.R` times single primitives against their base R equivalents (`simd_sum()`
against `sum()`). The scenarios here answer a different question: **which statistical
calculations that base R has no function for get faster when written with rsimd's
primitives?** Each scenario is a short computation a statistician writes by hand, such as
skewness, a weighted variance, a Box-Cox profile likelihood or a particle filter step,
written twice: as base R code and with `simd_*` calls.

The numbers are informational. They vary between machines, and nothing checks them.

## Running

```sh
Rscript bench/scenarios/run.R --quick        # each scenario's two smallest sizes, ~30 s
Rscript bench/scenarios/run.R                # all sizes, a few minutes
Rscript bench/scenarios/run.R --only moments,boxcox
Rscript bench/scenarios/run.R --family ml    # descriptive, likelihood, ml, montecarlo,
                                             # kernel, distance, rolling, accuracy
Rscript bench/scenarios/run.R --flagship --out /tmp/scen
```

The runner needs rsimd and bench. It writes
`bench/results/scenarios-<time>-<machine>-<sha>.{rds,csv,md}` and
`bench/results/scenarios-latest.md` (git ignores `bench/results/`). The `.md` file opens
with the same machine description as `bench/run.R`'s results. It then has a flagship
summary, and for each scenario a table with one row per n and these columns:

- `base`: the fair base R code (see the rules below).
- `idiomatic`: the usual base R form, where it differs from the fair one.
- Reference implementations from CRAN packages, where installed (matrixStats, moments,
  MASS, data.table, ...). Their errors are shown but not checked, because a package may
  define the statistic slightly differently.
- One `rsimd` column per tier in `simd_available()`, `none` included. When a scenario
  sweeps the precision modes, each column is `tier/mode`.

Each cell shows the median time and the speedup over `base`, the memory allocated per call,
and the error against `base`, or against the exact answer where the scenario has one. The
error is normwise relative (`max|got - want| / max|want|`) unless the scenario defines its
own. Each table also states the crossover n: the smallest size at which the best tier beats
base. The file ends with the cost of composing the primitives that rsimd lacks (see
[Gaps](#gaps)).

If a column's error exceeds the scenario's tolerance, the runner exits with status 1 and
lists the cells. Those numbers must not be quoted.

On arm64 Linux, base R's `long double` is 128-bit IEEE and is computed in software, so base
`sum()`, `mean()`, `var()` and `cumsum()` are much slower there than on x86-64. Speedups of
reduction-heavy scenarios are inflated on that platform. Quote x86-64 numbers, which the
`benchmarks` workflow produces with `suite: scenarios`.

## Each scenario is a sample

Every `NN-*.R` file runs on its own with only rsimd and base R installed. It loads rsimd and
defines a minimal `scenario()` when `common.R` has not been sourced. Its value is the
scenario, or a list of scenarios:

```r
s <- source("bench/scenarios/01-moments.R")$value
d <- s$setup(1e5)
s$rsimd(d)   # the rsimd version
s$base(d)    # the same computation in base R
```

Three files also run a full-scale demo with `Rscript`:

- `07-boxcox.R`: the Box-Cox profile log-likelihood at n = 1e6.
- `10-smc.R`: a particle filter with 500 steps and 1e5 particles.
- `13-tanimoto.R`: Tanimoto search over databases of 1e3 to 1e5 fingerprints.

## Fairness rules

1. **Same algebra on both sides.** When the rsimd version uses a cheaper formulation, the
   `base` column uses it too, and the idiomatic form is shown as an extra column. Examples:
   a clamp instead of `ifelse()`, `var(x^λ)/λ²` instead of `var((x^λ - 1)/λ)`, and `d2 * d`
   instead of `d^3`.
2. **Vectorised base code.** Base R loops only where a competent user would loop: per point,
   for O(n²) statistics whose n × n matrix would not fit. The `outer()` form is shown too
   where it fits.
3. **Results must agree.** Every cell's error is checked against the scenario's tolerance.
4. **The platform is stated** in every results file (see the caveat above).
5. **No timing assertions.** Timings are never checked.

## Scenarios

★ marks a flagship. Controls are scenarios where rsimd is not the right tool; they are kept
to show that.

| File | Id | What |
|---|---|---|
| 01 | `moments` ★ | Skewness and kurtosis (ref: moments) |
| 02 | `winsorized` ★, `winsorized_q` | Winsorized mean and variance, without and with the quantile() call |
| 03 | `weighted_var`, `gini`, `gini_sort` | Weighted mean and variance (ref: matrixStats); Gini coefficient, without and with the sort |
| 04 | `means` | Geometric and harmonic means, coefficient of variation |
| 05 | `logistic_nll`, `logistic_optim` | Logistic regression negative log-likelihood, once and inside an `optim()` fit |
| 06 | `metropolis` | Metropolis sampler for a Student-t location and scale, 200 iterations |
| 07 | `boxcox` ★ | Box-Cox profile log-likelihood over 41 values of λ (ref: MASS) |
| 08 | `logsumexp` | log-sum-exp (ref: matrixStats) |
| 09 | `logloss`, `kl`, `huber` | Log loss and Brier score; KL divergence and entropy; Huber loss and IRLS weights |
| 10 | `smc_weights` ★, `resample`, `particle_filter` ★ | Normalise log-weights and compute the ESS; systematic resampling; a bootstrap particle filter |
| 11 | `kde_cv`, `rbf_matrix` (control) | KDE leave-one-out likelihood CV; an RBF kernel matrix |
| 12 | `energy` | Energy distance between two samples |
| 13 | `tanimoto` ★ | Tanimoto search over 2048-bit fingerprints |
| 14 | `knn_columns`, `knn_blas` (controls) | kNN distances by a loop over columns, and by the BLAS identity |
| 15 | `rolling`, `auc` | Rolling mean and sd (refs: stats::filter, data.table); trapezoid AUC |
| 16 | `sum_cancel`, `var_offset` ★ | Accuracy at scale, across the precision modes |

## Gaps

Some scenarios have to compose several primitives where one fused primitive would make a
single pass. The results file times each of these compositions against the passes a fused
primitive would make, timed as `simd_sum()` (memory traffic only, no arithmetic):

| Gap | Missing primitive | Scenarios |
|---|---|---|
| G1 | A vector applied against each column of a matrix (`X[, j]` copies) | `knn_columns` |
| G2 | Keeping `dim` on elementwise results | `rbf_matrix` |
| G3 | A fused log-sum-exp | `logsumexp` |
| G4 | Central moments Σ(x − c)ᵏ without temporaries | `moments` |
| G5 | A weighted variance | `weighted_var` |
| G6 | Reductions of a transform (Σ exp, Σ log, Σ log1p, exp(a·x)) without the temporary | `means`, `logistic_nll`, `metropolis`, `kde_cv` |
| G8 | Rolling-window statistics (Welford) and exponentially weighted scans | `rolling` |
| G9 | Offsets or differences without copying slices | `auc` |
| G10 | The population count of an AND, per block | `tanimoto` |

Whether to add any of these primitives is a separate decision; the scenarios only measure
what composing them costs today.
