# S11: normalising importance (log-)weights and their effective sample
# size, the step every particle-filter iteration runs; S11b: systematic
# resampling; and the bootstrap particle filter for a stochastic-volatility
# model that uses both. Full demo: Rscript bench/scenarios/10-smc.R
#
# Sample:
#   s <- source("bench/scenarios/10-smc.R")$value
#   d <- s[[3]]$setup(1e4); s[[3]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# x_t = phi x_{t-1} + sigma e_t, y_t ~ N(0, exp(x_t)); weights from the
# log-density of y_t given x_t (up to a constant): -x/2 - y^2 exp(-x) / 2.
sv_simulate <- function(steps, phi = 0.95, sigma = 0.25) {
  x <- numeric(steps)
  x[1] <- rnorm(1, 0, sigma / sqrt(1 - phi^2))
  for (t in 2:steps) x[t] <- phi * x[t - 1] + rnorm(1, 0, sigma)
  rnorm(steps) * exp(x / 2)
}

# Systematic resampling: the indices of the cumulative weights reached by
# (u + 0:(N - 1)) / N.
resample_base <- function(w, u) findInterval((u + seq_along(w) - 1) / length(w), cumsum(w)) + 1L
resample_rsimd <- function(w, u) {
  findInterval((u + seq_along(w) - 1) / length(w), simd_cumsum(w)) + 1L
}

# The filter's log-likelihood estimate. The propagation and resampling
# draws are the same on both sides (seeded inside); `step` computes the
# log-weights, normalised weights, ESS and log of the mean weight.
particle_filter <- function(y, n, step, resample, phi = 0.95, sigma = 0.25) {
  set.seed(7)
  x <- rnorm(n, 0, sigma / sqrt(1 - phi^2))
  ll <- 0
  for (t in seq_along(y)) {
    x <- phi * x + sigma * rnorm(n)
    s <- step(x, y[t])
    ll <- ll + s$lmean
    if (s$ess < n / 2) x <- x[resample(s$w, runif(1))]
  }
  ll
}
pf_step_base <- function(x, yt) {
  lw <- -0.5 * x - 0.5 * yt^2 * exp(-x)
  m <- max(lw)
  e <- exp(lw - m)
  se <- sum(e)
  w <- e / se
  list(w = w, ess = 1 / sum(w * w), lmean = m + log(se / length(x)))
}
pf_step_rsimd <- function(x, yt) {
  lw <- simd_mul_add(simd_exp(simd_neg(x)), -0.5 * yt^2, simd_mul(x, -0.5))
  w <- simd_softmax(lw)
  # The log of the mean weight from the softmax normaliser, lw[k] - log w[k]
  # at the largest weight (which cannot underflow).
  k <- simd_which_max(lw)
  list(w = w, ess = 1 / simd_sum_sq(w), lmean = lw[k] - log(w[k]) - log(length(x)))
}

pf <- scenario(
  id = "particle_filter", family = "montecarlo", flagship = TRUE,
  title = "Bootstrap particle filter, stochastic volatility, 100 steps (n = particles)",
  why = "Every SMC step normalises weights, checks the ESS and resamples; nothing in base R does it.",
  sizes = c(1e3, 1e4, 1e5),
  quick_sizes = c(1e3, 1e4),
  setup = function(n) list(y = sv_simulate(100), n = n),
  base = function(d) particle_filter(d$y, d$n, pf_step_base, resample_base),
  rsimd = function(d) particle_filter(d$y, d$n, pf_step_rsimd, resample_rsimd),
  # The two sides resample from weights equal to rounding; an index that a
  # rounding difference moves changes the path, so agreement is to 1e-3.
  tolerance = 1e-3
)

out <- list(
  scenario(
    id = "smc_weights", family = "montecarlo", flagship = TRUE,
    title = "Normalise log-weights and compute the effective sample size",
    why = "The inner step of importance sampling and every particle filter.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) list(lw = rnorm(n, -1000, 5)),
    base = function(d) {
      w <- exp(d$lw - max(d$lw))
      w <- w / sum(w)
      c(w, ess = 1 / sum(w^2))
    },
    rsimd = function(d) {
      w <- simd_softmax(d$lw)
      c(w, ess = 1 / simd_sum_sq(w))
    },
    precision = c("fast", "pairwise"),
    tolerance = c(fast = 1e-12, pairwise = 1e-13)
  ),
  scenario(
    id = "resample", family = "montecarlo",
    title = "Systematic resampling (cumulative weights + findInterval)",
    why = "findInterval() is the same on both sides and limits the gain.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) {
      w <- rexp(n)
      list(w = w / sum(w), u = runif(1))
    },
    base = function(d) resample_base(d$w, d$u),
    rsimd = function(d) resample_rsimd(d$w, d$u),
    # The share of indices that differ (a rounding difference in the
    # cumulative sum can move a boundary).
    error = function(got, want) mean(got != want),
    tolerance = 1e-4
  ),
  pf
)

# Full demo: T = 500 steps, 1e5 particles.
if (sys.nframe() == 0L) {
  set.seed(20261006)
  y <- sv_simulate(500)
  tb <- system.time(lb <- particle_filter(y, 1e5, pf_step_base, resample_base))[["elapsed"]]
  tr <- system.time(lr <- particle_filter(y, 1e5, pf_step_rsimd, resample_rsimd))[["elapsed"]]
  cat(sprintf(
    "Particle filter, T = 500, N = 1e5: base %.2f s, rsimd (%s) %.2f s, %.1fx; loglik %.4f / %.4f\n",
    tb, simd_current(), tr, tb / tr, lb, lr
  ))
}

invisible(out)
