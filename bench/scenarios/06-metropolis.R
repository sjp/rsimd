# S6: a random-walk Metropolis sampler for the location and scale of a
# Student-t (3 df) sample. A Normal likelihood would reduce to sufficient
# statistics; this one needs a pass over the data per iteration.
#
# Sample:
#   s <- source("bench/scenarios/06-metropolis.R")$value
#   d <- s$setup(1e4); colMeans(s$rsimd(d))
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# log L(mu, s) up to a constant: -n log(s) - 2 sum(log1p(z^2 / 3)), z = (x - mu) / s.
t3_loglik_base <- function(x, mu, s) {
  z <- x * (1 / s) + (-mu / s)
  -length(x) * log(s) - 2 * sum(log1p(z * z * (1 / 3)))
}
t3_loglik_rsimd <- function(x, mu, s) {
  z <- simd_mul_add(x, 1 / s, -mu / s)
  -length(x) * log(s) - 2 * simd_sum(simd_log1p(simd_mul(simd_mul(z, z), 1 / 3)))
}

# The same proposals and uniforms on both sides (seeded inside), so the
# chains agree while the log-likelihoods agree to rounding.
metropolis <- function(x, loglik, iter = 200L) {
  set.seed(42)
  prop <- matrix(rnorm(2 * iter, sd = 0.01), iter, 2)
  u <- log(runif(iter))
  th <- c(median(x), log(mad(x)))
  ll <- loglik(x, th[1], exp(th[2]))
  out <- matrix(0, iter, 2)
  for (k in seq_len(iter)) {
    cand <- th + prop[k, ]
    llc <- loglik(x, cand[1], exp(cand[2]))
    if (u[k] < llc - ll) {
      th <- cand
      ll <- llc
    }
    out[k, ] <- th
  }
  out
}

invisible(scenario(
  id = "metropolis", family = "likelihood",
  title = "Metropolis sampler, Student-t location and scale, 200 iterations",
  why = "Hand-written MCMC evaluates a non-Gaussian likelihood over all data every iteration.",
  sizes = c(1e3, 1e4, 1e5),
  quick_sizes = c(1e3, 1e4),
  setup = function(n) list(x = 5 + 2 * rt(n, 3)),
  base = function(d) metropolis(d$x, t3_loglik_base),
  rsimd = function(d) metropolis(d$x, t3_loglik_rsimd),
  tolerance = 1e-10,
  # One pass over the data per likelihood evaluation (201 of them).
  gap = list(id = "G6", ideal = function(d) for (k in 1:201) simd_sum(d$x))
))
