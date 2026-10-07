# S9: log loss (with the Brier score), S9b: KL divergence (with Shannon
# entropy), S10: Huber loss and the IRLS weights of an M-estimator.
#
# Sample:
#   s <- source("bench/scenarios/09-losses.R")$value
#   d <- s[[1]]$setup(1e5); s[[1]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(list(
  scenario(
    id = "logloss", family = "ml",
    title = "Log loss of probabilities clamped to [eps, 1 - eps], and the Brier score",
    why = "Scoring a classifier's probabilities; base R has no loss functions.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) {
      p <- runif(n)
      list(p = p, y = as.numeric(runif(n) < p), eps = 1e-15)
    },
    base = function(d) {
      n <- length(d$p)
      q <- pmin(pmax(d$p, d$eps), 1 - d$eps)
      c(
        logloss = -(sum(d$y * log(q)) + sum((1 - d$y) * log1p(-q))) / n,
        brier = sum((d$p - d$y)^2) / n
      )
    },
    rsimd = function(d) {
      n <- length(d$p)
      q <- simd_clamp(d$p, d$eps, 1 - d$eps)
      c(
        logloss = -(simd_dot(d$y, simd_log(q)) +
          simd_dot(simd_sub(1, d$y), simd_log1p(simd_neg(q)))) / n,
        brier = simd_dist(d$p, d$y)^2 / n
      )
    },
    tolerance = 1e-12
  ),
  scenario(
    id = "kl", family = "ml",
    title = "KL divergence of two distributions, and the entropy of the first",
    why = "Base R has no divergence or entropy functions.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) {
      p <- runif(n)
      q <- runif(n)
      list(p = p / sum(p), q = q / sum(q))
    },
    base = function(d) c(kl = sum(d$p * log(d$p / d$q)), entropy = -sum(d$p * log(d$p))),
    rsimd = function(d) {
      c(
        kl = simd_dot(d$p, simd_log(simd_div(d$p, d$q))),
        entropy = -simd_dot(d$p, simd_log(d$p))
      )
    },
    tolerance = 1e-12
  ),
  # Huber loss sum(rho(r)) with rho(r) = r^2/2 for |r| <= k, k(|r| - k/2)
  # beyond, written as clamp(r) * (r - clamp(r)/2); the IRLS weights are
  # k / max(|r|, k).
  scenario(
    id = "huber", family = "ml",
    title = "Huber loss and IRLS weights (one M-estimation step)",
    why = "Robust regression by hand: MASS::rlm() hides the loss; the weights are needed per step.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) list(r = rt(n, 2), k = 1.345),
    base = function(d) {
      cr <- pmin(pmax(d$r, -d$k), d$k)
      c(loss = sum(cr * (d$r - cr * 0.5)), weights = d$k / pmax(abs(d$r), d$k))
    },
    base_idiomatic = function(d) {
      a <- abs(d$r)
      c(
        loss = sum(ifelse(a <= d$k, d$r^2 / 2, d$k * (a - d$k / 2))),
        weights = ifelse(a <= d$k, 1, d$k / a)
      )
    },
    rsimd = function(d) {
      cr <- simd_clamp(d$r, -d$k, d$k)
      c(
        loss = simd_dot(cr, simd_sub(d$r, simd_mul(cr, 0.5))),
        weights = simd_div(d$k, simd_pmax(simd_abs(d$r), d$k))
      )
    },
    tolerance = 1e-12
  )
))
