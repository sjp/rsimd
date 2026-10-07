# S5: the negative log-likelihood of a logistic regression, once for given
# linear predictors and inside an optim() fit.
#
# Sample:
#   s <- source("bench/scenarios/05-logistic.R")$value
#   d <- s[[1]]$setup(1e5); s[[1]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# The stable form on both sides: log(1 + exp(eta)) = max(eta, 0) + log1p(exp(-|eta|)).
logistic_nll_base <- function(eta, y) {
  sum(pmax(eta, 0) + log1p(exp(-abs(eta)))) - sum(y * eta)
}
logistic_nll_rsimd <- function(eta, y) {
  simd_sum(simd_add(simd_pmax(eta, 0), simd_log1p(simd_exp(simd_neg(simd_abs(eta)))))) -
    simd_dot(y, eta)
}

logistic_setup <- function(n) {
  X <- cbind(1, matrix(rnorm(n * 4), n, 4))
  beta <- c(-0.5, 1, -1, 0.5, 0.25)
  eta <- drop(X %*% beta)
  list(X = X, eta = eta, y = as.numeric(runif(n) < plogis(eta)))
}

# A fixed budget of BFGS iterations; the gradient X'(p - y) is the same base
# R code on both sides, as is X %*% beta (BLAS).
logistic_fit <- function(d, nll) {
  fn <- function(b) nll(drop(d$X %*% b), d$y)
  gr <- function(b) drop(crossprod(d$X, plogis(drop(d$X %*% b)) - d$y))
  optim(rep(0, ncol(d$X)), fn, gr, method = "BFGS", control = list(maxit = 15))$par
}

invisible(list(
  scenario(
    id = "logistic_nll", family = "likelihood",
    title = "Logistic regression negative log-likelihood, linear predictor given",
    why = "Writing an MLE or MCMC by hand means evaluating this sum many times.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = logistic_setup,
    base = function(d) logistic_nll_base(d$eta, d$y),
    base_idiomatic = function(d) -sum(dbinom(d$y, 1, plogis(d$eta), log = TRUE)),
    rsimd = function(d) logistic_nll_rsimd(d$eta, d$y),
    tolerance = 1e-12,
    gap = list(id = "G6", ideal = function(d) simd_dot(d$y, d$eta))
  ),
  scenario(
    id = "logistic_optim", family = "likelihood",
    title = "Logistic regression fitted by optim() (15 BFGS iterations, 5 coefficients)",
    why = "The likelihood is one part of the work: X %*% beta (BLAS) and the gradient are the same on both sides.",
    sizes = c(1e3, 1e5),
    setup = logistic_setup,
    base = function(d) logistic_fit(d, logistic_nll_base),
    rsimd = function(d) logistic_fit(d, logistic_nll_rsimd),
    tolerance = 1e-8
  )
))
