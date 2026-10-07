# S2: winsorized mean and variance: values beyond the limits are set to
# the limits, not dropped.
#
# Sample:
#   s <- source("bench/scenarios/02-winsorized.R")$value
#   d <- s[[1]]$setup(1e5); s[[1]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

winsor_setup <- function(n) {
  x <- rgamma(n, 2)
  lim <- quantile(x, c(0.05, 0.95), names = FALSE)
  list(x = x, lo = lim[1], hi = lim[2])
}

invisible(list(
  scenario(
    id = "winsorized", family = "descriptive", flagship = TRUE,
    title = "Winsorized mean and variance (limits given)",
    why = "Base R trims (mean(trim =)) but does not winsorize.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = winsor_setup,
    base = function(d) {
      w <- pmin(pmax(d$x, d$lo), d$hi)
      c(mean(w), var(w))
    },
    rsimd = function(d) {
      w <- simd_clamp(d$x, d$lo, d$hi)
      c(simd_mean(w), simd_var(w))
    },
    tolerance = 1e-12
  ),
  # The limits are usually quantiles; quantile() sorts, and that cost is
  # the same on both sides.
  scenario(
    id = "winsorized_q", family = "descriptive",
    title = "Winsorized mean and variance (5% and 95% quantiles as limits)",
    why = "As winsorized, with the quantile() call both sides make.",
    sizes = c(1e3, 1e5, 1e6),
    setup = function(n) list(x = rgamma(n, 2)),
    base = function(d) {
      lim <- quantile(d$x, c(0.05, 0.95), names = FALSE)
      w <- pmin(pmax(d$x, lim[1]), lim[2])
      c(mean(w), var(w))
    },
    rsimd = function(d) {
      lim <- quantile(d$x, c(0.05, 0.95), names = FALSE)
      w <- simd_clamp(d$x, lim[1], lim[2])
      c(simd_mean(w), simd_var(w))
    },
    tolerance = 1e-12
  )
))
