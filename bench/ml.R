# Smoke benchmark of simd_sigmoid(), simd_softmax() and simd_log_softmax()
# against the base R one-liners, for n = 1e3, 1e5 and 1e7 doubles, once per
# available implementation: a sanity check that the SIMD tiers beat the
# scalar one and base R, not a benchmark suite.
#
# Usage: Rscript bench/ml.R   (with rsimd and bench installed)

library(rsimd)

# Median time of `expr` in milliseconds (the expression is re-evaluated on
# every iteration).
time_ms <- function(expr, iterations) {
  call <- substitute(bench::mark(expr, iterations = iterations, check = FALSE))
  m <- eval(call, parent.frame())
  round(as.numeric(m$median) * 1e3, 4)
}

base_sigmoid <- function(x) 1 / (1 + exp(-x))
base_softmax <- function(x) exp(x - max(x)) / sum(exp(x - max(x)))
base_log_softmax <- function(x) x - max(x) - log(sum(exp(x - max(x))))

rows <- new.env()
rows$list <- list()
add <- function(op, n, impl, ms) {
  rows$list[[length(rows$list) + 1L]] <- data.frame(op = op, n = n, impl = impl, median_ms = ms)
}
for (n in c(1e3, 1e5, 1e7)) {
  set.seed(1)
  x <- stats::runif(n, -10, 10)
  it <- if (n >= 1e7) 10 else 100
  for (tier in simd_available()) {
    simd_use(tier)
    add("sigmoid", n, tier, time_ms(simd_sigmoid(x), it))
    add("softmax", n, tier, time_ms(simd_softmax(x), it))
    add("log_softmax", n, tier, time_ms(simd_log_softmax(x), it))
  }
  simd_use("auto")
  add("sigmoid", n, "base", time_ms(base_sigmoid(x), it))
  add("sigmoid", n, "plogis", time_ms(stats::plogis(x), it))
  add("softmax", n, "base", time_ms(base_softmax(x), it))
  add("log_softmax", n, "base", time_ms(base_log_softmax(x), it))
}

res <- do.call(rbind, rows$list)
key <- paste(res$op, res$n)
base <- res$impl == "base"
res$speedup_vs_base <- round(res$median_ms[base][match(key, key[base])] / res$median_ms, 2)
print(res, row.names = FALSE)
