# Smoke benchmark of simd_dot() and simd_cumsum() against base R for 1e7
# doubles and integers, once per available implementation: a sanity check
# that the SIMD tiers beat the scalar one, not a benchmark suite. base R's
# dot product is sum(x * y), which also allocates x * y.
#
# Usage: Rscript bench/linalg.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
inputs <- list(
  double = list(x = stats::runif(n), y = stats::runif(n)),
  integer = list(x = sample.int(100L, n, replace = TRUE), y = sample.int(100L, n, replace = TRUE))
)

# Median time of `expr` in milliseconds (the expression is re-evaluated on
# every iteration).
time_ms <- function(expr) {
  call <- substitute(bench::mark(expr, iterations = 20, check = FALSE))
  m <- eval(call, parent.frame())
  round(as.numeric(m$median) * 1e3, 2)
}

rows <- new.env()
rows$list <- list()
add <- function(...) rows$list[[length(rows$list) + 1L]] <- data.frame(...)
for (type in names(inputs)) {
  x <- inputs[[type]]$x
  y <- inputs[[type]]$y
  for (tier in simd_available()) {
    simd_use(tier)
    add(op = "dot", type = type, impl = tier, median_ms = time_ms(simd_dot(x, y)))
    add(op = "cumsum", type = type, impl = tier, median_ms = time_ms(simd_cumsum(x)))
  }
  add(op = "dot", type = type, impl = "base sum(x * y)", median_ms = time_ms(sum(x * y)))
  add(op = "cumsum", type = type, impl = "base::cumsum", median_ms = time_ms(cumsum(x)))
}
simd_use("auto")

print(do.call(rbind, rows$list), row.names = FALSE)
