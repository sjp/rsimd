# Smoke benchmark of simd_exp(), simd_log(), simd_pow(), simd_sin() and
# simd_tanh() against base R for 1e6 doubles, once per available
# implementation: a sanity check that the SIMD tiers beat the scalar one and
# base R, not a benchmark suite.
#
# Usage: Rscript bench/math.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e6
set.seed(1)
x <- stats::runif(n, -10, 10)
pos <- stats::runif(n, 1e-3, 1e3)
y <- stats::runif(n, -3, 3)

# Median time of `expr` in milliseconds (the expression is re-evaluated on
# every iteration).
time_ms <- function(expr) {
  call <- substitute(bench::mark(expr, iterations = 30, check = FALSE))
  m <- eval(call, parent.frame())
  round(as.numeric(m$median) * 1e3, 3)
}

rows <- new.env()
rows$list <- list()
add <- function(op, impl, ms) {
  rows$list[[length(rows$list) + 1L]] <- data.frame(op = op, impl = impl, median_ms = ms)
}
for (tier in simd_available()) {
  simd_use(tier)
  add("exp", tier, time_ms(simd_exp(x)))
  add("log", tier, time_ms(simd_log(pos)))
  add("pow", tier, time_ms(simd_pow(pos, y)))
  add("sin", tier, time_ms(simd_sin(x)))
  add("tanh", tier, time_ms(simd_tanh(x)))
}
simd_use("auto")
add("exp", "base", time_ms(exp(x)))
add("log", "base", time_ms(log(pos)))
add("pow", "base", time_ms(pos^y))
add("sin", "base", time_ms(sin(x)))
add("tanh", "base", time_ms(tanh(x)))

res <- do.call(rbind, rows$list)
base_ms <- res$median_ms[res$impl == "base"][match(res$op, res$op[res$impl == "base"])]
res$speedup_vs_base <- round(base_ms / res$median_ms, 2)
print(res, row.names = FALSE)
