# Smoke benchmark of simd_min() and simd_any_na() against base R for 1e7
# doubles and integers, once per available implementation, plus the integer
# simd_min() with na_check = FALSE against the checked default. A sanity
# check that the SIMD tiers beat the scalar one and that skipping the NA
# check costs nothing, not a benchmark suite. (simd_sum() is in sum.R.)
#
# Usage: Rscript bench/reductions.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
inputs <- list(double = stats::runif(n), integer = sample.int(1000L, n, replace = TRUE))

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
  x <- inputs[[type]]
  for (tier in simd_available()) {
    simd_use(tier)
    add(op = "min", type = type, impl = tier, median_ms = time_ms(simd_min(x)))
    add(op = "any_na", type = type, impl = tier, median_ms = time_ms(simd_any_na(x)))
    if (type == "integer") {
      add(
        op = "min, na_check = FALSE", type = type, impl = tier,
        median_ms = time_ms(simd_min(x, na_check = FALSE))
      )
    }
  }
  add(op = "min", type = type, impl = "base::min", median_ms = time_ms(min(x)))
  add(op = "any_na", type = type, impl = "base::anyNA", median_ms = time_ms(anyNA(x)))
}
simd_use("auto")

print(do.call(rbind, rows$list), row.names = FALSE)
