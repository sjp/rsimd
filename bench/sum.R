# Smoke benchmark of simd_sum() against base sum() for 1e7 doubles and
# integers, once per available implementation and precision mode. A sanity
# check that the SIMD tiers beat the scalar one, not a benchmark suite.
#
# Usage: Rscript bench/sum.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
inputs <- list(double = stats::runif(n), integer = sample.int(1000L, n, replace = TRUE))

rows <- list()
for (type in names(inputs)) {
  x <- inputs[[type]]
  modes <- if (type == "double") c("fast", "pairwise", "compensated") else "fast"
  for (mode in modes) {
    old <- simd_precision(mode)
    for (tier in simd_available()) {
      simd_use(tier)
      m <- bench::mark(simd_sum(x), iterations = 20, check = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        type = type, mode = mode, impl = tier, median_ms = as.numeric(m$median) * 1e3
      )
    }
    simd_precision(old)
  }
  m <- bench::mark(sum(x), iterations = 20)
  rows[[length(rows) + 1L]] <- data.frame(
    type = type, mode = "-", impl = "base::sum", median_ms = as.numeric(m$median) * 1e3
  )
}
simd_use("auto")

res <- do.call(rbind, rows)
res$median_ms <- round(res$median_ms, 2)
print(res, row.names = FALSE)
