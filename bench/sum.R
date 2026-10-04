# Smoke benchmark of simd_sum() against base sum() for 1e7 doubles,
# integers and bit64 integer64 values, once per available implementation and
# precision mode. A sanity check that the SIMD tiers beat the scalar one, not
# a benchmark suite. integer64 is compared with bit64's sum() when bit64 is
# installed.
#
# Usage: Rscript bench/sum.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
inputs <- list(
  double = stats::runif(n), integer = sample.int(1000L, n, replace = TRUE),
  integer64 = simd_as_integer64(round(stats::runif(n, -2^40, 2^40)))
)

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
  if (type == "integer64") {
    if (!requireNamespace("bit64", quietly = TRUE)) next
    loadNamespace("bit64")
    m <- bench::mark(sum(x), iterations = 20)
    impl <- "bit64::sum"
  } else {
    m <- bench::mark(sum(x), iterations = 20)
    impl <- "base::sum"
  }
  rows[[length(rows) + 1L]] <- data.frame(
    type = type, mode = "-", impl = impl, median_ms = as.numeric(m$median) * 1e3
  )
}
simd_use("auto")

res <- do.call(rbind, rows)
res$median_ms <- round(res$median_ms, 2)
print(res, row.names = FALSE)
