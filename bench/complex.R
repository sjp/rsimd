# Smoke benchmark of the complex functions against base R for 1e7 complex
# numbers, once per available implementation: a sanity check that the SIMD
# tiers beat the scalar one, not a benchmark suite.
#
# Usage: Rscript bench/complex.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
x <- complex(real = stats::runif(n), imaginary = stats::runif(n))
y <- complex(real = stats::runif(n), imaginary = stats::runif(n))

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
for (tier in simd_available()) {
  simd_use(tier)
  add(op = "add", impl = tier, median_ms = time_ms(simd_add(x, y)))
  add(op = "add scalar", impl = tier, median_ms = time_ms(simd_add(x, 1i)))
  add(op = "neg", impl = tier, median_ms = time_ms(simd_neg(x)))
  add(op = "conj", impl = tier, median_ms = time_ms(simd_conj(x)))
  add(op = "re", impl = tier, median_ms = time_ms(simd_re(x)))
  add(op = "is_na", impl = tier, median_ms = time_ms(simd_is_na(x)))
  for (mode in c("fast", "pairwise", "compensated")) {
    simd_precision(mode)
    add(op = paste("sum", mode), impl = tier, median_ms = time_ms(simd_sum(x)))
  }
  simd_precision("fast")
}
simd_use("auto")
add(op = "add", impl = "base x + y", median_ms = time_ms(x + y))
add(op = "add scalar", impl = "base x + 1i", median_ms = time_ms(x + 1i))
add(op = "neg", impl = "base -x", median_ms = time_ms(-x))
add(op = "conj", impl = "base Conj", median_ms = time_ms(Conj(x)))
add(op = "re", impl = "base Re", median_ms = time_ms(Re(x)))
add(op = "is_na", impl = "base is.na", median_ms = time_ms(is.na(x)))
add(op = "sum", impl = "base sum", median_ms = time_ms(sum(x)))

print(do.call(rbind, rows$list), row.names = FALSE)
