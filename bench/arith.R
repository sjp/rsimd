# Smoke benchmark of simd_add(), simd_fma(), simd_pmax() and simd_idiv()
# against base R for 1e7 doubles and integers, once per available
# implementation: a sanity check that the SIMD tiers beat the scalar one,
# not a benchmark suite. base R has no fused multiply-add; its row is
# x * y + z. simd_idiv() is timed with a vector divisor and with a scalar
# one (integers take a division-free path then). Base R's double %/% uses
# long double, which is very slow where that is software quad precision
# (arm64 Linux), so it is timed with fewer iterations.
#
# Usage: Rscript bench/arith.R   (with rsimd and bench installed)

library(rsimd)

n <- 1e7
set.seed(1)
inputs <- list(
  double = list(x = stats::runif(n), y = stats::runif(n), z = stats::runif(n)),
  integer = list(
    x = sample.int(1000L, n, replace = TRUE), y = sample.int(1000L, n, replace = TRUE),
    z = sample.int(1000L, n, replace = TRUE)
  )
)

# Median time of `expr` in milliseconds (the expression is re-evaluated on
# every iteration).
time_ms <- function(expr, iterations = 20) {
  call <- substitute(bench::mark(expr, iterations = iterations, check = FALSE))
  m <- eval(call, parent.frame())
  round(as.numeric(m$median) * 1e3, 2)
}

rows <- new.env()
rows$list <- list()
add <- function(...) rows$list[[length(rows$list) + 1L]] <- data.frame(...)
for (type in names(inputs)) {
  x <- inputs[[type]]$x
  y <- inputs[[type]]$y
  z <- inputs[[type]]$z
  d <- if (type == "double") 7 else 7L
  for (tier in simd_available()) {
    simd_use(tier)
    add(op = "add", type = type, impl = tier, median_ms = time_ms(simd_add(x, y)))
    if (type == "double") {
      add(op = "fma", type = type, impl = tier, median_ms = time_ms(simd_fma(x, y, z)))
    } else {
      add(op = "mul_add", type = type, impl = tier, median_ms = time_ms(simd_mul_add(x, y, z)))
    }
    add(op = "pmax", type = type, impl = tier, median_ms = time_ms(simd_pmax(x, y)))
    add(op = "idiv", type = type, impl = tier, median_ms = time_ms(simd_idiv(x, y)))
    add(op = "idiv by 7", type = type, impl = tier, median_ms = time_ms(simd_idiv(x, d)))
  }
  add(op = "add", type = type, impl = "base x + y", median_ms = time_ms(x + y))
  add(
    op = if (type == "double") "fma" else "mul_add", type = type, impl = "base x * y + z",
    median_ms = time_ms(x * y + z)
  )
  add(op = "pmax", type = type, impl = "base pmax", median_ms = time_ms(pmax(x, y)))
  base_iter <- if (type == "double") 3 else 20
  add(op = "idiv", type = type, impl = "base %/%", median_ms = time_ms(x %/% y, base_iter))
  add(op = "idiv by 7", type = type, impl = "base %/%", median_ms = time_ms(x %/% d, base_iter))
}
simd_use("auto")

print(do.call(rbind, rows$list), row.names = FALSE)
