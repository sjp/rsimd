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
# Modulus 1, so products stay finite.
u <- complex(modulus = 1, argument = stats::runif(n, -pi, pi))

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
  add(op = "mul", impl = tier, median_ms = time_ms(simd_mul(x, y)))
  add(op = "mul scalar", impl = tier, median_ms = time_ms(simd_mul(x, 2 - 1i)))
  add(op = "div", impl = tier, median_ms = time_ms(simd_div(x, y)))
  add(op = "div real", impl = tier, median_ms = time_ms(simd_div(x, 3)))
  add(op = "abs", impl = tier, median_ms = time_ms(simd_abs(x)))
  add(op = "eq", impl = tier, median_ms = time_ms(simd_eq(x, y)))
  add(op = "cumprod", impl = tier, median_ms = time_ms(simd_cumprod(u)))
  add(op = "sqrt", impl = tier, median_ms = time_ms(simd_sqrt(x)))
  add(op = "exp", impl = tier, median_ms = time_ms(simd_exp(x)))
  add(op = "log", impl = tier, median_ms = time_ms(simd_log(x)))
  add(op = "sin", impl = tier, median_ms = time_ms(simd_sin(x)))
  add(op = "asin", impl = tier, median_ms = time_ms(simd_asin(x)))
  add(op = "pow 3", impl = tier, median_ms = time_ms(simd_pow(x, 3)))
  add(op = "pow 0.5+0.5i", impl = tier, median_ms = time_ms(simd_pow(x, 0.5 + 0.5i)))
  for (mode in c("fast", "pairwise", "compensated")) {
    simd_precision(mode)
    add(op = paste("sum", mode), impl = tier, median_ms = time_ms(simd_sum(x)))
    add(op = paste("prod", mode), impl = tier, median_ms = time_ms(simd_prod(u)))
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
add(op = "mul", impl = "base x * y", median_ms = time_ms(x * y))
add(op = "mul scalar", impl = "base x * (2-1i)", median_ms = time_ms(x * (2 - 1i)))
add(op = "div", impl = "base x / y", median_ms = time_ms(x / y))
add(op = "div real", impl = "base x / 3", median_ms = time_ms(x / 3))
add(op = "abs", impl = "base Mod", median_ms = time_ms(Mod(x)))
add(op = "eq", impl = "base x == y", median_ms = time_ms(x == y))
add(op = "cumprod", impl = "base cumprod", median_ms = time_ms(cumprod(u)))
add(op = "prod", impl = "base prod", median_ms = time_ms(prod(u)))
add(op = "sqrt", impl = "base sqrt", median_ms = time_ms(sqrt(x)))
add(op = "exp", impl = "base exp", median_ms = time_ms(exp(x)))
add(op = "log", impl = "base log", median_ms = time_ms(log(x)))
add(op = "sin", impl = "base sin", median_ms = time_ms(sin(x)))
add(op = "asin", impl = "base asin", median_ms = time_ms(asin(x)))
add(op = "pow 3", impl = "base x^3", median_ms = time_ms(x^3))
add(op = "pow 0.5+0.5i", impl = "base x^(0.5+0.5i)", median_ms = time_ms(x^(0.5 + 0.5i)))

print(do.call(rbind, rows$list), row.names = FALSE)
