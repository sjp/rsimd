# Checked and wrapping integer arithmetic, integer sums and %/% / %%,
# checked through the self-test kernels.

imax <- .Machine$integer.max
overflow_msg <- "NAs produced by integer overflow"

# Base R's checked result for op, with the warning muffled.
base_op <- function(op, x, y) {
  suppressWarnings(switch(op,
    add = x + y,
    sub = x - y,
    mul = x * y
  ))
}

test_that("checked add, sub and mul give NA and one warning, as base R", {
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_warning(r <- .debug_arith("add", imax, 1L), overflow_msg, fixed = TRUE)
      expect_identical(r, NA_integer_)
      # 1000 overflows, exactly one warning.
      w <- 0L
      r <- withCallingHandlers(.debug_arith("add", rep(imax, 1000L), 1:1000),
        warning = function(c) {
          w <<- w + 1L
          invokeRestart("muffleWarning")
        }
      )
      expect_identical(w, 1L, info = tier)
      expect_identical(r, rep(NA_integer_, 1000L))
      expect_warning(r <- .debug_arith("mul", 46341L, 46341L), overflow_msg, fixed = TRUE)
      expect_identical(r, NA_integer_)
      # -INT_MAX - 1 wraps to INT_MIN without two's complement overflow, but
      # INT_MIN is NA: base R warns, so does rsimd.
      expect_warning(r <- .debug_arith("add", -imax, -1L), overflow_msg, fixed = TRUE)
      expect_identical(r, NA_integer_)
      expect_warning(r <- .debug_arith("sub", -imax, 1L), overflow_msg, fixed = TRUE)
      expect_identical(r, NA_integer_)
      expect_warning(r <- .debug_arith("mul", -65536L, 32768L), overflow_msg, fixed = TRUE)
      expect_identical(r, NA_integer_)
      # NA operands give NA without a warning.
      expect_no_warning(r <- .debug_arith("add", c(NA, 1L), c(1L, NA)))
      expect_identical(r, c(NA_integer_, NA_integer_))
      # In-range results are exact.
      expect_identical(.debug_arith("add", imax - 1L, 1L), imax)
      expect_identical(.debug_arith("mul", 46340L, 46340L), 2147395600L)
      expect_identical(.debug_arith("mul", -46340L, 46340L), -2147395600L)
      expect_identical(.debug_arith("sub", -imax + 1L, 1L), -imax)
    })
  }
  expect_warning(base_check <- -imax - 1L, overflow_msg, fixed = TRUE)
  expect_identical(base_check, NA_integer_)
})

test_that("checked ops match base R on random operands", {
  set.seed(3)
  big <- function(n) as.integer(sample(c(-1, 1), n, TRUE) * stats::runif(n, 0, imax))
  small <- function(n) sample(-70000L:70000L, n, TRUE)
  for (n in c(1, 3, 4, 7, 8, 15, 16, 17, 33, 100)) {
    for (op in c("add", "sub", "mul")) {
      gen <- if (op == "mul") small else big
      x <- gen(n)
      y <- gen(n)
      x[seq_len(n) %% 9 == 4] <- NA
      want <- base_op(op, x, y)
      warns <- anyNA(want[!is.na(x) & !is.na(y)])
      res <- expect_tiers_identical(function() suppressWarnings(.debug_arith(op, x, y)),
        label = paste(op, n)
      )
      expect_identical(res$none, want, info = paste(op, n))
      for_each_tier(function(tier) {
        got <- tryCatch(.debug_arith(op, x, y), warning = function(w) "warned")
        expect_identical(identical(got, "warned"), warns, info = paste(op, n, tier))
      })
      # Scalar operands broadcast on either side.
      expect_tiers_identical(function() {
        suppressWarnings(list(.debug_arith(op, x[1], y), .debug_arith(op, x, y[1])))
      }, label = paste(op, n, "scalar"))
      expect_identical(suppressWarnings(.debug_arith(op, x, y[1])), base_op(op, x, y[1]))
    }
  }
})

test_that("wrapping ops wrap, honour NA and never warn", {
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      info <- tier
      expect_no_warning(r <- .debug_arith("add_wrap", c(imax, imax, NA), c(2L, 1L, 1L)))
      # imax + 1 wraps to INT_MIN, which reads as NA: the documented trap.
      expect_identical(r, c(-imax, NA, NA), info = info)
      expect_identical(.debug_arith("sub_wrap", -imax, 2L), imax, info = info)
      expect_identical(.debug_arith("mul_wrap", 65537L, 65537L), 131073L, info = info)
      expect_identical(.debug_arith("mul_wrap", 46341L, 46341L), -2147479015L, info = info)
      expect_identical(.debug_arith("neg_wrap", c(1L, -imax, NA)), c(-1L, imax, NA), info = info)
      expect_identical(.debug_arith("abs_wrap", c(-1L, -imax, NA)), c(1L, imax, NA), info = info)
    })
  }
  set.seed(5)
  x <- as.integer(stats::runif(50, -imax, imax))
  y <- as.integer(stats::runif(50, -imax, imax))
  x[c(3, 20)] <- NA
  for (op in c("add_wrap", "sub_wrap", "mul_wrap")) {
    expect_tiers_identical(function() .debug_arith(op, x, y), label = op)
  }
})

test_that("neg and abs map NA to NA", {
  expect_tiers_identical(function() {
    list(
      .debug_arith("neg", c(NA, 5L, -imax, imax, 0L)),
      .debug_arith("abs", c(NA, -5L, -imax, imax, 0L))
    )
  })
  expect_identical(.debug_arith("neg", NA_integer_), NA_integer_)
  expect_identical(.debug_arith("abs", NA_integer_), NA_integer_)
  expect_identical(.debug_arith("neg", c(5L, -imax)), c(-5L, imax))
})

test_that("%/% and %% follow base R, including zero divisors", {
  vals <- c(0L, 1L, -1L, 2L, -2L, 3L, 7L, -7L, 100L, -100L, imax, -imax, 46341L, NA)
  grid <- expand.grid(x = vals, y = vals)
  for (op in c("idiv", "mod")) {
    want <- if (op == "idiv") grid$x %/% grid$y else grid$x %% grid$y
    res <- expect_tiers_identical(function() .debug_arith(op, grid$x, grid$y), label = op)
    expect_identical(res$none, want, info = op)
    for_each_tier(function(tier) expect_no_warning(.debug_arith(op, grid$x, grid$y)))
    want3 <- if (op == "idiv") grid$x %/% 3L else grid$x %% 3L
    expect_identical(.debug_arith(op, grid$x, 3L), want3)
  }
  set.seed(9)
  x <- as.integer(stats::runif(500, -imax, imax))
  y <- as.integer(stats::runif(500, -1000, 1000))
  expect_identical(.debug_arith("idiv", x, y), x %/% y)
  expect_identical(.debug_arith("mod", x, y), x %% y)
})

test_that("integer sums accumulate exactly and overflow into double", {
  res <- expect_tiers_identical(function() {
    list(
      fold(rep(imax, 3L)), fold(1:10), fold(c(-imax, -imax)), fold(seq_len(1e5)),
      fold(c(imax, 1L, -2L))
    )
  })
  expect_identical(res$none[[1]], 6442450941)
  expect_identical(res$none[[1]], sum(rep(imax, 3L)))
  expect_identical(res$none[[2]], 55L)
  expect_identical(res$none[[3]], -2 * imax)
  expect_identical(res$none[[4]], sum(seq_len(1e5)))
  expect_identical(res$none[[5]], imax - 1L)
  # A compact sequence goes through the region path.
  expect_identical(.debug_regions(1:1e6)$path, "regions")
  expect_identical(fold(1:1e6), sum(1:1e6))
})
