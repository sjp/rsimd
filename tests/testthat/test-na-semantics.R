# The NA/NaN rules of reductions and elementwise operations, checked through
# the self-test kernels on every available tier.

missing_patterns <- list(
  na = NA_real_, nan = NaN, na_nan = c(NA, NaN), nan_na = c(NaN, NA),
  mixed = c(1, NA, NaN, 2)
)

# A pattern placed after `before` ones and followed by `after` twos, so the
# missing values land in every lane and in the vector body and the tail.
embed <- function(p, before, after) c(rep(1, before), p, rep(2, after))

test_that("double reductions: NA wins over NaN in any order, na.rm removes both", {
  terms <- c("x", "sq", "abs", "xy")
  for (name in names(missing_patterns)) {
    p <- missing_patterns[[name]]
    for (shift in c(0, 1, 3, 7, 17, 40)) {
      x <- embed(p, shift, (shift * 3) %% 11)
      kept <- x[!is.na(x)]
      want_rm <- list(x = sum(kept), sq = sum(kept^2), abs = sum(abs(kept)), xy = sum(kept^2))
      for (prec in prec_codes) {
        res <- expect_tiers_identical(function() {
          lapply(stats::setNames(nm = terms), function(t) {
            y <- if (t == "xy") x else NULL
            list(
              keep = .debug_fold(x, y, t, FALSE, TRUE, prec),
              rm = .debug_fold(x, y, t, TRUE, TRUE, prec)
            )
          })
        }, label = paste(name, shift, prec))
        for (t in terms) {
          keep <- res$none[[t]]$keep
          rm <- res$none[[t]]$rm
          info <- paste(name, shift, prec, t)
          if (anyNA(p) && any(is.na(p) & !is.nan(p))) {
            expect_identical(keep$value, NA_real_, info = info)
          } else {
            expect_identical(keep$value, NaN, info = info)
          }
          expect_identical(keep$count, as.double(length(x)), info = info)
          expect_identical(rm$value, want_rm[[t]], info = info)
          expect_identical(rm$count, as.double(length(kept)), info = info)
          expect_true(rm$saw_nan, info = info)
          expect_identical(rm$saw_na, any(is.na(p) & !is.nan(p)), info = info)
        }
      }
    }
  }
  # Where base R is deterministic (one kind of missing value) it agrees.
  for (p in list(NA_real_, NaN, c(1, NA, 2), c(3, NaN))) {
    expect_identical(fold(p), sum(p))
    expect_identical(fold(p, na_rm = TRUE), sum(p, na.rm = TRUE))
  }
  # Base R is order dependent here; rsimd is not.
  expect_identical(sum(c(NaN, NA)), NaN)
  expect_identical(fold(c(NaN, NA)), NA_real_)
})

test_that("a missing value in y alone makes x * y missing", {
  x <- c(1, 2, 3, 4, 5)
  expect_tiers_identical(function() {
    list(
      fold(x, c(1, NA, 1, 1, 1), "xy"), fold(x, c(1, NaN, 1, 1, 1), "xy"),
      fold(x, c(1, NA, 1, NaN, 1), "xy", na_rm = TRUE)
    )
  })
  expect_identical(fold(x, c(1, NA, 1, 1, 1), "xy"), NA_real_)
  expect_identical(fold(x, c(1, NaN, 1, 1, 1), "xy"), NaN)
  expect_identical(fold(x, c(1, NA, 1, NaN, 1), "xy", na_rm = TRUE), 9)
})

test_that("integer and logical sums: NA unless na.rm, exact 64-bit accumulation", {
  for (shift in c(0, 1, 2, 5, 9, 33)) {
    xi <- c(rep(1L, shift), NA, 2L, rep(3L, shift))
    xl <- c(rep(TRUE, shift), NA, FALSE, rep(TRUE, shift))
    expect_tiers_identical(function() {
      list(
        fold(xi), fold(xi, na_rm = TRUE), fold(xl), fold(xl, na_rm = TRUE),
        .debug_fold(xi, na_rm = TRUE)$count
      )
    }, label = shift)
    expect_identical(fold(xi), sum(xi))
    expect_identical(fold(xi, na_rm = TRUE), sum(xi, na.rm = TRUE))
    expect_identical(fold(xl), sum(xl))
    expect_identical(fold(xl, na_rm = TRUE), sum(xl, na.rm = TRUE))
  }
  expect_identical(fold(NA_integer_, na_rm = TRUE), 0L)
  expect_identical(fold(integer(0)), 0L)
  expect_identical(fold(numeric(0)), 0)
})

test_that("na_check = FALSE skips detection but stays memory safe", {
  x <- c(1L, NA, 3L)
  for_each_tier(function(tier) {
    v <- fold(x, na_check = FALSE)
    # NA counts as INT_MIN: an unspecified but well-defined number.
    expect_true(is.numeric(v) && length(v) == 1L, info = tier)
    expect_true(is.na(fold(c(1, NA, 3), na_check = FALSE)), info = tier)
  })
  # na.rm = TRUE detects missing values whatever na_check says.
  expect_tiers_identical(function() {
    list(fold(x, na_rm = TRUE, na_check = FALSE), fold(c(1, NA, 3), na_rm = TRUE, na_check = FALSE))
  })
  expect_identical(fold(x, na_rm = TRUE, na_check = FALSE), 4L)
})

test_that("finished reductions apply the NA, NaN and empty-input rules", {
  fin <- function(op, type, fields = list(), na_rm = FALSE, n = 10) {
    .debug_finish(op, type, n, fields, 0L, na_rm)$value
  }
  real_ops <- c(
    "sum", "prod", "mean", "min", "max", "sum_sq", "sum_abs", "dot", "norm",
    "dist", "cosine", "var", "sd"
  )
  for (op in real_ops) {
    full <- list(count = 10, f64 = 2)
    expect_identical(fin(op, "double", c(full, saw_na = 1, saw_nan = 1)), NA_real_, info = op)
    # var and sd give NA for NaN too, as base R.
    want_nan <- if (op %in% c("var", "sd")) NA_real_ else NaN
    expect_identical(fin(op, "double", c(full, saw_nan = 1)), want_nan, info = op)
    # na.rm = TRUE: the kernel already removed them.
    expect_identical(fin(op, "double", c(full, saw_na = 1, saw_nan = 1), na_rm = TRUE), 2,
      info = op
    )
  }
  # NA takes the result type.
  expect_identical(fin("sum", "integer", list(saw_na = 1, count = 1)), NA_integer_)
  expect_identical(fin("min", "logical", list(saw_na = 1, count = 1)), NA_integer_)
  expect_identical(fin("mean", "integer", list(saw_na = 1, count = 1)), NA_real_)
  skip_if_not_installed("bit64")
  expect_identical(fin("sum", "integer64", list(saw_na = 1, count = 1)), bit64::NA_integer64_)
})

test_that("empty input after na.rm follows base R", {
  fin <- function(op, type = "double", count = 0, na_rm = TRUE) {
    .debug_finish(op, type, 3, list(count = count, saw_na = 1, saw_nan = 1), 0L, na_rm)$value
  }
  expect_identical(fin("sum"), 0)
  expect_identical(fin("prod"), 1)
  expect_identical(fin("mean"), NaN)
  expect_identical(fin("var", count = 1), NA_real_)
  expect_identical(fin("sd", count = 1), NA_real_)
  expect_identical(fin("var", count = 0), NA_real_)
  expect_identical(fin("which_min"), integer(0))
  expect_identical(fin("which_max", "integer"), integer(0))
  for (type in c("double", "integer", "logical")) {
    expect_warning(v <- fin("min", type), "no non-missing arguments to min; returning Inf",
      fixed = TRUE
    )
    expect_identical(v, Inf)
    expect_warning(v <- fin("max", type), "no non-missing arguments to max; returning -Inf",
      fixed = TRUE
    )
    expect_identical(v, -Inf)
  }
  # The same texts as base R.
  expect_warning(min(numeric(0)), "no non-missing arguments to min; returning Inf", fixed = TRUE)
  # Without na.rm, an empty input is empty too.
  expect_identical(.debug_finish("mean", "double", 0, list(), 0L)$value, NaN)
  expect_warning(.debug_finish("min", "double", 0, list(), 0L), "returning Inf", fixed = TRUE)
  # var of one element is NA, as var(1).
  expect_identical(.debug_finish("var", "double", 1, list(count = 1, f64 = 0), 0L)$value, var(1))
})

test_that("any and all use three-valued logic", {
  cases <- list(
    c(NA, TRUE), c(NA, FALSE), c(TRUE, NA), c(FALSE, NA), NA, TRUE, FALSE, logical(0),
    c(TRUE, FALSE), c(NA, NA)
  )
  for (x in cases) {
    for (shift in c(0, 3, 15, 40)) {
      # Padding with the op's neutral value does not change the answer.
      xa <- c(rep(FALSE, shift), x, rep(FALSE, shift %% 7))
      xb <- c(rep(TRUE, shift), x, rep(TRUE, shift %% 7))
      expect_tiers_identical(function() {
        list(
          .debug_lgl(xa, "any")$value, .debug_lgl(xa, "any", na_rm = TRUE)$value,
          .debug_lgl(xb, "all")$value, .debug_lgl(xb, "all", na_rm = TRUE)$value
        )
      }, label = deparse(x))
      info <- paste(deparse(x), shift)
      expect_identical(.debug_lgl(xa, "any")$value, any(xa), info = info)
      expect_identical(.debug_lgl(xa, "any", na_rm = TRUE)$value, any(xa, na.rm = TRUE),
        info = info
      )
      expect_identical(.debug_lgl(xb, "all")$value, all(xb), info = info)
      expect_identical(.debug_lgl(xb, "all", na_rm = TRUE)$value, all(xb, na.rm = TRUE),
        info = info
      )
    }
  }
  # Integers: any non-zero, non-NA value is TRUE.
  expect_identical(.debug_lgl(c(0L, NA, 5L))$value, TRUE)
  expect_identical(.debug_lgl(c(0L, NA), "all")$value, FALSE)
})

test_that("any and all may stop at the first deciding element", {
  x <- c(rep(FALSE, 100), TRUE, NA, rep(FALSE, 5000))
  for_each_tier(function(tier) {
    r <- .debug_lgl(x, "any")
    expect_identical(r$value, TRUE, info = tier)
    full <- .debug_lgl(x, "scan")
    expect_true(full$saw_na && full$any_true && full$any_false, info = tier)
  })
  # The NA after the deciding TRUE does not matter, whether or not it was seen.
  expect_identical(.debug_lgl(c(TRUE, rep(NA, 10000)))$value, TRUE)
  expect_identical(.debug_lgl(c(FALSE, rep(NA, 10000)), "all")$value, FALSE)
})

test_that("elementwise doubles keep NA and NaN apart", {
  ops <- c("add", "sub", "mul", "div", "pmin", "pmax")
  for (op in ops) {
    for_each_tier(function(tier) {
      info <- paste(op, tier)
      r <- .debug_arith(op, c(NA, NaN, 1), 1)
      expect_identical(is.na(r), c(TRUE, TRUE, FALSE), info = info)
      expect_identical(is.nan(r), c(FALSE, TRUE, FALSE), info = info)
      # NA wins over NaN in either operand position.
      r <- .debug_arith(op, c(NA, NaN, NA, NaN), c(NaN, NA, NA, 1))
      expect_identical(is.nan(r), c(FALSE, FALSE, FALSE, TRUE), info = info)
      expect_identical(is.na(r), rep(TRUE, 4), info = info)
      r <- .debug_arith(op, 2, c(1, NA, NaN))
      expect_identical(is.nan(r), c(FALSE, FALSE, TRUE), info = info)
      expect_identical(is.na(r), c(FALSE, TRUE, TRUE), info = info)
      # Without the check NA may degrade to NaN but stays missing.
      r <- .debug_arith(op, c(NA, NaN, 1), 1, na_check = FALSE)
      expect_identical(is.na(r), c(TRUE, TRUE, FALSE), info = info)
    })
  }
  expect_identical(.debug_arith("add", c(NA, NaN, 1), 1), c(NA, NaN, 2))
  expect_identical(.debug_arith("pmin", c(1, 5, NaN), c(2, 3, 0)), c(1, 3, NaN))
  expect_identical(.debug_arith("pmax", c(1, 5, 0), c(2, NaN, 0)), c(2, NaN, 0))
})

test_that("elementwise doubles match the none tier in every lane", {
  set.seed(1)
  for (n in c(0, 1, 2, 3, 5, 8, 15, 16, 17, 33, 70)) {
    x <- stats::runif(n, -5, 5)
    y <- stats::runif(n, -5, 5)
    x[seq_len(n) %% 5 == 1] <- NA
    y[seq_len(n) %% 7 == 2] <- NaN
    x[seq_len(n) %% 11 == 3] <- Inf
    for (op in c("add", "sub", "mul", "div", "pmin", "pmax")) {
      res <- expect_tiers_identical(function() .debug_arith(op, x, y), label = paste(op, n))
      # NA exactly where an operand is NA; missing wherever an operand is.
      is_na_only <- function(v) is.na(v) & !is.nan(v)
      expect_identical(is_na_only(res$none), is_na_only(x) | is_na_only(y), info = paste(op, n))
      expect_true(all(is.na(res$none)[is.na(x) | is.na(y)]), info = paste(op, n))
    }
  }
})
