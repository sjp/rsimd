# Helpers for the integer64 tests. They do not need bit64: integer64
# vectors are built from, and read back as, their exact bit patterns.

# The 8 little-endian bytes of each element of an integer64 (or any double)
# vector, as a raw matrix with one column per element.
i64_bytes <- function(x) {
  matrix(writeBin(as.vector(unclass(x)), raw(), size = 8, endian = "little"), nrow = 8)
}

# An integer64 vector from a raw matrix of little-endian bytes.
i64_from_bytes <- function(b) {
  n <- length(b) %/% 8
  structure(readBin(as.vector(b), "double", n = n, size = 8, endian = "little"),
    class = "integer64"
  )
}

# The 64 bits of each element, least significant first, as a 64-row 0/1
# integer matrix; and back.
i64_bits <- function(x) matrix(as.integer(rawToBits(as.vector(i64_bytes(x)))), nrow = 64)
i64_from_bits <- function(bits) i64_from_bytes(packBits(as.integer(bits), "raw"))

# An integer64 vector from the signed high and unsigned low 32-bit halves
# (doubles); missing halves give NA.
i64_halves <- function(hi, lo) {
  n <- max(length(hi), length(lo))
  hi <- rep_len(hi, n) %% 2^32
  lo <- rep_len(lo, n)
  bytes <- function(u) vapply(0:3, function(k) floor(u / 256^k) %% 256, numeric(n))
  b <- cbind(matrix(bytes(lo), n), matrix(bytes(hi), n))
  i64_from_bytes(as.raw(t(b)))
}

# The halves of each element: list(hi = signed high word, lo = unsigned
# low word), as doubles.
i64_parts <- function(x) {
  m <- matrix(as.integer(i64_bytes(x)), nrow = 8)
  w <- 256^(0:3)
  hi <- colSums(m[5:8, , drop = FALSE] * w)
  list(hi = ifelse(hi >= 2^31, hi - 2^32, hi), lo = colSums(m[1:4, , drop = FALSE] * w))
}

# integer64 from integer-valued doubles of magnitude below 2^63 (exact for
# every such double), NA giving NA.
i64 <- function(x) {
  x <- as.double(x)
  na <- is.na(x)
  x[na] <- 0
  hi <- floor(x / 2^32)
  lo <- x - hi * 2^32
  hi[na] <- -2^31
  i64_halves(hi, lo)
}

# integer64 from decimal strings, exactly ("NA" or NA gives NA).
i64_dec <- function(s) {
  one <- function(s) {
    if (is.na(s) || s == "NA") {
      return(c(-2^31, 0))
    }
    neg <- startsWith(s, "-")
    hi <- 0
    lo <- 0
    for (d in as.integer(strsplit(sub("^[-+]", "", s), "")[[1L]])) {
      lo <- lo * 10 + d
      hi <- (hi * 10 + floor(lo / 2^32)) %% 2^32
      lo <- lo %% 2^32
    }
    if (neg) {
      hi <- (2^32 - 1 - hi + (lo == 0)) %% 2^32
      lo <- (2^32 - lo) %% 2^32
    }
    c(hi, lo)
  }
  p <- vapply(as.character(s), one, numeric(2), USE.NAMES = FALSE)
  i64_halves(p[1, ], p[2, ])
}

# Exact decimal strings of an integer64 vector (NA for NA).
i64_str <- function(x) {
  p <- i64_parts(x)
  vapply(seq_along(p$hi), function(i) {
    hi <- p$hi[i] %% 2^32
    lo <- p$lo[i]
    if (hi == 2^31 && lo == 0) {
      return(NA_character_)
    }
    neg <- hi >= 2^31
    if (neg) {
      hi <- (2^32 - 1 - hi + (lo == 0)) %% 2^32
      lo <- (2^32 - lo) %% 2^32
    }
    digits <- character(0)
    repeat {
      r <- hi %% 10
      hi <- hi %/% 10
      t <- r * 2^32 + lo
      lo <- t %/% 10
      digits <- c(as.character(t %% 10), digits)
      if (hi == 0 && lo == 0) break
    }
    paste0(if (neg) "-", paste(digits, collapse = ""))
  }, character(1))
}

# The double nearest each element (hi * 2^32 + lo rounds once, to even),
# NA for NA: the reference for simd_as_double().
i64_value <- function(x) {
  p <- i64_parts(x)
  v <- p$hi * 2^32 + p$lo
  v[p$hi == -2^31 & p$lo == 0] <- NA
  v
}

# Concatenation keeping the class (c() drops it without bit64).
i64c <- function(...) {
  structure(unlist(lapply(list(...), function(v) as.vector(unclass(v)))), class = "integer64")
}

# Elements idx of an integer64 vector, keeping the class.
i64_at <- function(x, idx) structure(as.vector(unclass(x))[idx], class = "integer64")

max64 <- i64_dec("9223372036854775807")
na64 <- i64_dec("NA")

# Random integer64 values: full 64-bit patterns ("full"), values in
# (-2^31, 2^31) ("small", products exact), or a mix with the boundaries;
# a fraction na_frac of them NA.
rand_i64 <- function(n, kind = c("mix", "full", "small"), na_frac = 0.05, seed = 1L) {
  kind <- match.arg(kind)
  with_seed(seed, {
    full <- i64_halves(floor(stats::runif(n) * 2^32) - 2^31, floor(stats::runif(n) * 2^32))
    small <- i64(round(stats::runif(n, -2^31 + 1, 2^31 - 1)))
    x <- switch(kind,
      full = full,
      small = small,
      mix = {
        v <- as.vector(unclass(full))
        pick <- stats::runif(n) < 0.5
        v[pick] <- as.vector(unclass(small))[pick]
        edges <- as.vector(unclass(i64_dec(c(
          "9223372036854775807", "-9223372036854775807", "0", "1", "-1",
          "4294967296", "-4294967296", "3037000500", "-3037000499", "9007199254740993"
        ))))
        put <- which(stats::runif(n) < 0.05)
        v[put] <- edges[(seq_along(put) - 1L) %% length(edges) + 1L]
        structure(v, class = "integer64")
      }
    )
    v <- as.vector(unclass(x))
    v[stats::runif(n) < na_frac] <- as.vector(unclass(na64))
    structure(v, class = "integer64")
  })
}

# Expects got to be the integer64 vector want, bit for bit, showing both
# in decimal when they differ.
expect_i64 <- function(got, want, info = NULL) {
  ok <- inherits(got, "integer64") && identical(attributes(got), attributes(want)) &&
    identical(as.vector(unclass(got)), as.vector(unclass(want)), num.eq = FALSE)
  testthat::expect(ok, sprintf(
    "%s: got %s, want %s", paste(info, collapse = " "),
    if (inherits(got, "integer64")) paste(i64_str(got), collapse = " ") else deparse1(got),
    paste(i64_str(want), collapse = " ")
  ))
  invisible(got)
}

# Expects f(...) to give the integer64 vector `want` with every tier
# selected, raising exactly the warnings `msgs`.
expect_tiers_i64 <- function(want, f, ..., msgs = character(0)) {
  res <- with_each_tier(function() value_and_warnings(f(...)))
  for (tier in names(res)) {
    expect_i64(res[[tier]]$value, want, info = paste("tier", tier))
    expect_identical(res[[tier]]$warnings, msgs, info = paste("tier", tier))
  }
  invisible(res)
}
