# S14: Tanimoto similarity of a query fingerprint (2048 bits, 64 integers)
# against a database of fingerprints stored as an integer matrix, one
# column per molecule. Full demo: Rscript bench/scenarios/13-tanimoto.R
#
# Sample:
#   s <- source("bench/scenarios/13-tanimoto.R")$value
#   d <- s$setup(1e4); head(s$rsimd(d))
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

words <- 64L

# Random fingerprints with about 10% of the bits set: each integer is the
# OR of three random patterns AND-ed down.
tanimoto_setup <- function(n) {
  rnd <- function(k) {
    v <- function() sample.int(65536L, k, TRUE) - 1L + (sample.int(32768L, k, TRUE) - 1L) * 65536L
    bitwAnd(bitwAnd(v(), v()), bitwOr(v(), v()))
  }
  fp <- matrix(rnd(words * n), words)
  q <- fp[, 1]
  # Bit counts of a 16-bit half, for the base R lookup.
  lut <- vapply(0:65535, function(v) sum(as.integer(intToBits(v))), 0L)
  pop <- function(m) lut[bitwAnd(m, 65535L) + 1L] + lut[bitwAnd(bitwShiftR(m, 16L), 65535L) + 1L]
  list(
    fp = fp, q = q, qrep = rep(q, n), lut = lut,
    # The bit count of every molecule, as a database index stores it.
    popdb = colSums(matrix(pop(fp), words)), popq = sum(pop(q))
  )
}

# |A and B| / (|A| + |B| - |A and B|).
tanimoto <- function(inter, d) inter / (d$popdb + d$popq - inter)

sc <- scenario(
  id = "tanimoto", family = "distance", flagship = TRUE,
  title = "Tanimoto search, 2048-bit fingerprints (n = database size)",
  why = "Base R has no population count; chemistry similarity search is all popcounts.",
  sizes = c(1e3, 1e4, 2e4),
  setup = tanimoto_setup,
  # Fair base: a 16-bit lookup table for the bit counts.
  base = function(d) {
    a <- bitwAnd(d$fp, d$qrep)
    pc <- d$lut[bitwAnd(a, 65535L) + 1L] + d$lut[bitwAnd(bitwShiftR(a, 16L), 65535L) + 1L]
    tanimoto(colSums(matrix(pc, words)), d)
  },
  base_idiomatic = function(d) {
    bits <- intToBits(bitwAnd(d$fp, d$qrep))
    tanimoto(colSums(matrix(as.integer(bits), 32L * words)), d)
  },
  rsimd = function(d) {
    tanimoto(colSums(matrix(simd_popcount(simd_bit_and(d$fp, d$qrep)), words)), d)
  },
  tolerance = 0,
  gap = list(id = "G10", ideal = function(d) simd_hamming_bits(d$fp, d$qrep))
)

# Full demo: database sizes 1e3 to 1e5.
if (sys.nframe() == 0L) {
  set.seed(20261006)
  for (n in c(1e3, 1e4, 1e5)) {
    d <- sc$setup(n)
    tb <- system.time(for (k in 1:5) sb <- sc$base(d))[["elapsed"]] / 5
    tr <- system.time(for (k in 1:5) sr <- sc$rsimd(d))[["elapsed"]] / 5
    cat(sprintf(
      "Tanimoto search, %g molecules: base %.1f ms, rsimd (%s) %.1f ms, %.1fx, identical %s\n",
      n, tb * 1e3, simd_current(), tr * 1e3, tb / tr, identical(sb, sr)
    ))
  }
}

invisible(sc)
