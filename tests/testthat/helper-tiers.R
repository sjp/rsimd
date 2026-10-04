# Implementation tiers to test: every available one, best first, always
# ending with the "none" oracle. RSIMD_TEST_TIERS narrows the list (see
# helper-subset.R).
tiers_to_test <- function() {
  avail <- simd_available()
  only <- Sys.getenv("RSIMD_TEST_TIERS", "")
  if (!nzchar(only)) {
    return(avail)
  }
  only <- trimws(strsplit(only, ",", fixed = TRUE)[[1]])
  missing <- setdiff(only, avail)
  if (length(missing)) {
    stop("RSIMD_TEST_TIERS names tiers that are not available: ", paste(missing, collapse = ", "))
  }
  avail[avail %in% c(only, "none")]
}

skip_if_no_tier <- function(tier) {
  if (!tier %in% simd_available()) testthat::skip(paste("tier", tier, "unavailable"))
}

# Calls fn() once per tier with that tier selected; a list named by tier.
with_each_tier <- function(fn) {
  lapply(stats::setNames(nm = tiers_to_test()), function(t) simd_with_impl(t, fn()))
}

# Runs `f(tier)` once per available implementation tier, with that tier
# selected, and returns the results in a list named by tier.
for_each_tier <- function(f) {
  tiers <- tiers_to_test()
  stats::setNames(lapply(tiers, function(tier) simd_with_impl(tier, f(tier))), tiers)
}
