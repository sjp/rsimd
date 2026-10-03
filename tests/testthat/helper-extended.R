# Long and randomised tests run only when RSIMD_EXTENDED_TESTS is "true",
# so that the default suite stays fast.
skip_unless_extended <- function() {
  if (!identical(Sys.getenv("RSIMD_EXTENDED_TESTS"), "true")) {
    testthat::skip("extended tests not enabled (set RSIMD_EXTENDED_TESTS=true)")
  }
}
