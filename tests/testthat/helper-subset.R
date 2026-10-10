# Reduced runs for slow environments such as emulators and valgrind.
#
# RSIMD_TEST_SUBSET selects one:
#   "quick"          input lengths are capped at 1000 (the tail and block
#                    cases around every vector width are kept);
#   "tier_emulation" as "quick", and tests that start a subprocess are
#                    skipped (QEMU user mode runs them natively, SDE
#                    slowly).
# RSIMD_TEST_TIERS (comma-separated, e.g. "avx512,none") restricts the
# tiers the kernel tests loop over; "none" is always added as the oracle.
test_subset <- function() {
  s <- Sys.getenv("RSIMD_TEST_SUBSET", "")
  if (!s %in% c("", "quick", "tier_emulation")) {
    stop("RSIMD_TEST_SUBSET must be empty, \"quick\" or \"tier_emulation\", not \"", s, "\"")
  }
  s
}

reduced_lengths <- function() test_subset() != ""

skip_if_no_subprocess <- function() {
  testthat::skip_if_not_installed("callr")
  if (identical(test_subset(), "tier_emulation")) {
    testthat::skip("subprocess tests are skipped when RSIMD_TEST_SUBSET is tier_emulation")
  }
}
