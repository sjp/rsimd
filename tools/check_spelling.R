# Spell-checks the Rd files, vignettes, README.md and NEWS.md with the
# spelling package against inst/WORDLIST.
#
# The benchmark snapshots in vignettes/benchmark-results are generated
# (commit hashes, compiler flags and CPU feature names), so words found only
# there are not reported.
#
# Usage: Rscript tools/check_spelling.R (from the package root)
# Prints each misspelt word with where it was found; exit status 1 if any.

snapshots <- basename(list.files("vignettes/benchmark-results", full.names = TRUE))
res <- spelling::spell_check_package(".")
in_snapshot <- vapply(res$found, function(where) {
  all(sub(":[0-9,]+$", "", where) %in% snapshots)
}, logical(1))
res <- res[!in_snapshot, , drop = FALSE]
if (nrow(res) > 0L) {
  print(res)
  cat("Add correctly spelt words to inst/WORDLIST.\n")
  quit(status = 1L)
}
cat("No spelling errors.\n")
