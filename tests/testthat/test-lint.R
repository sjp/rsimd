# The C lint (tools/lint_c.sh) needs the package sources, which are not part
# of the built package, so this test runs only from a source checkout.

test_that("only the access layer touches R vector data directly", {
  skip_on_cran()
  script <- normalizePath(test_path("..", "..", "tools", "lint_c.sh"), mustWork = FALSE)
  skip_if_not(file.exists(script), "needs a source checkout")
  out <- suppressWarnings(system2("sh", script, stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), label = paste(out, collapse = "\n"))
})
