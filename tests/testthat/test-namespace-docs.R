# Consistency guard for the hand-written NAMESPACE and Rd files. Uses only
# the installed package, so it also runs under R CMD check.

# Classes whose S3 methods must be registered with S3method() in NAMESPACE.
s3_classes <- "simd_vec"

# Rd topics that document no function (options, concepts).
doc_topics <- "rsimd_options"

pkg_dir <- function() dirname(system.file(package = "rsimd"))

rd_aliases <- function(rd) {
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
  vapply(rd[tags == "\\alias"], function(x) as.character(x[[1L]]), character(1))
}

# R's Rd processing drops \docType{package}, so recognise package-level
# documentation by its "<pkg>-package" name.
rd_is_package_doc <- function(rd) {
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
  name <- as.character(rd[tags == "\\name"][[1L]][[1L]])
  grepl("-package$", name)
}

test_that("every export is a function defined in the namespace", {
  ns <- asNamespace("rsimd")
  exports <- getNamespaceExports("rsimd")
  is_fun <- vapply(exports, function(n) {
    exists(n, envir = ns, inherits = FALSE) && is.function(get(n, envir = ns))
  }, logical(1))
  expect_true(all(is_fun), info = paste("not functions:", toString(exports[!is_fun])))
})

test_that("every export has an Rd alias", {
  aliases <- unlist(lapply(tools::Rd_db("rsimd"), rd_aliases), use.names = FALSE)
  undocumented <- setdiff(getNamespaceExports("rsimd"), aliases)
  expect_identical(undocumented, character(0),
    info = paste("exports without \\alias:", toString(undocumented))
  )
})

test_that("every function-like Rd alias is exported, an S3 method or a topic", {
  nsinfo <- parseNamespaceFile("rsimd", pkg_dir())
  s3 <- nsinfo$S3methods
  s3_names <- if (length(s3)) paste(s3[, 1L], s3[, 2L], sep = ".") else character(0)
  exports <- getNamespaceExports("rsimd")

  db <- tools::Rd_db("rsimd")
  aliases <- unlist(lapply(db[!vapply(db, rd_is_package_doc, logical(1))], rd_aliases),
    use.names = FALSE
  )
  function_like <- aliases[grepl("^[A-Za-z.][A-Za-z0-9._]*$", aliases)]
  stray <- setdiff(function_like, c(exports, s3_names, doc_topics))
  expect_identical(stray, character(0),
    info = paste("aliases that are not exports or S3 methods:", toString(stray))
  )
})

test_that("every S3 method in the namespace is registered", {
  nsinfo <- parseNamespaceFile("rsimd", pkg_dir())
  s3 <- nsinfo$S3methods
  registered <- if (length(s3)) paste(s3[, 1L], s3[, 2L], sep = ".") else character(0)

  pattern <- paste0("^[A-Za-z.]+\\.(", paste(s3_classes, collapse = "|"), ")$")
  ns <- asNamespace("rsimd")
  candidates <- grep(pattern, ls(ns, all.names = TRUE), value = TRUE)
  methods <- candidates[vapply(
    candidates, function(n) is.function(get(n, envir = ns)),
    logical(1)
  )]
  unregistered <- setdiff(methods, registered)
  expect_identical(unregistered, character(0),
    info = paste("S3 methods without S3method():", toString(unregistered))
  )
})
