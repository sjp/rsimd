.onLoad <- function(libname, pkgname) {
  impl <- Sys.getenv("RSIMD_IMPL", "")
  if (nzchar(impl)) {
    options(rsimd.impl = impl)
  } else if (is.null(getOption("rsimd.impl"))) {
    options(rsimd.impl = "auto")
  }
  if (is.null(getOption("rsimd.precision"))) {
    options(rsimd.precision = "fast")
  }
  if (is.null(getOption("rsimd.na_check"))) {
    options(rsimd.na_check = TRUE)
  }
  invisible()
}
