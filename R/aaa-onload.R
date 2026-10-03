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
  warn_unknown_cpu_mask()
  invisible()
}

# The C side ignores names it does not know in RSIMD_CPU_FEATURES_MASK; say so
# here, because a typo would otherwise silently leave a feature switched on.
warn_unknown_cpu_mask <- function() {
  mask <- Sys.getenv("RSIMD_CPU_FEATURES_MASK", "")
  if (!nzchar(mask)) {
    return(invisible())
  }
  tokens <- tolower(trimws(strsplit(mask, ",", fixed = TRUE)[[1L]]))
  unknown <- setdiff(tokens[nzchar(tokens)], names(simd_cpu_features()$features))
  if (length(unknown)) {
    warning(
      "unknown CPU feature name(s) in RSIMD_CPU_FEATURES_MASK ignored: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  invisible()
}
