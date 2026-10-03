.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("rsimd.precision"))) {
    options(rsimd.precision = "fast")
  }
  if (is.null(getOption("rsimd.na_check"))) {
    options(rsimd.na_check = TRUE)
  }
  warn_unknown_cpu_mask()
  warn_bad_debug_stride()
  init_impl()
  invisible()
}

# Initial implementation: the rsimd.impl option if set before loading, else
# the RSIMD_IMPL environment variable, else "auto". A value that cannot be
# selected gives a warning and "auto"; loading never fails because of it.
init_impl <- function() {
  impl <- getOption("rsimd.impl")
  source <- "option rsimd.impl"
  if (is.null(impl)) {
    impl <- Sys.getenv("RSIMD_IMPL", "")
    source <- "RSIMD_IMPL"
    if (!nzchar(impl)) impl <- "auto"
  }
  problem <- .impl_problem(impl)
  if (!is.null(problem)) {
    warning(source, ": ", problem, "; using \"auto\"", call. = FALSE)
    impl <- "auto"
  }
  .impl_state$requested <- NULL
  simd_use(impl)
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

# The C side ignores an RSIMD_DEBUG_STRIDE that is not a positive integer
# and keeps the default interrupt stride; say so.
warn_bad_debug_stride <- function() {
  stride <- Sys.getenv("RSIMD_DEBUG_STRIDE", "")
  if (nzchar(stride) && !grepl("^[0-9]*[1-9][0-9]*$", stride)) {
    warning(
      "RSIMD_DEBUG_STRIDE must be a positive integer; ignoring \"", stride, "\"",
      call. = FALSE
    )
  }
  invisible()
}
