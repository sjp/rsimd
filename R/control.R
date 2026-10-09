simd_version <- function() {
  .Call(C_simd_version)
}

simd_cpu_features <- function() {
  .Call(C_simd_cpu_features)
}

# Named logical over simd_tiers() ids: can this CPU and OS run the tier,
# regardless of whether it was compiled in. Internal, used by tests.
simd_cpu_tiers <- function() {
  .Call(C_simd_cpu_tiers)
}

# Tiers compiled into this build (named by the tier objects), with
# configure's list, the tiers built with SLEEF elementary functions and the
# tiers disabled at build time as attributes "configured", "sleef" and
# "disabled" (comma-separated strings). Internal, used by tests.
simd_compiled_tiers <- function() {
  .Call(C_simd_compiled_tiers)
}

# R-side record of the requested implementation ("auto" or a tier id), set
# by simd_use(); NULL until .onLoad has run.
.impl_state <- new.env(parent = emptyenv())
.impl_state$requested <- NULL

simd_tiers <- function() {
  c("none", "sse2", "avx2", "avx512", "neon", "sve", "sve2", "rvv", "wasm128")
}

simd_available <- function() {
  .Call(C_simd_available)
}

simd_current <- function() {
  .sync_impl()
  .Call(C_simd_current)
}

simd_complex_variants <- function() {
  .Call(C_simd_c128_variants)
}

simd_use <- function(impl) {
  problem <- .impl_problem(impl)
  if (!is.null(problem)) {
    .stop(problem)
  }
  previous <- .previous_impl()
  status <- .Call(C_simd_select, impl)
  if (status != 0L) {
    # Only reachable if the R and C views of availability disagree.
    .stop("could not select implementation '", impl, "'")
  }
  .impl_state$requested <- impl
  options(rsimd.impl = impl)
  invisible(previous)
}

simd_precision <- function(mode) {
  if (missing(mode)) {
    return(.precision_mode())
  }
  problem <- .precision_problem(mode)
  if (!is.null(problem)) {
    .stop(problem)
  }
  # An invalid current value is replaced, not reported, so that this call
  # repairs it; the default stands in as the previous mode.
  old <- getOption("rsimd.precision", "fast")
  if (!is.null(.precision_problem(old))) old <- "fast"
  options(rsimd.precision = mode)
  invisible(old)
}

.precision_modes <- c("fast", "pairwise", "compensated")

# For a candidate value of rsimd.precision: NULL if valid, else the error.
.precision_problem <- function(mode) {
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
    !mode %in% .precision_modes) {
    return(paste0(
      "precision mode must be one of ",
      paste0("\"", .precision_modes, "\"", collapse = ", ")
    ))
  }
  NULL
}

# The precision mode in effect: the rsimd.precision option, which the user
# may have set directly, validated.
.precision_mode <- function() {
  mode <- getOption("rsimd.precision", "fast")
  problem <- .precision_problem(mode)
  if (!is.null(problem)) {
    .stop(
      "invalid option rsimd.precision: ", problem,
      "; reset it with simd_precision()"
    )
  }
  mode
}

# For a candidate value of rsimd.na_check: NULL if valid, else the error.
.na_check_problem <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    return("must be TRUE or FALSE")
  }
  NULL
}

simd_math_accuracy <- function(mode) {
  if (missing(mode)) {
    return(.math_accuracy_mode())
  }
  problem <- .math_accuracy_problem(mode)
  if (!is.null(problem)) {
    .stop(problem)
  }
  # As in simd_precision(): replace an invalid current value.
  old <- getOption("rsimd.math_accuracy", "accurate")
  if (!is.null(.math_accuracy_problem(old))) old <- "accurate"
  options(rsimd.math_accuracy = mode)
  invisible(old)
}

.math_accuracy_modes <- c("accurate", "fast")

# For a candidate value of rsimd.math_accuracy: NULL if valid, else the
# error.
.math_accuracy_problem <- function(mode) {
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
    !mode %in% .math_accuracy_modes) {
    return(paste0(
      "math accuracy mode must be one of ",
      paste0("\"", .math_accuracy_modes, "\"", collapse = ", ")
    ))
  }
  NULL
}

# The math accuracy mode in effect: the rsimd.math_accuracy option, which
# the user may have set directly, validated.
.math_accuracy_mode <- function() {
  mode <- getOption("rsimd.math_accuracy", "accurate")
  problem <- .math_accuracy_problem(mode)
  if (!is.null(problem)) {
    .stop(
      "invalid option rsimd.math_accuracy: ", problem,
      "; reset it with simd_math_accuracy()"
    )
  }
  mode
}

simd_with_impl <- function(impl, expr) {
  old <- simd_use(impl)
  on.exit(simd_use(old))
  force(expr)
}

# For a candidate value of rsimd.impl: NULL if it can be selected, otherwise
# the error message.
.impl_problem <- function(impl) {
  if (!is.character(impl) || length(impl) != 1L || is.na(impl)) {
    return("'impl' must be a single string")
  }
  if (identical(impl, "auto")) {
    return(NULL)
  }
  if (!impl %in% simd_tiers()) {
    return(paste0(
      "unknown implementation '", impl, "'; use \"auto\" or one of: ",
      paste(simd_tiers(), collapse = ", ")
    ))
  }
  available <- simd_available()
  if (!impl %in% available) {
    return(paste0(
      "implementation '", impl, "' is not available on this machine; available: ",
      paste(available, collapse = ", ")
    ))
  }
  NULL
}

# The request simd_use() replaces: the rsimd.impl option if the user set it
# directly to a usable value since the last simd_use(), else the last
# request, so simd_with_impl() restores what was in effect.
.previous_impl <- function() {
  opt <- getOption("rsimd.impl")
  if (!is.null(opt) && is.null(.impl_problem(opt))) opt else .impl_state$requested
}

# Honours a change to the rsimd.impl option made without simd_use(). The C
# side (rsimd_entry() in src/rvec.c) does the same at the start of a call
# when the option differs from the value it last saw, so the R wrappers
# need not.
.sync_impl <- function() {
  problem <- .sync_option()
  if (!is.null(problem)) .sync_error(problem)
  invisible()
}

# Selects the rsimd.impl option if it differs from the request, or returns
# why it cannot be selected (NULL when it was). The C side evaluates this
# and raises the error itself once its sync is over, so that a calling
# handler of the error can call rsimd and sync again.
.sync_option <- function() {
  impl <- getOption("rsimd.impl", "auto")
  if (identical(impl, .impl_state$requested)) return(NULL)
  problem <- .impl_problem(impl)
  if (is.null(problem)) simd_use(impl)
  problem
}

.sync_error <- function(problem) {
  .stop("invalid option rsimd.impl: ", problem, "; reset it with simd_use()")
}

# Internal, for tests: for an available tier (default: the active one), the
# tier whose kernel each dispatch slot actually runs after fill-down.
simd_kernel_tiers <- function(tier = NULL) {
  .Call(C_simd_kernel_tiers, tier)
}

# Internal, for tests: calls the internal tier_name and fill_probe slots
# through the active table.
simd_probe_slots <- function() {
  .Call(C_simd_probe_slots)
}
