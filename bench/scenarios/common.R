# Helpers for the statistical scenarios (bench/scenarios/*.R): the
# scenario() constructor, loading, timing, error measures and the markdown
# writer. Sourced by run.R; the scenario files do not need it (each defines
# a minimal scenario() when this file has not been sourced), so a scenario
# file runs on its own as a sample.

# A scenario: a calculation base R has no function for, as a fair base R
# function and an rsimd function with the same signature.
#   id, family, title, why  identification and the one-line reason it exists
#   flagship                TRUE for the headline scenarios
#   control                 TRUE for the scenarios kept to show where rsimd
#                           does not help
#   sizes, quick_sizes      the n swept (quick_sizes for --quick)
#   setup(n)                the data, a list, from n (called once per n,
#                           after set.seed)
#   base(d), rsimd(d)       the computation, same algebra on both sides
#   base_idiomatic(d)       the idiomatic base R form, when the fair base
#                           uses a cheaper formulation (fairness rule 1)
#   refs                    named list of implementations from CRAN (or
#                           base) packages, each list(pkg = , fun = ) and
#                           optionally part(base_value, d), the part of the
#                           base result it computes; run only when the
#                           package is installed, and their errors are shown
#                           but never checked
#   reference(d)            the value errors are measured against; default
#                           base(d)
#   error(got, want)        the error measure; default the normwise
#                           relative error max|got - want| / max|want|
#   tolerance               the largest error allowed for the base and rsimd
#                           columns; a named vector by precision mode for a
#                           scenario that sweeps them
#   precision               precision modes the rsimd columns run under
#   gap                     list(id = "G4", ideal = function(d) ...): a
#                           single pass over the data standing in for the
#                           missing primitive, timed at the largest n
scenario <- function(id, family, title, why, setup, base, rsimd, sizes,
                     quick_sizes = sizes[seq_len(min(2L, length(sizes)))],
                     base_idiomatic = NULL, refs = list(), reference = NULL,
                     error = NULL, tolerance = 1e-12, precision = "fast",
                     flagship = FALSE, control = FALSE, gap = NULL) {
  stopifnot(is.function(setup), is.function(base), is.function(rsimd))
  list(
    id = id, family = family, title = title, why = why, setup = setup, base = base,
    rsimd = rsimd, sizes = sizes, quick_sizes = quick_sizes,
    base_idiomatic = base_idiomatic, refs = refs, reference = reference, error = error,
    tolerance = tolerance, precision = precision, flagship = flagship, control = control,
    gap = gap
  )
}

# The scenarios defined by the files NN-*.R in `dir`, in file order. Each
# file's last value is its scenario (or a list of them).
load_scenarios <- function(dir) {
  files <- sort(list.files(dir, pattern = "^[0-9][0-9]-.*[.]R$", full.names = TRUE))
  out <- list()
  for (f in files) {
    v <- source(f, local = new.env(parent = globalenv()))$value
    if (!is.null(v$id)) v <- list(v)
    for (s in v) out[[s$id]] <- s
  }
  out
}

# Normwise relative error of numeric results (lists are flattened).
rel_error <- function(got, want) {
  got <- as.numeric(unlist(got))
  want <- as.numeric(unlist(want))
  if (length(got) != length(want)) {
    return(Inf)
  }
  if (!identical(is.na(got), is.na(want))) {
    return(Inf)
  }
  ok <- !is.na(want)
  if (!any(ok)) {
    return(0)
  }
  scale <- max(abs(want[ok]))
  if (scale == 0) {
    return(max(abs(got[ok])))
  }
  max(abs(got[ok] - want[ok])) / scale
}

# Timing budget per cell, as bench/run.R's.
cell_budget <- function(quick) {
  if (quick) {
    list(min_time = 0.05, min_iterations = 3L, max_iterations = 50L)
  } else {
    list(min_time = 0.15, min_iterations = 5L, max_iterations = 200L)
  }
}

# Times f(d) with bench::mark: list(median_us, mem_bytes, value).
time_cell <- function(f, d, quick) {
  b <- cell_budget(quick)
  value <- f(d)
  m <- bench::mark(f(d),
    check = FALSE, memory = TRUE, filter_gc = FALSE, min_time = b$min_time,
    min_iterations = b$min_iterations, max_iterations = b$max_iterations,
    time_unit = "us"
  )
  list(
    median_us = as.numeric(m$median), mem_bytes = as.numeric(m$mem_alloc),
    value = value
  )
}

fmt_time <- function(us) {
  ifelse(is.na(us), "",
    ifelse(us < 1e3, sprintf("%s µs", signif(us, 3)),
      ifelse(us < 1e6, sprintf("%s ms", signif(us / 1e3, 3)), sprintf("%s s", signif(us / 1e6, 3)))
    )
  )
}

fmt_mem <- function(b) {
  ifelse(is.na(b), "",
    ifelse(b < 1024, sprintf("%.0f B", b),
      ifelse(b < 1024^2, sprintf("%.0f KB", b / 1024), sprintf("%.1f MB", b / 1024^2))
    )
  )
}

fmt_n <- function(n) sub("e\\+?0*([0-9])", "e\\1", format(n, scientific = TRUE))

fmt_err <- function(e) ifelse(is.na(e), "", ifelse(e == 0, "0", sprintf("%.1e", e)))
