# Machine and build metadata shared by the benchmark runners (bench/run.R
# and bench/scenarios/run.R). Needs rsimd loaded.

run_cmd <- function(cmd, args) {
  out <- tryCatch(
    suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
    return(character())
  }
  out
}

cpu_model <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Linux" && file.exists("/proc/cpuinfo")) {
    info <- readLines("/proc/cpuinfo", warn = FALSE)
    field <- function(name) {
      hit <- grep(paste0("^", name, "\\s*:"), info, value = TRUE)
      if (length(hit) == 0L) NA_character_ else trimws(sub("^[^:]*:", "", hit[1L]))
    }
    model <- field("model name")
    if (!is.na(model)) {
      return(model)
    }
    # arm64 Linux has no model name; the implementer/part codes identify the core.
    impl <- field("CPU implementer")
    part <- field("CPU part")
    if (!is.na(impl)) {
      return(sprintf("implementer %s part %s", impl, part))
    }
  } else if (sys == "Darwin") {
    model <- run_cmd("sysctl", c("-n", "machdep.cpu.brand_string"))
    if (length(model) > 0L) {
      return(model[1L])
    }
  } else if (sys == "Windows") {
    model <- Sys.getenv("PROCESSOR_IDENTIFIER")
    if (nzchar(model)) {
      return(model)
    }
  }
  NA_character_
}

git_sha <- function(dir) {
  sha <- run_cmd("git", c("-C", shQuote(dir), "rev-parse", "HEAD"))
  if (length(sha) == 0L) {
    return(NA_character_)
  }
  dirty <- run_cmd(
    "git", c("-C", shQuote(dir), "status", "--porcelain", "--untracked-files=no")
  )
  if (length(dirty) > 0L) paste0(sha[1L], "-dirty") else sha[1L]
}

r_config <- function(var) {
  out <- run_cmd(file.path(R.home("bin"), "R"), c("CMD", "config", var))
  if (length(out) == 0L) NA_character_ else paste(out, collapse = " ")
}

# The machine description both runners record: CPU features (with the SVE
# length), OS, CPU model, git SHA of the checkout at `dir`, compiler and
# flags, and the time.
machine_info <- function(dir) {
  cpu <- simd_cpu_features()
  features <- names(cpu$features)[cpu$features]
  if (cpu$sve_vector_length_bits > 0L) {
    features <- c(features, paste0("sve", cpu$sve_vector_length_bits))
  }
  list(
    cpu = cpu,
    sysinfo = Sys.info()[c("sysname", "release", "machine")],
    features = features,
    model = cpu_model(),
    sha = git_sha(dir),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    cc = r_config("CC"),
    cflags = r_config("CFLAGS"),
    cores = parallel::detectCores()
  )
}
