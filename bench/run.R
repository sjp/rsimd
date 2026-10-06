# Benchmark suite: times representative rsimd operations once per available
# implementation tier and against the base R equivalent, for several input
# sizes, and writes a machine-stamped results file (rds and csv) plus a
# markdown summary (latest.md). Informational only: nothing here is asserted
# anywhere, and results from different machines are not comparable.
#
# Usage:
#   Rscript bench/run.R [--quick] [--ops sum,add] [--sizes 1e3,1e5] [--out dir]
#
#   --quick   sizes 1e3 and 1e5 only, shorter timing budgets (a local sanity
#             run of a minute or two)
#   --ops     comma-separated subset of the operations below
#   --sizes   comma-separated input lengths (default 1e3,1e5,1e7)
#   --out     output directory (default: the results/ folder next to this
#             script)
#
# Needs rsimd (installed) and bench. See bench/README.md for how to read the
# output.

suppressPackageStartupMessages(library(rsimd))
if (!requireNamespace("bench", quietly = TRUE)) {
  stop("the bench package is required: install.packages(\"bench\")", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Command line

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(file.path(getwd(), "bench"))
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
}

parse_args <- function(args) {
  opts <- list(quick = FALSE, ops = NULL, sizes = NULL, out = NULL)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    value <- function() {
      if (grepl("=", a, fixed = TRUE)) {
        return(sub("^[^=]*=", "", a))
      }
      i <<- i + 1L
      if (i > length(args)) stop("missing value for ", a, call. = FALSE)
      args[[i]]
    }
    key <- sub("=.*$", "", a)
    if (key == "--quick") {
      opts$quick <- TRUE
    } else if (key == "--ops") {
      opts$ops <- strsplit(value(), ",", fixed = TRUE)[[1L]]
    } else if (key == "--sizes") {
      opts$sizes <- as.numeric(strsplit(value(), ",", fixed = TRUE)[[1L]])
    } else if (key == "--out") {
      opts$out <- value()
    } else if (key %in% c("-h", "--help")) {
      cat("Usage: Rscript bench/run.R [--quick] [--ops a,b] [--sizes 1e3,1e5] [--out dir]\n")
      quit(status = 0L)
    } else {
      stop("unknown argument '", a, "'", call. = FALSE)
    }
    i <- i + 1L
  }
  opts
}

opts <- parse_args(commandArgs(TRUE))

# Input sizes as written on the command line (1e5, not 1e+05).
fmt_n <- function(n) sub("e\\+?0*([0-9])", "e\\1", format(n, scientific = TRUE))

# ---------------------------------------------------------------------------
# Operations: one row per op x type. `simd` and `base` are expressions over
# the inputs x, y, z (and xi, yi for integer); `ref` is an optional extra
# reference (BLAS for dot). `bound` marks the compute-bound ops for which a
# SIMD tier barely faster than `none` is suspicious.

ops <- list(
  list(
    op = "sum", type = "double", simd = quote(simd_sum(x)), base = quote(sum(x)), bound = TRUE
  ),
  list(
    op = "sum", type = "integer", simd = quote(simd_sum(xi)), base = quote(sum(xi)), bound = TRUE
  ),
  list(
    op = "mean", type = "double", simd = quote(simd_mean(x)), base = quote(mean(x)), bound = TRUE
  ),
  list(
    op = "dot", type = "double", simd = quote(simd_dot(x, y)), base = quote(sum(x * y)),
    ref = quote(crossprod(x, y)[1L]), ref_name = "blas", bound = TRUE
  ),
  list(
    op = "add", type = "double", simd = quote(simd_add(x, y)), base = quote(x + y), bound = FALSE
  ),
  list(
    op = "add", type = "integer", simd = quote(simd_add(xi, yi)), base = quote(xi + yi),
    bound = FALSE
  ),
  list(
    op = "fma", type = "double", simd = quote(simd_fma(x, y, z)), base = quote(x * y + z),
    bound = FALSE
  ),
  list(
    op = "pmax", type = "double", simd = quote(simd_pmax(x, y)), base = quote(pmax(x, y)),
    bound = TRUE
  ),
  list(
    op = "exp", type = "double", simd = quote(simd_exp(xe)), base = quote(exp(xe)), bound = TRUE
  ),
  list(
    op = "any_na", type = "double", simd = quote(simd_any_na(x)), base = quote(anyNA(x)),
    bound = TRUE
  ),
  list(
    op = "any_na", type = "integer", simd = quote(simd_any_na(xi)), base = quote(anyNA(xi)),
    bound = TRUE
  ),
  list(
    op = "is_na", type = "double", simd = quote(simd_is_na(x)), base = quote(is.na(x)),
    bound = FALSE
  ),
  list(
    op = "which", type = "logical", simd = quote(simd_which(xl)), base = quote(which(xl)),
    bound = FALSE
  ),
  list(
    op = "count", type = "logical", simd = quote(simd_count(xl)), base = quote(sum(xl)),
    bound = TRUE
  ),
  list(
    op = "max_abs", type = "double", simd = quote(simd_max_abs(x)), base = quote(max(abs(x))),
    bound = TRUE
  ),
  list(
    op = "which_max_abs", type = "double", simd = quote(simd_which_max_abs(x)),
    base = quote(which.max(abs(x))), bound = TRUE
  ),
  list(
    op = "as_integer", type = "double", simd = quote(simd_as_integer(xc, mode = "truncating")),
    base = quote(as.integer(xc)), bound = FALSE
  ),
  list(
    op = "hamming", type = "double", simd = quote(simd_hamming(x, xh)), base = quote(sum(x != xh)),
    bound = TRUE
  ),
  list(
    op = "hamming", type = "integer", simd = quote(simd_hamming(xi, xih)),
    base = quote(sum(xi != xih)), bound = TRUE
  ),
  list(
    op = "is_whole", type = "double", simd = quote(simd_is_whole(xw)),
    base = quote(xw == trunc(xw)), bound = FALSE
  ),
  list(
    op = "is_pow2", type = "double", simd = quote(simd_is_pow2(xw)),
    base = quote(xw > 0 & log2(abs(xw)) == trunc(log2(abs(xw)))), bound = FALSE
  ),
  list(
    op = "recip_approx", type = "double", simd = quote(simd_recip_approx(xp)),
    base = quote(1 / xp), bound = FALSE
  ),
  list(
    op = "rsqrt", type = "double", simd = quote(simd_rsqrt(xp)), base = quote(1 / sqrt(xp)),
    bound = FALSE
  ),
  list(
    op = "rsqrt_approx", type = "double", simd = quote(simd_rsqrt_approx(xp)),
    base = quote(1 / sqrt(xp)), bound = TRUE
  ),
  list(
    op = "rootn", type = "double", simd = quote(simd_rootn(xp, 3L)), base = quote(xp^(1 / 3)),
    bound = TRUE
  ),
  list(
    op = "mul", type = "complex", simd = quote(simd_mul(cx, cy)), base = quote(cx * cy),
    bound = FALSE
  ),
  list(
    op = "div", type = "complex", simd = quote(simd_div(cx, cy)), base = quote(cx / cy),
    bound = TRUE
  ),
  list(
    op = "prod", type = "complex", simd = quote(simd_prod(cu)), base = quote(prod(cu)),
    bound = TRUE
  ),
  list(
    op = "abs", type = "complex", simd = quote(simd_abs(cx)), base = quote(Mod(cx)), bound = TRUE
  ),
  list(
    op = "pow_int", type = "complex", simd = quote(simd_pow(cs, 3)), base = quote(cs^3),
    bound = FALSE
  )
)

# Elementary functions timed in both accuracy modes (simd_math_accuracy())
# at the largest size: the "math" table. Not in the main table.
math_ops <- list(
  list(op = "sin", type = "double", simd = quote(simd_sin(x)), base = quote(sin(x))),
  list(op = "log", type = "double", simd = quote(simd_log(xp)), base = quote(log(xp))),
  list(op = "tanh", type = "double", simd = quote(simd_tanh(xt)), base = quote(tanh(xt))),
  list(
    op = "atan2", type = "double", simd = quote(simd_atan2(y, x)), base = quote(atan2(y, x))
  ),
  list(
    op = "hypot", type = "double", simd = quote(simd_hypot(x, y)),
    base = quote(sqrt(x * x + y * y))
  ),
  # Two of the functions the neon tier hands to the C math library (issue
  # 036); pow has no fast variant.
  list(op = "pow", type = "double", simd = quote(simd_pow(xp, xt)), base = quote(xp^xt)),
  list(op = "asinh", type = "double", simd = quote(simd_asinh(x)), base = quote(asinh(x))),
  list(op = "sqrt", type = "complex", simd = quote(simd_sqrt(cs)), base = quote(sqrt(cs))),
  list(op = "exp", type = "complex", simd = quote(simd_exp(cs)), base = quote(exp(cs))),
  list(op = "log", type = "complex", simd = quote(simd_log(cs)), base = quote(log(cs))),
  list(op = "sin", type = "complex", simd = quote(simd_sin(cs)), base = quote(sin(cs))),
  list(op = "asin", type = "complex", simd = quote(simd_asin(cs)), base = quote(asin(cs))),
  # On base R's branch cut: real parts beyond +-1, imaginary parts zero.
  list(
    op = "asin_cut", type = "complex", simd = quote(simd_asin(ccut)),
    base = quote(asin(ccut))
  ),
  list(
    op = "pow", type = "complex", simd = quote(simd_pow(cs, 0.5 + 0.5i)),
    base = quote(cs^(0.5 + 0.5i))
  )
)
math_modes <- c("accurate", "fast")

# Per-call overhead: small inputs, timed on the auto tier only against base
# R, with a bare .Call() of the sum entry point as the floor (the
# "overhead" table). Not in the main table.
overhead_ops <- list(
  list(
    op = "sum", type = "double", simd = quote(simd_sum(x)), base = quote(sum(x)),
    ref = quote(.Call(rsimd:::C_simd_sum, x, FALSE, NULL, NULL)), ref_name = "bare"
  ),
  list(op = "add", type = "double", simd = quote(simd_add(x, y)), base = quote(x + y)),
  list(op = "exp", type = "double", simd = quote(simd_exp(xe)), base = quote(exp(xe))),
  list(op = "eq", type = "double", simd = quote(simd_eq(x, y)), base = quote(x == y)),
  list(
    op = "as_integer", type = "double", simd = quote(simd_as_integer(xc)),
    base = quote(as.integer(xc))
  ),
  list(op = "dot", type = "double", simd = quote(simd_dot(x, y)), base = quote(sum(x * y))),
  list(op = "mul", type = "complex", simd = quote(simd_mul(cx, cy)), base = quote(cx * cy))
)
overhead_sizes <- c(1, 10, 100)
op_names <- unique(vapply(c(ops, math_ops, overhead_ops), `[[`, "", "op"))

if (!is.null(opts$ops)) {
  unknown <- setdiff(opts$ops, op_names)
  if (length(unknown) > 0L) {
    stop("unknown op(s): ", paste(unknown, collapse = ", "),
      "; available: ", paste(op_names, collapse = ", "),
      call. = FALSE
    )
  }
  ops <- Filter(function(o) o$op %in% opts$ops, ops)
  math_ops <- Filter(function(o) o$op %in% opts$ops, math_ops)
  overhead_ops <- Filter(function(o) o$op %in% opts$ops, overhead_ops)
}

sizes <- opts$sizes
if (is.null(sizes)) sizes <- if (opts$quick) c(1e3, 1e5) else c(1e3, 1e5, 1e7)
if (anyNA(sizes) || any(sizes < 1)) stop("--sizes must be positive numbers", call. = FALSE)

# Tables beyond the main one: 1% NA inputs for sum and any_na, and the
# precision modes for sum and dot at the largest size.
na_ops <- c("sum", "any_na")
precision_ops <- c("sum", "dot")
precision_modes <- c("fast", "pairwise", "compensated")

tiers <- simd_available()
compute_bound_flag <- 1.1 # a tier below this multiple of `none` gets a warning sign

# ---------------------------------------------------------------------------
# Inputs: generated once per size from a fixed seed.

make_inputs <- function(n) {
  set.seed(20261003)
  env <- new.env(parent = globalenv())
  env$x <- stats::runif(n, -100, 100)
  env$y <- stats::runif(n, -100, 100)
  env$z <- stats::runif(n, -100, 100)
  env$xe <- stats::runif(n, -50, 50)
  env$xc <- stats::runif(n, -1e6, 1e6)
  env$xi <- sample.int(2e6L, n, replace = TRUE) - 1000000L
  env$yi <- sample.int(2e6L, n, replace = TRUE) - 1000000L
  # Drawn last so that the inputs above stay the same as in older runs.
  env$xp <- 10^stats::runif(n, -3, 3)
  env$xt <- stats::runif(n, -5, 5)
  # Half the pairs equal (Hamming); half whole numbers, some powers of two.
  half <- stats::runif(n) < 0.5
  env$xh <- ifelse(half, env$x, env$y)
  env$xih <- ifelse(half, env$xi, env$yi)
  env$xw <- ifelse(stats::runif(n) < 0.5, round(env$x), env$x)
  p2 <- stats::runif(n) < 0.1
  env$xw[p2] <- 2^sample(-20:20, sum(p2), replace = TRUE)
  # Complex inputs from the doubles above (no further draws); cu has
  # modulus 1, so its product stays finite.
  env$cx <- complex(real = env$x, imaginary = env$y)
  env$cy <- complex(real = env$z, imaginary = env$xt)
  env$cu <- complex(modulus = 1, argument = env$x)
  # Parts within +-5 for the elementary functions (no overflow).
  env$cs <- complex(real = env$xt, imaginary = env$y / 20)
  # A logical vector, TRUE where x > 0 (no further draws).
  env$xl <- env$x > 0
  # |x| in [1.01, 10] on the real axis (base R's branch cut of asin).
  env$ccut <- as.complex(sign(env$x) * (1.01 + abs(env$x) * 0.0899))
  env
}

sprinkle_na <- function(env) {
  n <- length(env$x)
  out <- new.env(parent = globalenv())
  for (nm in ls(env)) assign(nm, get(nm, env), envir = out)
  k <- max(1L, as.integer(round(n / 100)))
  out$x[sample.int(n, k)] <- NA_real_
  out$xi[sample.int(n, k)] <- NA_integer_
  out
}

# ---------------------------------------------------------------------------
# Timing

budget <- function(n) {
  min_time <- if (n >= 1e7) 0.5 else if (n >= 1e5) 0.2 else 0.1
  if (opts$quick) {
    list(min_time = min_time / 2, min_iterations = 10L, max_iterations = 100L)
  } else {
    list(min_time = min_time, min_iterations = 20L, max_iterations = 200L)
  }
}

# Times `expr` (evaluated in `env`) under implementation `impl` (NULL for
# base R code). bench::mark quotes its arguments, so the call is built.
time_expr <- function(expr, env, n, impl = NULL) {
  b <- budget(n)
  call <- as.call(c(
    list(quote(bench::mark), expr),
    list(
      check = FALSE, min_time = b$min_time, min_iterations = b$min_iterations,
      max_iterations = b$max_iterations, filter_gc = TRUE, memory = FALSE,
      time_unit = "us"
    )
  ))
  run <- function() suppressWarnings(eval(call, env))
  if (is.null(impl)) run() else simd_with_impl(impl, run())
}

acc <- new.env()
acc$marks <- list()
acc$rows <- list()

record <- function(table, spec, mode, n, impl, m) {
  key <- paste(table, spec$op, spec$type, mode, fmt_n(n), impl, sep = "/")
  acc$marks[[key]] <- m
  acc$rows[[length(acc$rows) + 1L]] <- data.frame(
    table = table, op = spec$op, type = spec$type, mode = mode, n = n, impl = impl,
    median_us = as.numeric(m$median), itr_sec = as.numeric(m$`itr/sec`),
    n_itr = as.integer(m$n_itr), stringsAsFactors = FALSE
  )
}

# Times spec on every tier and base R with `mode` set by `setter`
# (simd_precision() or simd_math_accuracy()).
run_spec <- function(table, spec, env, n, mode = "fast", setter = simd_precision) {
  old <- setter(mode)
  on.exit(setter(old))
  for (tier in tiers) {
    record(table, spec, mode, n, tier, time_expr(spec$simd, env, n, tier))
  }
  record(table, spec, mode, n, "base", time_expr(spec$base, env, n))
  if (!is.null(spec$ref)) {
    record(table, spec, mode, n, spec$ref_name, time_expr(spec$ref, env, n))
  }
}

started <- Sys.time()
for (n in sizes) {
  env <- make_inputs(n)
  message(sprintf("n = %s: main table", fmt_n(n)))
  for (spec in ops) run_spec("main", spec, env, n)

  na_specs <- Filter(function(o) o$op %in% na_ops, ops)
  if (length(na_specs) > 0L) {
    message(sprintf("n = %s: 1%% NA table", fmt_n(n)))
    na_env <- sprinkle_na(env)
    for (spec in na_specs) run_spec("na", spec, na_env, n)
  }

  prec_specs <- Filter(function(o) o$op %in% precision_ops && o$type == "double", ops)
  if (n == max(sizes) && length(prec_specs) > 0L) {
    message(sprintf("n = %s: precision-mode table", fmt_n(n)))
    for (spec in prec_specs) {
      for (mode in precision_modes) run_spec("precision", spec, env, n, mode)
    }
  }

  if (n == max(sizes) && length(math_ops) > 0L) {
    message(sprintf("n = %s: math accuracy table", fmt_n(n)))
    for (spec in math_ops) {
      for (mode in math_modes) run_spec("math", spec, env, n, mode, simd_math_accuracy)
    }
  }
}

if (length(overhead_ops) > 0L) {
  message("overhead table")
  for (n in overhead_sizes) {
    env <- make_inputs(n)
    for (spec in overhead_ops) {
      record("overhead", spec, "fast", n, "auto", time_expr(spec$simd, env, n, "auto"))
      record("overhead", spec, "fast", n, "base", time_expr(spec$base, env, n))
      if (!is.null(spec$ref)) {
        record("overhead", spec, "fast", n, spec$ref_name, time_expr(spec$ref, env, n))
      }
    }
  }
}
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

# ---------------------------------------------------------------------------
# Speedups and the accidentally-scalar flag

res <- do.call(rbind, acc$rows)
group <- paste(res$table, res$op, res$type, res$mode, res$n)
ref_time <- function(impl) {
  t <- res$median_us[res$impl == impl]
  t[match(group, group[res$impl == impl])]
}
res$speedup_none <- ref_time("none") / res$median_us
res$speedup_base <- ref_time("base") / res$median_us
bound_ops <- unique(vapply(Filter(function(o) isTRUE(o$bound), ops), `[[`, "", "op"))
# Main table only: with NAs present the early-exit scans stop after a few
# hundred elements, so every tier takes the same time.
res$flag <- res$table == "main" & res$impl %in% setdiff(tiers, "none") & res$op %in% bound_ops &
  res$n >= 1e5 & !is.na(res$speedup_none) & res$speedup_none < compute_bound_flag

# ---------------------------------------------------------------------------
# Metadata

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

git_sha <- function() {
  sha <- run_cmd("git", c("-C", shQuote(script_dir()), "rev-parse", "HEAD"))
  if (length(sha) == 0L) {
    return(NA_character_)
  }
  dirty <- run_cmd(
    "git", c("-C", shQuote(script_dir()), "status", "--porcelain", "--untracked-files=no")
  )
  if (length(dirty) > 0L) paste0(sha[1L], "-dirty") else sha[1L]
}

r_config <- function(var) {
  out <- run_cmd(file.path(R.home("bin"), "R"), c("CMD", "config", var))
  if (length(out) == 0L) NA_character_ else paste(out, collapse = " ")
}

cpu <- simd_cpu_features()
sysinfo <- Sys.info()[c("sysname", "release", "machine")]
features <- names(cpu$features)[cpu$features]
if (cpu$sve_vector_length_bits > 0L) {
  features <- c(features, paste0("sve", cpu$sve_vector_length_bits))
}
model <- cpu_model()
sha <- git_sha()
timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

metadata <- list(
  timestamp = timestamp,
  rsimd_version = as.character(utils::packageVersion("rsimd")),
  bench_version = as.character(utils::packageVersion("bench")),
  git_sha = sha,
  r_version = R.version.string,
  cc = r_config("CC"),
  cflags = r_config("CFLAGS"),
  sysinfo = sysinfo,
  cpu_model = model,
  cores = parallel::detectCores(),
  cpu_features = cpu,
  available = tiers,
  auto = simd_with_impl("auto", as.character(simd_current())),
  machine_key = paste(
    sysinfo[["sysname"]], sysinfo[["machine"]], model,
    paste(features, collapse = ","),
    sep = " | "
  ),
  sizes = sizes,
  quick = opts$quick,
  seed = 20261003L,
  elapsed_sec = elapsed
)

# ---------------------------------------------------------------------------
# Writers

out_dir <- if (is.null(opts$out)) file.path(script_dir(), "results") else opts$out
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
slug <- tolower(gsub("[^A-Za-z0-9]+", "-", paste(sysinfo[["sysname"]], sysinfo[["machine"]])))
stem <- sprintf(
  "%s-%s-%s", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), slug,
  if (is.na(sha)) "nogit" else substr(sha, 1L, 7L)
)

saveRDS(
  list(results = res, marks = acc$marks, metadata = metadata),
  file.path(out_dir, paste0(stem, ".rds"))
)

csv <- res
csv$rsimd_version <- metadata$rsimd_version
csv$timestamp <- timestamp
csv$machine <- slug
utils::write.csv(csv, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)

fmt_time <- function(us) {
  ifelse(us < 1e3, sprintf("%s µs", signif(us, 3)),
    ifelse(us < 1e6, sprintf("%s ms", signif(us / 1e3, 3)), sprintf("%s s", signif(us / 1e6, 3)))
  )
}
fmt_x <- function(r) {
  if (r < 0.01) "<0.01×" else sprintf(if (r < 0.95) "%.2f×" else "%.1f×", r)
}

cell <- function(r) {
  t <- fmt_time(r$median_us)
  if (r$impl == "base") {
    return(t)
  }
  parts <- character()
  if (r$impl %in% tiers && r$impl != "none") parts <- c(parts, paste(fmt_x(r$speedup_none), "n"))
  parts <- c(parts, paste(fmt_x(r$speedup_base), "b"))
  out <- sprintf("%s (%s)", t, paste(parts, collapse = ", "))
  if (isTRUE(r$flag)) out <- paste(out, "⚠")
  out
}

md_table <- function(sub, row_key) {
  impls <- c(tiers, "base", setdiff(unique(sub$impl), c(tiers, "base")))
  impls <- impls[impls %in% sub$impl]
  keys <- unique(sub[[row_key]])
  head <- c(
    paste0("| ", row_key, " | ", paste(impls, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(impls) + 1L), collapse = "|"), "|")
  )
  body <- vapply(keys, function(k) {
    cells <- vapply(impls, function(im) {
      r <- sub[sub[[row_key]] == k & sub$impl == im, , drop = FALSE]
      if (nrow(r) == 0L) "" else cell(r[1L, ])
    }, "")
    label <- if (row_key == "n") fmt_n(k) else k
    paste0("| ", label, " | ", paste(cells, collapse = " | "), " |")
  }, "")
  c(head, body)
}

# Keyed "op type", "overhead:op type" for the overhead table, whose
# expressions differ.
base_names <- vapply(c(ops, math_ops, overhead_ops), function(o) {
  ref <- if (is.null(o$ref)) "" else sprintf("; %s: `%s`", o$ref_name, deparse(o$ref))
  sprintf("`%s` vs base `%s`%s", deparse(o$simd), deparse(o$base), ref)
}, "")
names(base_names) <- c(
  vapply(c(ops, math_ops), function(o) paste(o$op, o$type), ""),
  vapply(overhead_ops, function(o) paste0("overhead:", o$op, " ", o$type), "")
)

md_section <- function(table, title, row_key) {
  sub_all <- res[res$table == table, , drop = FALSE]
  if (nrow(sub_all) == 0L) {
    return(character())
  }
  out <- c(paste("##", title), "")
  for (g in unique(paste(sub_all$op, sub_all$type))) {
    sub <- sub_all[paste(sub_all$op, sub_all$type) == g, , drop = FALSE]
    if (row_key == "mode") sub$mode <- as.character(sub$mode)
    out <- c(
      out, sprintf("### %s (%s)", sub$op[1L], sub$type[1L]), "",
      base_names[[if (table == "overhead") paste0("overhead:", g) else g]], "",
      md_table(sub, row_key), ""
    )
  }
  out
}

meta_lines <- c(
  "# rsimd benchmark results", "",
  "| | |", "|---|---|",
  sprintf("| Timestamp (UTC) | %s |", timestamp),
  sprintf("| rsimd | %s (git %s) |", metadata$rsimd_version, if (is.na(sha)) "unknown" else sha),
  sprintf("| R | %s |", R.version.string),
  sprintf("| Compiler | `%s` |", metadata$cc),
  sprintf("| CFLAGS | `%s` |", metadata$cflags),
  sprintf("| OS | %s %s (%s) |", sysinfo[["sysname"]], sysinfo[["release"]], sysinfo[["machine"]]),
  sprintf("| CPU | %s, %s cores |", if (is.na(model)) "unknown" else model, metadata$cores),
  sprintf("| CPU features | %s |", paste(features, collapse = ", ")),
  sprintf("| Tiers | %s (auto: %s) |", paste(tiers, collapse = ", "), metadata$auto),
  sprintf("| bench | %s%s |", metadata$bench_version, if (opts$quick) ", quick run" else ""),
  sprintf("| Run time | %.0f s |", elapsed), "",
  paste0(
    "Cells show the median time per call, then in parentheses the speedup versus the ",
    "`none` tier (`n`) and versus base R (`b`); above 1× is faster. ",
    sprintf(
      "⚠ marks a SIMD tier less than %.1f× faster than `none` on a compute-bound op at n ≥ 1e5, ",
      compute_bound_flag
    ),
    "which suggests it is running scalar code. Timings are informational: they vary between ",
    "machines and runs, and at n = 1e3 the fixed per-call overhead dominates."
  ),
  ""
)

md <- c(
  meta_lines,
  md_section("main", "Main table (precision \"fast\", no NAs)", "n"),
  md_section("na", "Inputs with 1% NA (precision \"fast\")", "n"),
  md_section("precision", sprintf("Precision modes (n = %s)", fmt_n(max(sizes))), "mode"),
  md_section("math", sprintf("Math accuracy modes (n = %s)", fmt_n(max(sizes))), "mode"),
  md_section("overhead", "Per-call overhead (auto tier)", "n")
)
writeLines(md, file.path(out_dir, paste0(stem, ".md")), useBytes = TRUE)
invisible(file.copy(
  file.path(out_dir, paste0(stem, ".md")), file.path(out_dir, "latest.md"),
  overwrite = TRUE
))

message(sprintf(
  "Wrote %s.{rds,csv,md} and latest.md to %s (%.0f s)",
  stem, normalizePath(out_dir), elapsed
))
if (any(res$flag)) {
  f <- res[res$flag, ]
  message(
    "Warning: possibly scalar tiers: ",
    paste(unique(paste(f$impl, f$op, f$type, fmt_n(f$n))), collapse = "; ")
  )
}
