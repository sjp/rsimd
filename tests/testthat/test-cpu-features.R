feature_names <- c(
  "sse2", "sse3", "ssse3", "sse4_1", "sse4_2", "avx", "avx2", "fma",
  "avx512f", "avx512bw", "avx512dq", "avx512vl", "os_avx", "os_avx512",
  "neon", "fp16", "dotprod", "i8mm", "bf16", "sve", "sve2"
)
x86_names <- feature_names[1:14]
arm_names <- feature_names[15:21]
tier_ids <- c("none", "sse2", "avx2", "avx512", "neon", "sve", "sve2", "rvv", "wasm128")

# Runs `f` in a fresh R process with RSIMD_CPU_FEATURES_MASK set, since the
# mask is only read when the package is loaded.
with_cpu_mask <- function(mask, f) {
  skip_if_no_subprocess()
  # RSIMD_IMPL is cleared so that a tier pinned for this session (as in the
  # per-tier CI runs) cannot add its own warning when the mask removes it.
  callr::r(f, env = c(callr::rcmd_safe_env(), RSIMD_CPU_FEATURES_MASK = mask, RSIMD_IMPL = ""))
}

# The unmasked features as a fresh R process sees them. Under an emulator the
# subprocess may run natively and see other features than this session, so
# masked results are compared with this, not with simd_cpu_features().
child_features <- function() with_cpu_mask("", function() rsimd::simd_cpu_features())

implies <- function(f, from, to) {
  if (f[[from]]) expect_true(f[[to]], label = paste(from, "->", to))
}

test_that("simd_cpu_features() has the documented structure", {
  cf <- simd_cpu_features()
  expect_type(cf, "list")
  expect_named(cf, c("arch", "os", "method", "features", "sve_vector_length_bits", "masked"))
  expect_true(cf$arch %in% c("x86_64", "i686", "aarch64", "armv7", "other"))
  expect_true(cf$os %in% c("linux", "darwin", "windows", "freebsd", "other"))
  expect_true(cf$method %in% c("cpuid", "getauxval", "sysctl", "win32", "elf_aux_info", "none"))
  expect_type(cf$features, "logical")
  expect_named(cf$features, feature_names)
  expect_false(anyNA(cf$features))
  expect_type(cf$sve_vector_length_bits, "integer")
  expect_length(cf$sve_vector_length_bits, 1L)
  expect_type(cf$masked, "character")
})

test_that("feature implications hold", {
  cf <- simd_cpu_features()
  skip_if(cf$arch == "other", "unknown architecture")
  f <- cf$features
  chain <- c("avx512f", "avx2", "avx", "sse4_2", "sse4_1", "ssse3", "sse3", "sse2")
  for (i in seq_len(length(chain) - 1L)) implies(f, chain[i], chain[i + 1L])
  implies(f, "avx2", "os_avx")
  implies(f, "fma", "avx")
  implies(f, "os_avx512", "os_avx")
  for (bit in c("avx512bw", "avx512dq", "avx512vl")) implies(f, bit, "avx512f")
  implies(f, "avx512f", "os_avx512")
  implies(f, "sve2", "sve")
  for (bit in c("fp16", "dotprod", "i8mm", "bf16", "sve")) implies(f, bit, "neon")
  if (f[["sve"]]) {
    expect_true(cf$sve_vector_length_bits %in% c(0L, 128L * 1:16))
  } else {
    expect_identical(cf$sve_vector_length_bits, 0L)
  }
})

test_that("features are consistent with the architecture", {
  cf <- simd_cpu_features()
  f <- cf$features
  if (cf$arch == "aarch64") {
    expect_true(f[["neon"]])
    expect_false(any(f[x86_names]))
  }
  if (cf$arch == "armv7") {
    expect_false(any(f[x86_names]))
    expect_false(any(f[c("fp16", "sve", "sve2")]))
  }
  if (cf$arch %in% c("x86_64", "i686")) {
    expect_false(any(f[arm_names]))
    expect_identical(cf$method, "cpuid")
  }
  if (cf$arch == "x86_64") {
    expect_true(f[["sse2"]])
  }
  if (cf$arch == "other") {
    expect_false(any(f))
    expect_identical(cf$method, "none")
  }
  if (cf$os == "darwin" && cf$arch == "aarch64") {
    expect_false(f[["sve"]])
    expect_identical(cf$method, "sysctl")
  }
  expect_identical(cf$masked, character(0))
})

test_that("CPU tier support follows the feature bits", {
  tiers <- simd_cpu_tiers()
  f <- simd_cpu_features()$features
  expect_type(tiers, "logical")
  expect_named(tiers, tier_ids)
  expect_true(tiers[["none"]])
  expect_false(tiers[["rvv"]])
  expect_false(tiers[["wasm128"]])
  expect_identical(tiers[["sse2"]], f[["sse2"]])
  expect_identical(tiers[["avx2"]], f[["avx2"]] && f[["fma"]] && f[["os_avx"]])
  expect_identical(
    tiers[["avx512"]],
    all(f[c("avx512f", "avx512bw", "avx512dq", "avx512vl", "os_avx512")])
  )
  expect_identical(tiers[["neon"]], f[["neon"]])
  expect_identical(tiers[["sve"]], f[["sve"]])
  expect_identical(tiers[["sve2"]], f[["sve2"]])
})

test_that("RSIMD_CPU_FEATURES_MASK switches features off, with their dependents", {
  before <- child_features()
  res <- with_cpu_mask("AVX2 , neon", function() {
    list(cf = rsimd::simd_cpu_features(), tiers = rsimd:::simd_cpu_tiers())
  })
  f <- res$cf$features
  off <- c("avx2", "avx512f", "avx512bw", "avx512dq", "avx512vl", arm_names)
  expect_false(any(f[off]))
  expect_false(any(res$tiers[c("avx2", "avx512", "neon", "sve", "sve2")]))
  expect_identical(res$cf$sve_vector_length_bits, 0L)

  # Only features that were on are reported as masked; everything else is unchanged.
  expect_setequal(res$cf$masked, intersect(off, names(before$features)[before$features]))
  keep <- setdiff(feature_names, off)
  expect_identical(f[keep], before$features[keep])
})

test_that("RSIMD_CPU_FEATURES_MASK is a no-op for unsupported features", {
  before <- child_features()
  absent <- names(before$features)[!before$features]
  skip_if(length(absent) == 0L, "every feature is supported")
  res <- with_cpu_mask(paste(absent, collapse = ","), function() rsimd::simd_cpu_features())
  expect_identical(res$features, before$features)
  expect_identical(res$masked, character(0))
})

test_that("RSIMD_CPU_FEATURES_MASK warns about unknown names", {
  msgs <- with_cpu_mask("avx-512, sse2", function() {
    w <- character(0)
    withCallingHandlers(
      loadNamespace("rsimd"),
      warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }
    )
    w
  })
  expect_length(msgs, 1L)
  expect_match(msgs, "avx-512", fixed = TRUE)
  expect_no_match(msgs, "sse2", fixed = TRUE)
})
