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
# configure's list and the tiers disabled at build time as attributes
# "configured" and "disabled". Internal, used by tests.
simd_compiled_tiers <- function() {
  .Call(C_simd_compiled_tiers)
}
