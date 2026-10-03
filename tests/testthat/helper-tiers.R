# Runs `f(tier)` once per available implementation tier, with that tier
# selected, and returns the results in a list named by tier.
for_each_tier <- function(f) {
  tiers <- simd_available()
  stats::setNames(lapply(tiers, function(tier) simd_with_impl(tier, f(tier))), tiers)
}
