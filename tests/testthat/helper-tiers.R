# Runs `f` once per implementation tier. Only the scalar tier exists for now;
# this is replaced by a loop over simd_available() once dispatch exists.
for_each_tier <- function(f) f("none")
