#!/bin/sh
# Runs a few small reductions on each available tier, each tier in its own
# R session under gdb, and prints the faulting instruction, registers and
# siginfo if one crashes. R's own SIGSEGV handler hides all of that (it
# prints only "address (nil), cause 'unknown'"), and gdb stops before it runs.
#
# Usage: sh tools/ci/debug_crash.sh [TIER ...]   (default: every available tier)

set -u

tmp=$(mktemp -d)
tiers=${*:-$(Rscript -e 'cat(rev(rsimd::simd_available()))')}

for tier in $tiers; do
  cat > "$tmp/repro.R" <<EOF
library(rsimd)
simd_use("$tier")
cat("tier:", attr(simd_current(), "requested"), "\n")
cat("sum:", simd_sum(c(1, 2)), "\n")
cat("popcount_total:", simd_popcount_total(c(1L, 3L)), "\n")
cat("dot:", simd_dot(1:3, c(0.5, 1, 2)), "\n")
EOF
  cat > "$tmp/cmds.gdb" <<EOF
set pagination off
set confirm off
handle SIGSEGV stop print nopass
run --vanilla --no-echo -f $tmp/repro.R
printf "---- faulting instruction\n"
x/i \$pc
printf "---- preceding instructions\n"
x/12i \$pc-48
printf "---- siginfo\n"
p \$_siginfo
printf "---- registers\n"
info registers
printf "---- backtrace\n"
bt 15
EOF
  echo "==== tier $tier"
  R -d gdb --debugger-args="-batch -x $tmp/cmds.gdb" 2>&1 | grep -v '^\[New Thread\|^\[Thread'
done
rm -rf "$tmp"
