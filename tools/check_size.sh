#!/bin/sh
# Report the size of vendored code and of the source tarball.
#
# Prints byte totals for src/vendor/simde, src/vendor/sleef (if present) and
# the source tarball. Fails if the tarball exceeds 5 MB (CRAN's guidance) or
# the SIMDe subset exceeds its budget.
#
# Usage: sh tools/check_size.sh [rsimd_<version>.tar.gz]
# Without an argument the package is built with R CMD build in a temporary
# directory and that tarball is measured.

set -eu

SIMDE_BUDGET_BYTES=3500000
TARBALL_LIMIT_BYTES=5000000

root=$(cd "$(dirname "$0")/.." && pwd -P)

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

dir_bytes() {
  if [ -d "$1" ]; then
    find "$1" -type f -exec cat {} + | wc -c | tr -d ' '
  else
    echo 0
  fi
}

status=0

simde_bytes=$(dir_bytes "$root/src/vendor/simde")
echo "src/vendor/simde: $simde_bytes bytes (budget $SIMDE_BUDGET_BYTES)"
if [ "$simde_bytes" -gt "$SIMDE_BUDGET_BYTES" ]; then
  echo "  FAIL: SIMDe subset exceeds its budget"
  status=1
fi

if [ -d "$root/src/vendor/sleef" ]; then
  echo "src/vendor/sleef: $(dir_bytes "$root/src/vendor/sleef") bytes"
else
  echo "src/vendor/sleef: not present"
fi

if [ $# -ge 1 ]; then
  tarball=$1
else
  echo "Building the source tarball with R CMD build"
  (cd "$tmp" && "${R_HOME:+$R_HOME/bin/}R" CMD build --no-build-vignettes --no-manual \
    "$root" > build.log 2>&1) || {
    cat "$tmp/build.log" >&2
    echo "R CMD build failed" >&2
    exit 1
  }
  tarball=$(ls "$tmp"/*.tar.gz)
fi
[ -f "$tarball" ] || { echo "no such tarball: $tarball" >&2; exit 1; }

tar_bytes=$(wc -c < "$tarball" | tr -d ' ')
echo "$(basename "$tarball"): $tar_bytes bytes (limit $TARBALL_LIMIT_BYTES)"
if [ "$tar_bytes" -gt "$TARBALL_LIMIT_BYTES" ]; then
  echo "  FAIL: source tarball exceeds 5 MB"
  status=1
fi

exit $status
