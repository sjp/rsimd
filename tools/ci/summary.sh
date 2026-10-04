#!/bin/sh
# Appends one row to the CI job summary table: label, platform, result and
# duration. In GitHub Actions the table goes to $GITHUB_STEP_SUMMARY (the
# header is written with the first row); elsewhere the row is printed.
#
# Usage: sh tools/ci/summary.sh <label> <result> <seconds>

set -eu

[ $# -eq 3 ] || { echo "usage: sh tools/ci/summary.sh <label> <result> <seconds>" >&2; exit 2; }
platform="${RUNNER_OS:-$(uname -s)} $(uname -m)"
out=${GITHUB_STEP_SUMMARY:-/dev/stdout}
if [ "$out" != /dev/stdout ] && ! grep -q '^| Run |' "$out" 2> /dev/null; then
  printf '| Run | Platform | Result | Duration |\n|-----|----------|--------|----------|\n' >> "$out"
fi
printf '| %s | %s | %s | %dm %02ds |\n' "$1" "$platform" "$2" $(($3 / 60)) $(($3 % 60)) >> "$out"
