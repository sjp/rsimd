#!/bin/sh
# Checks that the tiers every machine of the architecture must have (sse2
# and avx2 on x86-64, neon on arm64, none everywhere) were built and
# detected, runs the whole test suite once, then runs the dispatch tests
# once per available tier, each time starting the session with RSIMD_IMPL
# set to that tier. The whole suite needs no run per tier: its kernel tests
# loop over every available tier themselves (tests/testthat/helper-tiers.R).
#
# Usage: sh tools/ci/test_tiers.sh [--expect "TIER ..."] [--out FILE] [--logs DIR] [--smoke-only]
#   --expect      required tiers besides "none" (default: "sse2 avx2" on x86-64,
#                 "neon" on arm64)
#   --out         where to write the simd_available() list (default: tiers.txt)
#   --logs        write each run's output to DIR/test-all.log and
#                 DIR/test-<tier>.log
#   --smoke-only  skip the whole-suite run, when R CMD check has run it
# Every run is made even after a failure; the exit status is non-zero if
# any run failed or a required tier is missing. RSIMD_TEST_SUBSET is passed
# through to the tests.

set -eu

root=$(cd "$(dirname "$0")/../.." && pwd -P)
case "$(uname -m)" in
  x86_64 | amd64 | AMD64) expect="sse2 avx2" ;;
  aarch64 | arm64) expect="neon" ;;
  *) expect= ;;
esac
out=tiers.txt logs='' smoke_only=false
# The test files that check the session runs the tier RSIMD_IMPL names and
# that selection and forced implementations work.
smoke='^(dispatch|onload|cpu-features|build-tiers)$'
while [ $# -gt 0 ]; do
  case "$1" in
    --expect) expect=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --logs) logs=$2; shift 2 ;;
    --smoke-only) smoke_only=true; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

Rscript -e 'cat(rsimd::simd_available(), sep = "\n")' | tr -d '\r' > "$out"
echo "Available tiers: $(tr '\n' ' ' < "$out")"
status=0
for t in none $expect; do
  if ! grep -qx "$t" "$out"; then
    echo "FAIL: required tier $t is not available (was it dropped by configure, or not detected?)"
    status=1
  fi
done

[ -z "$logs" ] || mkdir -p "$logs"
if [ "$smoke_only" = false ]; then
  sh "$root/tools/ci/run_tests.sh" --label "all tiers" \
    ${logs:+--log "$logs/test-all.log"} || status=1
fi
# shellcheck disable=SC2013 # one tier name per line, no spaces
for t in $(cat "$out"); do
  sh "$root/tools/ci/run_tests.sh" --impl "$t" --expect-impl "$t" --filter "$smoke" \
    --label "$t (dispatch)" ${logs:+--log "$logs/test-$t.log"} || status=1
done
exit $status
