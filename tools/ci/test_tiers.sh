#!/bin/sh
# Runs the whole test suite once per available tier, each time starting the
# session with RSIMD_IMPL set to that tier, after checking that the tiers
# every machine of this architecture must have were built and detected.
#
# Usage: sh tools/ci/test_tiers.sh [--expect "TIER ..."] [--out FILE] [--logs DIR]
#   --expect  required tiers besides "none" (default: "sse2 avx2" on x86-64,
#             "neon" on arm64)
#   --out     where to write the simd_available() list (default: tiers.txt)
#   --logs    write each run's output to DIR/test-<tier>.log
# Every tier is run even after a failure; the exit status is non-zero if
# any run failed or a required tier is missing.

set -eu

root=$(cd "$(dirname "$0")/../.." && pwd -P)
case "$(uname -m)" in
  x86_64 | amd64 | AMD64) expect="sse2 avx2" ;;
  aarch64 | arm64) expect="neon" ;;
  *) expect= ;;
esac
out=tiers.txt logs=''
while [ $# -gt 0 ]; do
  case "$1" in
    --expect) expect=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --logs) logs=$2; shift 2 ;;
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
# shellcheck disable=SC2013 # one tier name per line, no spaces
for t in $(cat "$out"); do
  sh "$root/tools/ci/run_tests.sh" --impl "$t" --expect-impl "$t" \
    ${logs:+--log "$logs/test-$t.log"} || status=1
done
exit $status
