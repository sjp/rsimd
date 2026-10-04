#!/bin/sh
# Runs the test suite against the installed rsimd once, reports the CPU
# features and tiers seen, and adds a row to the job summary.
#
# Usage: sh tools/ci/run_tests.sh [options]
#   --impl TIER          start the session with RSIMD_IMPL=TIER (default: unset, auto)
#   --label NAME         row label in the summary (default: the implementation)
#   --runner "CMD ..."   run R under an emulator, e.g. "sde64 -spr --" or
#                        "qemu-aarch64 -cpu max,sve-default-vector-length=32"
#   --expect-impl TIER   fail unless the session's implementation is TIER
#   --expect-tiers "T.." fail unless simd_available() lists every tier named
#   --log FILE           also write the output to FILE
#
# RSIMD_TEST_SUBSET and RSIMD_TEST_TIERS are passed through to the tests (see
# tests/README.md). With --runner, R's own executable (R_HOME/bin/exec/R) is
# started directly under the emulator with the environment R's front-end
# script would set, because emulators do not reliably follow the exec chain
# from Rscript. Subprocesses started by the tests then run natively.

set -eu

root=$(cd "$(dirname "$0")/../.." && pwd -P)
impl='' label='' runner='' expect_impl='' expect_tiers='' log=''
while [ $# -gt 0 ]; do
  case "$1" in
    --impl) impl=$2; shift 2 ;;
    --label) label=$2; shift 2 ;;
    --runner) runner=$2; shift 2 ;;
    --expect-impl) expect_impl=$2; shift 2 ;;
    --expect-tiers) expect_tiers=$2; shift 2 ;;
    --log) log=$2; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[ -n "$label" ] || label=${impl:-auto}
tmplog=
[ -n "$log" ] || { log=$(mktemp); tmplog=$log; }

script=$(mktemp)
trap 'rm -f "$script" $tmplog' EXIT HUP INT TERM
cat > "$script" <<'RCODE'
library(rsimd)
cat("rsimd", format(utils::packageVersion("rsimd")), "on", R.version.string, "\n")
str(simd_cpu_features())
avail <- simd_available()
cat("simd_available():", avail, "\nsimd_current():", as.vector(simd_current()), "\n")
cat("RSIMD_TEST_SUBSET:", Sys.getenv("RSIMD_TEST_SUBSET"), " RSIMD_TEST_TIERS:", Sys.getenv("RSIMD_TEST_TIERS"), "\n")
want <- Sys.getenv("RSIMD_CI_EXPECT_IMPL")
if (nzchar(want) && !identical(as.vector(simd_current()), want)) {
  stop("expected implementation ", want, ", got ", simd_current())
}
want <- strsplit(Sys.getenv("RSIMD_CI_EXPECT_TIERS"), "[ ,]+")[[1]]
if (length(miss <- setdiff(want[nzchar(want)], avail))) {
  stop("simd_available() lacks expected tiers: ", paste(miss, collapse = ", "))
}
options(testthat.progress.max_fails = 1e6)
testthat::test_dir(Sys.getenv("RSIMD_CI_TEST_DIR"),
  load_package = "installed", package = "rsimd", stop_on_failure = TRUE
)
RCODE

# Native Windows R (from Git Bash) needs Windows paths.
if command -v cygpath > /dev/null 2>&1; then
  root=$(cygpath -m "$root")
  script=$(cygpath -m "$script")
fi

RSIMD_CI_EXPECT_IMPL=$expect_impl
RSIMD_CI_EXPECT_TIERS=$expect_tiers
RSIMD_CI_TEST_DIR=$root/tests/testthat
NOT_CRAN=true
export RSIMD_CI_EXPECT_IMPL RSIMD_CI_EXPECT_TIERS RSIMD_CI_TEST_DIR NOT_CRAN
if [ -n "$impl" ]; then RSIMD_IMPL=$impl; export RSIMD_IMPL; else unset RSIMD_IMPL; fi

if [ -n "$runner" ]; then
  # The variables R's front-end script exports, as an R session sees them.
  # shellcheck disable=SC2016
  eval "$(Rscript -e 'for (v in c("R_HOME", "R_SHARE_DIR", "R_INCLUDE_DIR", "R_DOC_DIR", "LD_LIBRARY_PATH", "R_ARCH"))
    cat(sprintf("%s=%s; export %s\n", v, shQuote(Sys.getenv(v)), v))
    cat(sprintf("R_LIBS=%s; export R_LIBS\n", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))))')"
  # shellcheck disable=SC2086 # the runner command is split into words
  set -- $runner "$R_HOME/bin${R_ARCH:-}/exec/R" --no-save --no-restore --no-echo -f "$script"
else
  set -- Rscript "$script"
fi

echo "== $label: $*"
start=$(date +%s)
status=0
{ "$@" 2>&1 || echo "rsimd-ci-exit-status $?"; } | tee "$log"
grep -q '^rsimd-ci-exit-status' "$log" && status=1
seconds=$(($(date +%s) - start))
if [ $status -eq 0 ]; then result=pass; else result=FAIL; fi
sh "$root/tools/ci/summary.sh" "$label" "$result" "$seconds"
exit $status
