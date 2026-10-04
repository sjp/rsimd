#!/bin/sh
# Checks an R CMD INSTALL log and the shared library it built:
#   - every tier object (tier_<tier>.c) is compiled with -ffp-contract=off;
#   - every tier configure reported as having SLEEF is compiled with
#     -DRSIMD_HAVE_SLEEF_<TIER>=1, and no other tier is;
#   - the shared library exports no SLEEF symbol (Sleef_*): the bundled
#     SLEEF functions are static and must stay internal.
#
# Usage: sh tools/check_build.sh <install log> [shared library]
# The log is the output of R CMD INSTALL (configure's messages and the
# compiler command lines). Without a library argument, the installed
# package's libs/rsimd.so (or .dll) is checked. NM overrides the nm program.
# Exit status is non-zero on any violation.

set -eu

LC_ALL=C
export LC_ALL

[ $# -ge 1 ] || { echo "usage: sh tools/check_build.sh <install log> [shared library]" >&2; exit 2; }
log=$1
[ -f "$log" ] || { echo "no such log: $log" >&2; exit 2; }
so=${2:-}
if [ -z "$so" ]; then
  libs=$(Rscript -e 'cat(system.file("libs", package = "rsimd"))' 2> /dev/null) || libs=
  for f in "$libs"/rsimd.so "$libs"/rsimd.dll "$libs"/*/rsimd.so "$libs"/*/rsimd.dll; do
    if [ -f "$f" ]; then so=$f; break; fi
  done
fi

bad=0

sleef=$(sed -n 's/^rsimd: tiers with SLEEF elementary functions: //p' "$log" | tail -n 1)
case "$sleef" in "(no tier)") sleef= ;; esac
sleef=" $(echo "$sleef" | tr ',' ' ') "

lines=$(grep -E ' -c tier_[a-z0-9]+\.c ' "$log" || true)
if [ -z "$lines" ]; then
  echo "FAIL no tier compile lines in $log (was the package compiled?)"
  exit 1
fi
echo "$lines" | while IFS= read -r line; do
  tier=$(echo "$line" | sed 's/.* -c tier_\([a-z0-9]*\)\.c .*/\1/')
  macro="-DRSIMD_HAVE_SLEEF_$(echo "$tier" | tr 'a-z' 'A-Z')=1"
  ok=yes
  case " $line " in *" -ffp-contract=off "*) ;; *)
    echo "FAIL tier_$tier.c is compiled without -ffp-contract=off"
    ok=no
    ;;
  esac
  case "$sleef" in
    *" $tier "*) want=yes ;;
    *) want=no ;;
  esac
  case " $line " in *" $macro "*) have=yes ;; *) have=no ;; esac
  if [ "$want" != "$have" ]; then
    echo "FAIL tier_$tier.c: SLEEF reported '$want' by configure but $macro present '$have'"
    ok=no
  fi
  if [ "$ok" = yes ]; then
    echo "ok   tier_$tier.c: -ffp-contract=off, SLEEF $have"
  else
    exit 1
  fi
done || bad=1

if [ -n "$so" ] && [ -f "$so" ]; then
  nm=${NM:-nm}
  case "$(uname -s)" in
    Darwin) exported=$("$nm" -gU "$so") ;;
    *) exported=$("$nm" -D --defined-only "$so" 2> /dev/null || "$nm" -g "$so") ;;
  esac
  found=$(echo "$exported" | grep -c 'Sleef_' || true)
  if [ "$found" -ne 0 ]; then
    echo "FAIL $so exports $found SLEEF symbols:"
    echo "$exported" | grep 'Sleef_' | head -n 10
    bad=1
  else
    echo "ok   $so exports no SLEEF symbol"
  fi
else
  echo "FAIL shared library not found (pass it as the second argument)"
  bad=1
fi

[ "$bad" -eq 0 ]
