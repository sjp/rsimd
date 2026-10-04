#!/bin/sh
# Checks the Makevars templates, the tier flags and any generated Makevars
# for flags and assignments CRAN does not accept:
#   - no optimisation level (-O...), -ffast-math, -funsafe-math-optimizations,
#     -march=native, -mtune=native, -mcpu=native, -w, -Werror or -Wno-*;
#   - no assignment to make variables other than PKG_CPPFLAGS, PKG_CFLAGS,
#     PKG_LIBS, OBJECTS and the package's own RSIMD_* variables (so no
#     CC=, CFLAGS= or CPPFLAGS= overrides);
#   - no instruction-set flag (-m..., -march=...) in PKG_CPPFLAGS or
#     PKG_CFLAGS: those apply to every file, and tier flags belong to the
#     per-tier rules only;
#   - every recipe line compiles with $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS).
# Comments are ignored.
#
# Usage: sh tools/check_makevars.sh [generated Makevars ...]
# Always checks src/Makevars.in, src/Makevars.win.in and the flags column of
# tools/tiers.txt; also checks src/Makevars and src/Makevars.win when they
# exist, and any file given as an argument (run ./configure first to check
# the generated file). Prints each offending line as file:line: text; exit
# status 1 if any.

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
cd "$root"

files="src/Makevars.in src/Makevars.win.in"
for f in src/Makevars src/Makevars.win; do
  [ -f "$f" ] && files="$files $f"
done
[ $# -eq 0 ] || files="$files $*"

forbidden='(^|[[:space:]])(-O[^[:space:]]*|-ffast-math|-funsafe-math-optimizations|-march=native|-mtune=native|-mcpu=native|-w|-Werror[^[:space:]]*|-Wno-[^[:space:]]*)([[:space:]]|$)'

status=0

report() {
  echo "$1"
  status=1
}

# Prints file:line:text for the lines of a file with comments removed.
strip() {
  awk -v f="$1" '{ sub(/#.*/, ""); if ($0 ~ /[^[:space:]]/) print f ":" NR ":" $0 }' "$1"
}

for f in $files; do
  [ -f "$f" ] || { report "$f: no such file"; continue; }
  stripped=$(strip "$f")

  hits=$(printf '%s\n' "$stripped" | grep -E "$forbidden" || true)
  [ -z "$hits" ] || report "$(printf '%s\n' "$hits" | sed 's/$/  [forbidden flag]/')"

  hits=$(printf '%s\n' "$stripped" |
    grep -E '^[^:]*:[0-9]+:[A-Za-z_][A-Za-z0-9_]*[[:space:]]*[:+?]?=' |
    grep -Ev '^[^:]*:[0-9]+:(PKG_CPPFLAGS|PKG_CFLAGS|PKG_LIBS|OBJECTS|RSIMD_[A-Z0-9_]*)[[:space:]]*[:+?]?=' || true)
  [ -z "$hits" ] || report "$(printf '%s\n' "$hits" | sed 's/$/  [assignment to a variable the package must not set]/')"

  hits=$(printf '%s\n' "$stripped" |
    grep -E '^[^:]*:[0-9]+:PKG_C(PP)?FLAGS[[:space:]]*[:+?]?=.*(^|[[:space:]])-m' || true)
  [ -z "$hits" ] || report "$(printf '%s\n' "$hits" | sed 's/$/  [instruction-set flag outside a tier rule]/')"

  hits=$(awk -v f="$f" '/^\t/ && $0 !~ /^\t\$\(CC\) \$\(ALL_CPPFLAGS\) \$\(ALL_CFLAGS\) / { print f ":" NR ":" $0 }' "$f")
  [ -z "$hits" ] || report "$(printf '%s\n' "$hits" | sed "s/\$/  [recipe not using \$(CC) \$(ALL_CPPFLAGS) \$(ALL_CFLAGS)]/")"
done

# Fields: tier | arch | os | flags | headers | enabled | sleef
hits=$(awk -F'|' '{ sub(/#.*/, "") } NF >= 4 { print "tools/tiers.txt:" NR ":" $4 }' tools/tiers.txt |
  grep -E "$forbidden" || true)
[ -z "$hits" ] || report "$(printf '%s\n' "$hits" | sed 's/$/  [forbidden flag]/')"

if [ "$status" -eq 0 ]; then
  echo "Makevars check passed: $files tools/tiers.txt"
fi
exit $status
