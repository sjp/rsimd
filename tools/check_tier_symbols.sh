#!/bin/sh
# Checks that every global symbol defined in a tier object (src/tier_<tier>.o)
# ends in _<tier>. Tier code is compiled with instruction-set flags and must
# only be reached through the dispatch tables after a CPU check; the suffix
# rule keeps tiers from clashing with each other or with the rest of the
# package, and makes any accidental direct call visible.
#
# Usage: sh tools/check_tier_symbols.sh [object directory]   (default: src)
# Run it after building the package in place (R CMD INSTALL --no-clean-on-error
# leaves the objects in src/, or use pkgbuild/devtools). NM overrides the nm
# program. Exit status is non-zero on a violation or if no tier object exists.

set -eu

LC_ALL=C
export LC_ALL

dir=${1:-$(cd "$(dirname "$0")/.." && pwd -P)/src}
nm=${NM:-nm}
found=0
bad=0

for obj in "$dir"/tier_*.o; do
  [ -f "$obj" ] || continue
  found=$((found + 1))
  tier=$(basename "$obj" .o)
  tier=${tier#tier_}
  # Defined global symbols: "address TYPE name" with an upper-case type other
  # than U. Mach-O prefixes C symbols with an underscore.
  syms=$("$nm" -g "$obj" | awk 'NF == 3 && $2 ~ /^[A-TV-Z]$/ { print $3 }' | sed 's/^_\(rsimd\)/\1/')
  n=0
  before=$bad
  for s in $syms; do
    n=$((n + 1))
    case "$s" in
      *_"$tier") ;;
      *)
        echo "FAIL tier_$tier.o: global symbol '$s' does not end in _$tier"
        bad=$((bad + 1))
        ;;
    esac
  done
  if [ "$bad" -eq "$before" ]; then
    echo "ok   tier_$tier.o: $n global symbols"
  fi
done

if [ "$found" -eq 0 ]; then
  echo "no tier objects found in $dir (build the package in place first)"
  exit 1
fi
[ "$bad" -eq 0 ]
