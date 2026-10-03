#!/bin/sh
# Vendor the pruned SIMDe subset into src/vendor/simde.
#
# Maintainer tool; it is not run at install time and needs network access.
# It fetches SIMDe at the pinned commit, computes the include closure of
# tools/simde_probe.c and replaces src/vendor/simde with exactly those headers,
# plus COPYING and a VERSION file. The hand-written README-rsimd.md in that
# directory is kept. It also refreshes the SIMDe section of inst/COPYRIGHTS.
#
# The closure is computed with `gcc -M` and with a textual #include scan; when
# both are available they must agree. Without gcc only the textual scan runs.
#
# To update SIMDe, change SIMDE_COMMIT deliberately, re-run this script, then
# run tools/check_simde_subset.sh and tools/check_size.sh.
#
# Usage: sh tools/vendor_simde.sh
# Environment: CC (default gcc), used only for the `-M` closure.

set -eu

SIMDE_REPO=https://github.com/simd-everywhere/simde-no-tests
SIMDE_COMMIT=a3fe1602e0dbad7de46f2563d79ffb132201db3b
SIMDE_BUDGET_BYTES=3500000

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
probe="$root/tools/simde_probe.c"
dest="$root/src/vendor/simde"
copyrights="$root/inst/COPYRIGHTS"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

die() {
  echo "vendor_simde.sh: $*" >&2
  exit 1
}

# Normalise relative paths read from stdin ("x86/../hedley.h" -> "hedley.h").
normpath_awk='
function normpath(p,    n, parts, out, k, i) {
  n = split(p, parts, "/")
  k = 0
  for (i = 1; i <= n; i++) {
    if (parts[i] == "" || parts[i] == ".") continue
    if (parts[i] == ".." && k > 0) { k--; continue }
    out[++k] = parts[i]
  }
  p = ""
  for (i = 1; i <= k; i++) p = p (i > 1 ? "/" : "") out[i]
  return p
}'

# Textual closure: follow every #include "..." regardless of preprocessor
# conditionals, so the result does not depend on the host architecture.
closure_textual() {
  awk -v src="$1" "$normpath_awk"'
    function dir(p) { return (p ~ /\//) ? substr(p, 1, match(p, /\/[^\/]*$/) - 1) : "" }
    function scan(file, from,    line, inc, path) {
      while ((getline line < file) > 0) {
        if (line !~ /^[ \t]*#[ \t]*include[ \t]*"/) continue
        inc = line
        sub(/^[^"]*"/, "", inc)
        sub(/".*$/, "", inc)
        path = normpath((from == "" ? "" : from "/") inc)
        if (!(path in seen)) {
          if ((getline line < (src "/" path)) < 0) {
            print "unresolved include \"" inc "\" in " file > "/dev/stderr"
            exit 1
          }
          close(src "/" path)
          seen[path] = 1
          queue[++nq] = path
        }
      }
      close(file)
    }
    BEGIN {
      nq = 0
      scan(ARGV[1], "")
      for (i = 1; i <= nq; i++) scan(src "/" queue[i], dir(queue[i]))
      for (p in seen) print p
      exit 0
    }
  ' "$2" | sort
}

# Compiler closure: the headers gcc actually opens for the probe on this host.
closure_gcc() {
  "${CC:-gcc}" -M -MT probe -std=c11 -I"$1" "$2" |
    tr ' \\' '\n\n' |
    awk -v src="$1/" "$normpath_awk"'
      index($0, src) == 1 { print normpath(substr($0, length(src) + 1)) }
    ' | sort -u
}

echo "Fetching SIMDe $SIMDE_COMMIT"
src="$tmp/simde"
git init -q "$src"
git -C "$src" fetch -q --depth 1 "$SIMDE_REPO" "$SIMDE_COMMIT"
git -C "$src" checkout -q FETCH_HEAD
[ "$(git -C "$src" rev-parse HEAD)" = "$SIMDE_COMMIT" ] || die "fetched the wrong commit"
commit_date=$(git -C "$src" log -1 --format=%cs)

version_of() {
  sed -n "s/^#define SIMDE_VERSION_$1 \\([0-9][0-9]*\\).*/\\1/p" "$src/simde-common.h"
}
simde_version="$(version_of MAJOR).$(version_of MINOR).$(version_of MICRO)"
case "$simde_version" in
  *[0-9].[0-9]*.[0-9]*) ;;
  *) die "could not read SIMDE_VERSION_* from simde-common.h" ;;
esac

echo "Computing the include closure of tools/simde_probe.c"
closure_textual "$src" "$probe" > "$tmp/closure.txt"
if command -v "${CC:-gcc}" > /dev/null 2>&1; then
  closure_gcc "$src" "$probe" > "$tmp/closure-gcc.txt"
  if ! cmp -s "$tmp/closure.txt" "$tmp/closure-gcc.txt"; then
    diff "$tmp/closure.txt" "$tmp/closure-gcc.txt" >&2 || true
    die "gcc -M and the textual include scan disagree (< textual, > gcc)"
  fi
else
  echo "  ${CC:-gcc} not found; using the textual include scan only"
fi

nfiles=$(wc -l < "$tmp/closure.txt" | tr -d ' ')
nbytes=$(cd "$src" && cat $(cat "$tmp/closure.txt") | wc -c | tr -d ' ')
echo "SIMDe subset: $nfiles headers, $nbytes bytes (budget $SIMDE_BUDGET_BYTES)"
[ "$nbytes" -le "$SIMDE_BUDGET_BYTES" ] || die "SIMDe subset exceeds its size budget"

# Stage the new tree, keeping the hand-written README.
stage="$tmp/stage"
mkdir -p "$stage"
while IFS= read -r f; do
  mkdir -p "$stage/$(dirname "$f")"
  cp "$src/$f" "$stage/$f"
done < "$tmp/closure.txt"
cp "$src/COPYING" "$stage/COPYING"
if [ -f "$dest/README-rsimd.md" ]; then
  cp "$dest/README-rsimd.md" "$stage/README-rsimd.md"
fi
cat > "$stage/VERSION" <<EOF
SIMDe subset vendored by tools/vendor_simde.sh
repository: $SIMDE_REPO
commit: $SIMDE_COMMIT
commit date: $commit_date
SIMDE_VERSION: $simde_version
headers: $nfiles
header bytes: $nbytes
EOF
find "$stage" -type f -exec chmod 644 {} +

rm -rf "$dest"
mkdir -p "$(dirname "$dest")"
cp -R "$stage" "$dest"

# Refresh the SIMDe section of inst/COPYRIGHTS.
holders="$tmp/holders.txt"
while IFS= read -r f; do
  awk '
    /^ \* Copyright:/ { on = 1; next }
    on && /^ \*[ \t]+[0-9][0-9][0-9][0-9]/ {
      line = $0
      sub(/^ \*[ \t]+[0-9][0-9-]*[ \t]+/, "", line)
      sub(/[ \t]*<.*$/, "", line)
      print line
      next
    }
    on { exit }
  ' "$src/$f"
done < "$tmp/closure.txt" | sort -u > "$holders"
cc0=$(while IFS= read -r f; do
  if grep -q 'SPDX-License-Identifier: CC0-1.0' "$src/$f"; then echo "$f"; fi
done < "$tmp/closure.txt")

section="$tmp/section.txt"
{
  echo "-- BEGIN SIMDe (generated by tools/vendor_simde.sh) --"
  echo
  echo "SIMDe (SIMD Everywhere)"
  echo "-----------------------"
  echo
  echo "Files:   src/vendor/simde (a subset of the headers, in the source package only)"
  echo "Source:  $SIMDE_REPO"
  echo "Commit:  $SIMDE_COMMIT ($commit_date)"
  echo "Version: $simde_version"
  echo "Licence: MIT, except the CC0-1.0 files listed below"
  echo
  echo "Copyright (c) 2017 Evan Nemerson."
  echo
  echo "Individual files also name these copyright holders:"
  echo
  sed 's/^/  /' "$holders"
  echo
  echo "These files are dedicated to the public domain under CC0-1.0"
  echo "(https://creativecommons.org/publicdomain/zero/1.0/):"
  echo
  for f in $cc0; do echo "  $f"; done
  echo
  echo "SIMDe licence text (src/vendor/simde/COPYING):"
  echo
  sed -e 's/^/  /' -e 's/[ \t]*$//' "$src/COPYING"
  echo
  echo "-- END SIMDe --"
} > "$section"

awk -v section="$section" '
  function emit(   line) {
    while ((getline line < section) > 0) print line
    close(section)
    done = 1
  }
  /^-- BEGIN SIMDe / { skip = 1; emit(); next }
  /^-- END SIMDe --$/ { skip = 0; next }
  skip { next }
  /^No third-party code is bundled yet\.$/ { next }
  { print }
  END { if (!done) { print ""; emit() } }
' "$copyrights" | awk 'NF { blank = 0; print; next } !blank++' > "$tmp/COPYRIGHTS"
cp "$tmp/COPYRIGHTS" "$copyrights"

echo "Vendored SIMDe $simde_version ($SIMDE_COMMIT) into src/vendor/simde"
