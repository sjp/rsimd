#!/bin/sh
# Report the size of vendored code, the source tarball and the shared library.
#
# Prints the byte total of src/vendor/simde, the size of every file in
# src/vendor/sleef, the source tarball and, when an installed rsimd is found,
# libs/rsimd.so (as installed and stripped of debugging information), and a
# table of the installed package's top-level entries with the library
# counted stripped, as CRAN measures it. Fails if the tarball exceeds 5 MB
# (CRAN's guidance), the SIMDe subset or the SLEEF headers exceed their
# budgets, the stripped library exceeds its, the installed package exceeds
# 5 MB, or an installed directory other than libs exceeds 1 MB (libs is
# held to the library budget instead).
#
# Usage: sh tools/check_size.sh [rsimd_<version>.tar.gz]
# Without an argument the package is built with R CMD build in a temporary
# directory and that tarball is measured. RSIMD_SO names the shared library
# to measure and RSIMD_PKG_DIR the installed package (defaults: the
# installed package's, found with Rscript; RSIMD_PKG_DIR defaults to the
# package containing RSIMD_SO when that is set).

set -eu

SIMDE_BUDGET_BYTES=3500000
SLEEF_BUDGET_BYTES=2500000
SO_BUDGET_BYTES=3000000
TARBALL_LIMIT_BYTES=5000000
INSTALLED_LIMIT_BYTES=5000000
SUBDIR_LIMIT_BYTES=1000000

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
  for f in "$root"/src/vendor/sleef/*; do
    echo "  src/vendor/sleef/$(basename "$f"): $(wc -c < "$f" | tr -d ' ') bytes"
  done
  sleef_bytes=$(cat "$root"/src/vendor/sleef/sleefinline_*.h | wc -c | tr -d ' ')
  echo "src/vendor/sleef: $(dir_bytes "$root/src/vendor/sleef") bytes, headers $sleef_bytes (budget $SLEEF_BUDGET_BYTES)"
  if [ "$sleef_bytes" -gt "$SLEEF_BUDGET_BYTES" ]; then
    echo "  FAIL: SLEEF headers exceed their budget"
    status=1
  fi
else
  echo "src/vendor/sleef: not present"
fi

so=${RSIMD_SO:-}
if [ -z "$so" ]; then
  so=$("${R_HOME:+$R_HOME/bin/}Rscript" -e 'cat(system.file("libs", package = "rsimd"))' 2> /dev/null) || so=
  for f in "$so"/rsimd.so "$so"/rsimd.dll "$so"/*/rsimd.so "$so"/*/rsimd.dll; do
    if [ -f "$f" ]; then so=$f; break; fi
  done
fi
if [ -n "$so" ] && [ -f "$so" ]; then
  so_bytes=$(wc -c < "$so" | tr -d ' ')
  if command -v "${STRIP:-strip}" > /dev/null 2>&1 &&
    "${STRIP:-strip}" -S -o "$tmp/stripped" "$so" 2> /dev/null; then
    stripped=$(wc -c < "$tmp/stripped" | tr -d ' ')
  else
    stripped=$so_bytes
  fi
  echo "$so: $so_bytes bytes, $stripped without debugging information (budget $SO_BUDGET_BYTES)"
  if [ "$stripped" -gt "$SO_BUDGET_BYTES" ]; then
    echo "  FAIL: the shared library exceeds its budget"
    status=1
  fi
else
  echo "shared library: rsimd is not installed (set RSIMD_SO)"
fi

pkg=${RSIMD_PKG_DIR:-}
if [ -z "$pkg" ] && [ -n "${RSIMD_SO:-}" ]; then
  d=$(dirname "$RSIMD_SO")
  [ "$(basename "$d")" = libs ] || d=$(dirname "$d")
  pkg=$(dirname "$d")
elif [ -z "$pkg" ]; then
  pkg=$("${R_HOME:+$R_HOME/bin/}Rscript" -e 'cat(system.file(package = "rsimd"))' 2> /dev/null) || pkg=
fi
if [ -n "$pkg" ] && [ -f "$pkg/DESCRIPTION" ]; then
  echo "Installed package $pkg (libs counted without debugging information):"
  total=0
  for entry in "$pkg"/*; do
    name=$(basename "$entry")
    if [ "$name" = libs ] && [ -n "${stripped:-}" ]; then
      bytes=$(( $(dir_bytes "$entry") - so_bytes + stripped ))
    elif [ -d "$entry" ]; then
      bytes=$(dir_bytes "$entry")
    else
      bytes=$(wc -c < "$entry" | tr -d ' ')
    fi
    total=$((total + bytes))
    printf '  %-14s %10s bytes\n' "$name" "$bytes"
    if [ -d "$entry" ] && [ "$name" != libs ] && [ "$bytes" -gt "$SUBDIR_LIMIT_BYTES" ]; then
      echo "  FAIL: installed directory $name exceeds 1 MB"
      status=1
    fi
  done
  printf '  %-14s %10s bytes (limit %s)\n' total "$total" "$INSTALLED_LIMIT_BYTES"
  if [ "$total" -gt "$INSTALLED_LIMIT_BYTES" ]; then
    echo "  FAIL: the installed package exceeds 5 MB"
    status=1
  fi
else
  echo "installed package: not found (set RSIMD_PKG_DIR)"
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
