#!/bin/sh
# Checks that only the R vector access layer (src/rvec.c) touches R vector
# data, lengths and interrupts directly. Kernels and .Call entry points go
# through src/rvec.h, which handles ALTREP, long vectors and interrupts.
#
# Flags, outside comments, any call of:
#   DATAPTR, DATAPTR_RO, DATAPTR_OR_NULL,
#   REAL, INTEGER, LOGICAL, RAW, COMPLEX and their _RO, _OR_NULL and
#   _GET_REGION forms, LENGTH, XLENGTH, ITERATE_BY_REGION,
#   R_CheckUserInterrupt
# in src/*.c, src/*.h and src/kernels/*. Exempt: rvec.c itself, and the
# files that only fill small result vectors for detection and control and
# never read user data:
allowed="rvec.c api_control.c cpu_features.c dispatch.c version.c init.c"
#
# Usage: sh tools/lint_c.sh
# Prints each offending line as file:line: text; exit status 1 if any.

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
cd "$root"

pattern='(^|[^A-Za-z0-9_])(DATAPTR(_RO|_OR_NULL)?|(REAL|INTEGER|LOGICAL|RAW|COMPLEX)(_RO|_OR_NULL|_GET_REGION)?|X?LENGTH|ITERATE_BY_REGION|R_CheckUserInterrupt)[[:space:]]*\('

status=0
for f in src/*.c src/*.h src/kernels/*; do
  [ -f "$f" ] || continue
  base=$(basename "$f")
  skip=0
  for a in $allowed; do
    if [ "$base" = "$a" ]; then skip=1; fi
  done
  [ "$skip" -eq 0 ] || continue

  # Blank out /* */ and // comments, keeping line numbers, then grep.
  hits=$(awk '
    {
      line = $0; out = ""
      while (length(line) > 0) {
        if (incomment) {
          e = index(line, "*/")
          if (e == 0) { line = ""; break }
          line = substr(line, e + 2); incomment = 0
        } else {
          s = index(line, "/*"); d = index(line, "//")
          if (d > 0 && (s == 0 || d < s)) { out = out substr(line, 1, d - 1); line = ""; break }
          if (s == 0) { out = out line; line = ""; break }
          out = out substr(line, 1, s - 1); line = substr(line, s + 2); incomment = 1
        }
      }
      print out
    }' "$f" | grep -nE "$pattern" || true)
  if [ -n "$hits" ]; then
    echo "$hits" | sed "s|^|$f:|"
    status=1
  fi
done

if [ "$status" -ne 0 ]; then
  echo "lint_c.sh: use the access layer (src/rvec.h) instead of the calls above" >&2
fi
exit "$status"
