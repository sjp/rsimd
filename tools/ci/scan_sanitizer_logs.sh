#!/bin/sh
# Fails if any of the given logs contains an AddressSanitizer or
# UndefinedBehaviorSanitizer report, printing the reports found.
#
# Usage: sh tools/ci/scan_sanitizer_logs.sh FILE...

set -eu

pattern='ERROR: AddressSanitizer|runtime error:|UndefinedBehaviorSanitizer'
status=0
for f in "$@"; do
  [ -f "$f" ] || continue
  if grep -Eq "$pattern" "$f"; then
    echo "Sanitizer report in $f:"
    grep -E -B2 -A12 "$pattern" "$f" | head -n 200
    status=1
  fi
done
[ $status -eq 0 ] && echo "No sanitizer reports in $# log(s)."
exit $status
