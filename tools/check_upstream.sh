#!/bin/sh
# Report whether the vendored SIMDe and SLEEF are behind upstream releases.
#
# Maintainer tool; needs network access. Reads the vendored versions from
# src/vendor/{simde,sleef}/VERSION and compares them with the newest release
# tag (pre-release tags such as v0.8.4-rc1 are ignored) of each upstream
# repository. Prints one line per library:
#
#   <name> <vendored version> <newest release> <repository> <current|behind>
#
# The upstream-check.yml workflow runs this weekly and opens an issue for
# each library that is behind. To update, follow tools/vendor_simde.sh or
# tools/vendor_sleef.sh.
#
# Usage: sh tools/check_upstream.sh

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)

die() {
  echo "check_upstream.sh: $*" >&2
  exit 1
}

# field <file> <key>: the value of a "key: value" line.
field() {
  sed -n "s/^$2: *//p" "$1" | head -n 1
}

# newest_release <repository>: the highest x.y[.z] tag, without a leading v.
newest_release() {
  git ls-remote --tags --refs "$1" |
    sed -n 's|.*refs/tags/v\{0,1\}\([0-9][0-9]*\(\.[0-9][0-9]*\)\{1,2\}\)$|\1|p' |
    sort -t . -k 1,1n -k 2,2n -k 3,3n |
    tail -n 1
}

# newer <a> <b>: whether version a is strictly newer than version b.
newer() {
  [ "$1" != "$2" ] &&
    [ "$(printf '%s\n%s\n' "$1" "$2" | sort -t . -k 1,1n -k 2,2n -k 3,3n | tail -n 1)" = "$1" ]
}

check() {
  name=$1
  file="$root/src/vendor/$2/VERSION"
  key=$3
  [ -f "$file" ] || die "missing $file"
  repo=$(field "$file" repository)
  vendored=$(field "$file" "$key")
  [ -n "$repo" ] || die "no repository line in $file"
  [ -n "$vendored" ] || die "no $key line in $file"
  latest=$(newest_release "$repo")
  [ -n "$latest" ] || die "no release tags found in $repo"
  if newer "$latest" "$vendored"; then
    status=behind
  else
    status=current
  fi
  echo "$name $vendored $latest $repo $status"
}

check SIMDe simde SIMDE_VERSION
check SLEEF sleef tag
