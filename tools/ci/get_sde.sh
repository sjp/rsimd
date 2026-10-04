#!/bin/sh
# Downloads Intel's Software Development Emulator (SDE) into a directory and
# prints the path of sde64. SDE's licence does not allow redistribution, so
# it is fetched from Intel at run time and never committed.
#
# Usage: sh tools/ci/get_sde.sh [DIR]   (default DIR: .sde)
# SDE_URL and SDE_SHA256 select another release (Intel's download page lists
# both); by default the pinned release below is used. An already downloaded
# tarball in DIR with the right checksum is reused.

set -eu

SDE_URL=${SDE_URL:-https://downloadmirror.intel.com/924984/sde-external-10.13.1-2026-07-28-lin.tar.xz}
SDE_SHA256=${SDE_SHA256:-94e97d623fec54385686e1e7ba65ebc9941748c05ee451423948334892bf2b50}

dir=${1:-.sde}
mkdir -p "$dir"
tarball=$dir/$(basename "$SDE_URL")
sum() { sha256sum "$1" 2> /dev/null | cut -d' ' -f1 || shasum -a 256 "$1" | cut -d' ' -f1; }
if [ ! -f "$tarball" ] || [ "$(sum "$tarball")" != "$SDE_SHA256" ]; then
  curl -fsSL --retry 3 -A "Mozilla/5.0" -o "$tarball" "$SDE_URL" >&2
fi
got=$(sum "$tarball")
if [ "$got" != "$SDE_SHA256" ]; then
  echo "SDE checksum mismatch: got $got, want $SDE_SHA256" >&2
  exit 1
fi
name=$(basename "$tarball" .tar.xz)
[ -x "$dir/$name/sde64" ] || tar -xJf "$tarball" -C "$dir"
echo "$(cd "$dir/$name" && pwd -P)/sde64"
