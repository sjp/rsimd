#!/bin/sh
# Behavioural check of the vector layer (src/kernels/common.inc.h), the
# vector forms of src/na.h and the reduction, scan, elementwise, predicate,
# comparison, bitwise and conversion kernels (src/kernels/reduce.inc.c,
# scan.inc.c, arith.inc.c, predicates.inc.c, compare.inc.c, bitwise.inc.c,
# convert.inc.c) on every tier, against plain C references
# (tools/vector_layer_test.c), and, on the tiers that have a SLEEF header
# (unless RSIMD_DISABLE_SLEEF=1), the SLEEF elementary-function wrappers
# against long double libm and the elementary-function kernels
# (src/kernels/math.inc.c) against the none tier's scalar forms, and the
# softmax passes (src/kernels/ml.inc.c) against libm and the scalar folds.
#
# Each tier is compiled with its flags from tools/tiers.txt plus
# -ffp-contract=off and -Wall -Wextra -pedantic -Werror (R's headers, which
# src/kernel_types.h needs for R_xlen_t, are included as system headers),
# linked statically, and run:
#   - natively when the tier's arch is the host's;
#   - under qemu-x86_64 -cpu max (sse2, avx2) or qemu-i386 -cpu max (i686
#     sse2) or qemu-arm (armv7 neon) when the cross compiler and qemu-user
#     are installed;
#   - under qemu-aarch64 -cpu max with SVE vector lengths of 128, 256, 512
#     and 2048 bits (sve, sve2);
#   - avx512 is compile-only: qemu cannot emulate AVX-512 (use Intel SDE).
# Tiers without a compiler are reported as skipped.
#
# Cross compilers are found as in tools/check_simde_subset.sh
# (x86_64-linux-gnu-gcc ..., or RSIMD_CC_<arch>).
#
# Usage: sh tools/check_vector_layer.sh
# Exit status is non-zero if any tier fails to compile or run.

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
tiers="$root/tools/tiers.txt"
test_src="$root/tools/vector_layer_test.c"

r_include=$(Rscript -e 'cat(R.home("include"))' 2> /dev/null) || r_include=
if [ -z "$r_include" ] || [ ! -f "$r_include/Rinternals.h" ]; then
  echo "check_vector_layer.sh: cannot find R's include directory (is Rscript on PATH?)" >&2
  exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

host_cc=${CC:-gcc}
case "$("$host_cc" -dumpmachine 2>/dev/null || uname -m)" in
  x86_64* | amd64*) host_arch=x86_64 ;;
  i?86*) host_arch=i686 ;;
  aarch64* | arm64*) host_arch=aarch64 ;;
  arm*) host_arch=armv7 ;;
  *) host_arch=unknown ;;
esac
echo "Host compiler: $host_cc ($host_arch)"

compiler_for() {
  if [ "$1" = "$host_arch" ]; then
    echo "$host_cc"
    return
  fi
  eval "override=\${RSIMD_CC_$1:-}"
  if [ -n "$override" ]; then
    echo "$override"
    return
  fi
  case "$1" in
    x86_64) echo x86_64-linux-gnu-gcc ;;
    aarch64) echo aarch64-linux-gnu-gcc ;;
    i686) echo i686-linux-gnu-gcc ;;
    armv7) echo arm-linux-gnueabihf-gcc ;;
    *) echo none ;;
  esac
}

# Prints the commands (one per line) that run the binary for a tier/arch,
# or nothing if it cannot be run here.
runners_for() {
  tier=$1
  arch=$2
  case "$tier/$arch" in
    avx512/*) return 0 ;;
    sve/aarch64 | sve2/aarch64)
      if command -v qemu-aarch64 > /dev/null 2>&1; then
        for vl in 128 256 512 2048; do
          echo "qemu-aarch64 -cpu max,sve-default-vector-length=$((vl / 8))"
        done
      fi
      return 0
      ;;
  esac
  if [ "$arch" = "$host_arch" ]; then
    echo ""
    return 0
  fi
  case "$arch" in
    x86_64) q=qemu-x86_64 ;;
    i686) q=qemu-i386 ;;
    aarch64) q=qemu-aarch64 ;;
    armv7) q=qemu-arm ;;
    *) return 0 ;;
  esac
  if command -v "$q" > /dev/null 2>&1; then
    echo "$q -cpu max"
  fi
}

pass=0
fail=0
skip=0
done_list=" "

# Fields: tier | arch | os | flags | headers | enabled | sleef
while IFS='|' read -r tier arch os flags headers enabled sleef; do
  tier=$(echo "$tier" | tr -d ' \t')
  case "$tier" in '' | '#'*) continue ;; esac
  arch=$(echo "$arch" | tr -d ' \t')
  enabled=$(echo "$enabled" | tr -d ' \t')
  flags=$(echo "$flags" | sed 's/^[ \t]*//; s/[ \t]*$//')
  [ "$flags" = "-" ] && flags=
  [ "$enabled" = yes ] || continue
  # Tiers with a SLEEF header also test its wrappers. The header is a
  # system include here because this build adds -Wextra.
  sleef=$(echo "$sleef" | tr -d ' \t')
  if [ "$sleef" != "-" ] && [ -n "$sleef" ] && [ "${RSIMD_DISABLE_SLEEF:-0}" != 1 ]; then
    flags="$flags -DRSIMD_HAVE_SLEEF_$(echo "$tier" | tr 'a-z' 'A-Z')=1 -isystem $root/src/vendor/sleef"
    # As configure does for GCC.
    case "$tier/${CC:-gcc}" in sve/*clang* | sve2/*clang*) ;; sve/* | sve2/*) flags="$flags -fno-tree-vrp" ;; esac
  fi

  # The none tier runs once per available architecture.
  if [ "$arch" = all ]; then
    arches="$host_arch x86_64"
  else
    arches=$arch
  fi
  for a in $arches; do
    label="$tier ($a)"
    # Only the first row that compiles counts for a tier/arch pair.
    case "$done_list" in *" $tier/$a "*) continue ;; esac

    cc=$(compiler_for "$a")
    if [ "$cc" = none ] || ! command -v "$cc" > /dev/null 2>&1; then
      echo "SKIP $label: no compiler for $a (tried '$cc'; set RSIMD_CC_$a)"
      skip=$((skip + 1))
      continue
    fi

    bin="$tmp/test_${tier}_$a"
    # shellcheck disable=SC2086 # $flags is a word list
    if ! "$cc" -std=gnu11 -O2 -Wall -Wextra -pedantic -Werror -ffp-contract=off $flags \
      -DRSIMD_TIER="$tier" -I"$root/src" -I"$root/src/vendor/simde" -isystem "$r_include" \
      "$test_src" -o "$bin" -static -lm > "$tmp/log" 2>&1; then
      echo "FAIL $label: does not compile with $cc${flags:+ $flags}"
      sed 's/^/    /' "$tmp/log" | head -40
      fail=$((fail + 1))
      continue
    fi
    done_list="$done_list$tier/$a "

    runners=$(runners_for "$tier" "$a")
    if [ -z "$runners" ] && [ "$a" != "$host_arch" ] || [ "$tier" = avx512 ]; then
      echo "PASS $label: compiled with $cc${flags:+ $flags} (not run here)"
      pass=$((pass + 1))
      continue
    fi
    ok=yes
    echo "$runners" | while IFS= read -r run; do
      # shellcheck disable=SC2086 # $run is a command prefix
      if out=$($run "$bin" 2>&1); then
        echo "     ${run:-native}: $out" | head -1
      else
        echo "     ${run:-native}: FAILED"
        echo "$out" | sed 's/^/    /' | head -40
        exit 1
      fi
    done || ok=no
    if [ "$ok" = yes ]; then
      echo "PASS $label: $cc${flags:+ $flags}"
      pass=$((pass + 1))
    else
      echo "FAIL $label: $cc${flags:+ $flags}"
      fail=$((fail + 1))
    fi
  done
done < "$tiers"

echo "Vector layer check: $pass passed, $fail failed, $skip skipped"
[ "$fail" -eq 0 ]
