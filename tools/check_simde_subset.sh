#!/bin/sh
# Compile-check the vendored SIMDe subset for every SIMDe-based tier.
#
# For each row of tools/tiers.txt whose headers come from SIMDe, compile a
# small translation unit that includes exactly that tier's headers, with the
# tier's flags and -Wall -Wextra -pedantic -Werror. The unit also fails to
# compile unless SIMDe selected native code generation for the tier.
#
# Rows for the host architecture use $CC (default gcc). Rows for another
# architecture use a cross compiler if one is on the PATH
# (x86_64-linux-gnu-gcc, aarch64-linux-gnu-gcc, i686-linux-gnu-gcc), or the
# one named by RSIMD_CC_<arch> (e.g. RSIMD_CC_x86_64=x86_64-w64-mingw32-gcc);
# otherwise the tier is reported as skipped. Tiers that do not use SIMDe
# (none, sve, sve2, armv7 neon) are not checked here.
#
# Usage: sh tools/check_simde_subset.sh
# Exit status is non-zero if any checked tier fails to compile.

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
tiers="$root/tools/tiers.txt"
probe="$root/tools/simde_probe.c"
simde="$root/src/vendor/simde"

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
    *) echo none ;;
  esac
}

# Representative operations and the SIMDe macros that prove native codegen.
tier_body() {
  case "$1" in
    sse2) cat <<'EOF'
#if !defined(SIMDE_X86_SSE2_NATIVE)
#  error "SIMDe is not using native SSE2"
#endif
void rsimd_probe(double *out, const double *in) {
  simde__m128d a = simde_mm_loadu_pd(in);
  simde_mm_storeu_pd(out, simde_mm_add_pd(a, simde_mm_set1_pd(1.0)));
}
EOF
      ;;
    neon) cat <<'EOF'
#if !defined(SIMDE_ARM_NEON_A64V8_NATIVE)
#  error "SIMDe is not using native AArch64 NEON"
#endif
void rsimd_probe(double *out, const double *in) {
  simde__m128d a = simde_mm_loadu_pd(in);
  simde__m128d m = simde_mm_cmpunord_pd(a, a);
  simde_mm_storeu_pd(out, simde_mm_blendv_pd(simde_mm_add_pd(a, a), a, m));
}
EOF
      ;;
    avx2) cat <<'EOF'
#if !defined(SIMDE_X86_AVX2_NATIVE) || !defined(SIMDE_X86_FMA_NATIVE)
#  error "SIMDe is not using native AVX2 and FMA"
#endif
void rsimd_probe(double *out, const double *in) {
  simde__m256d a = simde_mm256_loadu_pd(in);
  simde_mm256_storeu_pd(out, simde_mm256_fmadd_pd(a, a, simde_mm256_set1_pd(1.0)));
}
EOF
      ;;
    avx512) cat <<'EOF'
#if !defined(SIMDE_X86_AVX512F_NATIVE) || !defined(SIMDE_X86_AVX512BW_NATIVE) || \
    !defined(SIMDE_X86_AVX512DQ_NATIVE) || !defined(SIMDE_X86_AVX512VL_NATIVE)
#  error "SIMDe is not using native AVX-512 F+BW+DQ+VL"
#endif
void rsimd_probe(double *out, const double *in) {
  simde__m512d a = simde_mm512_loadu_pd(in);
  simde__mmask8 m = simde_mm512_cmp_pd_mask(a, a, SIMDE_CMP_UNORD_Q);
  a = simde_mm512_mask_blend_pd(m, simde_mm512_fmadd_pd(a, a, a), a);
  simde_mm512_storeu_pd(out, a);
}
EOF
      ;;
    *) return 1 ;;
  esac
}

pass=0
fail=0
skip=0

# Fields: tier | arch | flags | headers | enabled
while IFS='|' read -r tier arch flags headers enabled; do
  tier=$(echo "$tier" | tr -d ' \t')
  case "$tier" in '' | '#'*) continue ;; esac
  arch=$(echo "$arch" | tr -d ' \t')
  enabled=$(echo "$enabled" | tr -d ' \t')
  flags=$(echo "$flags" | sed 's/^[ \t]*//; s/[ \t]*$//')
  [ "$flags" = "-" ] && flags=
  headers=$(echo "$headers" | sed 's/^[ \t]*//; s/[ \t]*$//')
  [ "$enabled" = yes ] || continue

  case " $headers" in
    *" x86/"* | *" @avx512"*) ;;
    *) continue ;;
  esac
  label="$tier ($arch)"

  cc=$(compiler_for "$arch")
  if [ "$cc" = none ] || ! command -v "$cc" > /dev/null 2>&1; then
    echo "SKIP $label: no compiler for $arch (tried '$cc'; set RSIMD_CC_$arch)"
    skip=$((skip + 1))
    continue
  fi

  tu="$tmp/check_${tier}_$arch.c"
  : > "$tu"
  for h in $headers; do
    if [ "$h" = "@avx512" ]; then
      grep '^#include "x86/avx512/' "$probe" >> "$tu"
    else
      echo "#include \"$h\"" >> "$tu"
    fi
  done
  if ! tier_body "$tier" >> "$tu"; then
    echo "FAIL $label: no probe body defined for this tier"
    fail=$((fail + 1))
    continue
  fi

  # shellcheck disable=SC2086 # $flags is a word list
  if "$cc" -std=c11 -O2 -Wall -Wextra -pedantic -Werror $flags -I"$simde" \
    -c "$tu" -o "$tmp/check.o" > "$tmp/log" 2>&1; then
    echo "PASS $label: $cc${flags:+ $flags}"
    pass=$((pass + 1))
  else
    echo "FAIL $label: $cc${flags:+ $flags}"
    sed 's/^/    /' "$tmp/log" | head -40
    fail=$((fail + 1))
  fi
done < "$tiers"

echo "SIMDe subset check: $pass passed, $fail failed, $skip skipped"
[ "$fail" -eq 0 ]
