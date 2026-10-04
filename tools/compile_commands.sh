#!/bin/sh
# Writes compile_commands.json at the package root for clangd (the C
# language server used by Claude Code's LSP tool and the VS Code clangd
# extension). Nothing is compiled; the commands mirror src/Makevars.
#
#   - src/*.c is compiled with R's include directory and the package's -I
#     flags, for the host target.
#   - src/tier_<tier>.c gets its tier's flags from tools/tiers.txt, the
#     -ffp-contract=off that configure adds, and -DRSIMD_HAVE_SLEEF_<TIER>=1
#     if the tier has a SLEEF header. Every x86_64 and aarch64 tier is listed
#     whatever the host is, with --target=<arch>-linux-gnu, so the x86 tiers
#     are indexed on their native SIMDe and SLEEF paths on an arm64 machine
#     and the other way round. This needs the other architecture's libc
#     headers (the gcc-x86-64-linux-gnu or gcc-aarch64-linux-gnu package);
#     without them clangd reports errors in system headers for those tiers.
#   - The kernel sources (src/kernels/*.inc.c) and the headers only compile
#     inside a tier file, after what it defines and includes before them; a
#     vector layer header is even included from the middle of an #elif chain
#     in kernels/common.inc.h. Each one is therefore given the command of a
#     tier file that includes it -- $RSIMD_LSP_TIER (by default the host's
#     baseline: sse2 on x86_64, neon on aarch64) when it does, else the
#     first such tier in tools/tiers.txt -- plus -include of a generated
#     prelude: the text of every file on the include chain up to the line
#     that includes the next one, cut before a conditional that is still
#     open there, with #line markers so that clangd points into the real
#     files. The preludes are written to .cache/rsimd-lsp/. Headers that no
#     tier file includes are compiled on their own for the host.
#   - tools/vector_layer_test.c is compiled for $RSIMD_LSP_TIER, as
#     tools/check_vector_layer.sh does, and tools/simde_probe.c for avx512,
#     whose headers it lists.
#
# src/rsimd_config.h comes from ./configure, which is run first if the
# header is missing. Which tier file includes which header is found with
# the C compiler R uses (cc -MM).
#
# Usage: sh tools/compile_commands.sh
# Re-run it after changing tools/tiers.txt or the files' #include lines.

set -eu

LC_ALL=C
export LC_ALL

root=$(cd "$(dirname "$0")/.." && pwd -P)
cd "$root"

r_include=$(Rscript -e 'cat(R.home("include"))' 2> /dev/null) || r_include=
if [ -z "$r_include" ] || [ ! -f "$r_include/Rinternals.h" ]; then
  echo "compile_commands.sh: cannot find R's include directory (is Rscript on PATH?)" >&2
  exit 1
fi
cc=$(R CMD config CC)

if [ ! -f src/rsimd_config.h ]; then
  echo "compile_commands.sh: src/rsimd_config.h is missing, running ./configure"
  ./configure > /dev/null
fi

case "$(uname -m)" in
  x86_64 | amd64) host_arch=x86_64 baseline=sse2 ;;
  aarch64 | arm64) host_arch=aarch64 baseline=neon ;;
  *)
    echo "compile_commands.sh: unsupported host architecture $(uname -m)" >&2
    exit 1
    ;;
esac
lsp_tier=${RSIMD_LSP_TIER:-$baseline}

src="$root/src"
pkg_flags="-I$src -I$src/vendor/simde -I$src/vendor/sleef -isystem $r_include"
host_flags="--target=$host_arch-linux-gnu $pkg_flags"

# One "tier|arch|flags" line for the first x86_64 or aarch64 row of every
# enabled tier (sse2 and neon also have i686 and armv7 rows, sve2 has a
# fallback row). none is built for every target and gets the host's.
tier_rows=$(sed -e 's/#.*//' -e '/^[[:space:]]*$/d' tools/tiers.txt |
  awk -F'|' -v host="$host_arch" '
    function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
    {
      tier = trim($1); arch = trim($2); flags = trim($4)
      enabled = trim($6); sleef = trim($7)
      if (enabled != "yes" || seen[tier]) next
      if (arch == "all") arch = host
      if (arch != "x86_64" && arch != "aarch64") next
      seen[tier] = 1
      if (flags == "-") flags = ""
      if (sleef != "-" && (getline line < ("src/vendor/sleef/" sleef)) > 0)
        flags = flags " -DRSIMD_HAVE_SLEEF_" toupper(tier) "=1"
      if (system("test -f src/tier_" tier ".c") == 0) print tier "|" arch "|" flags
    }')
tiers=$(echo "$tier_rows" | cut -d'|' -f1)

tier_flags() {
  echo "$tier_rows" | awk -F'|' -v t="$1" \
    '$1 == t { print "--target=" $2 "-linux-gnu " $3 " -ffp-contract=off" }'
}

case " $(echo $tiers) " in
  *" $lsp_tier "*) ;;
  *)
    echo "compile_commands.sh: RSIMD_LSP_TIER: unknown tier '$lsp_tier'" >&2
    exit 1
    ;;
esac

cache="$root/.cache/rsimd-lsp"
rm -rf "$cache"
mkdir -p "$cache"

# The package files each tier file includes, one "src/..." path per line.
# Missing system headers (<arm_sve.h> on x86_64 ...) and the #error that
# common.inc.h raises for a tier file compiled without its flags do not stop
# the listing.
for t in $tiers; do
  (cd src && $cc -MM -MG -I. -Ivendor/simde -Ivendor/sleef -isystem "$r_include" \
    "tier_$t.c" 2> /dev/null || true) |
    tr ' \\' '\n\n' | grep -E '^[a-z_0-9/.]+\.(h|c)$' | grep -v '^vendor/' |
    sed 's|^|src/|' | sort -u > "$cache/deps_$t"
done

# write_prelude <tier> <file> <out>: writes the prelude for <file> in the
# tier file of <tier> to <out>. A depth-first walk in include order over the
# tier's package files finds the chain from the tier file to <file>; every
# file on it contributes its lines before the include of the next one,
# without its include guard's #ifndef (its #define stays, as after a real
# include) and without the outermost conditional still open there.
write_prelude() {
  awk -v target="$2" -v start="src/tier_$1.c" -v root="$root" '
    function dir(p) { sub(/\/[^\/]*$/, "", p); return p }
    function resolve(from, name,   p) {
      p = dir(from) "/" name
      if (p in dep) return p
      p = "src/" name
      return (p in dep) ? p : ""
    }
    function visit(f, d,   n, i, line, name, p) {
      n = 0
      while ((getline line < f) > 0) text[d, ++n] = line
      close(f)
      for (i = 1; i <= n; i++) {
        line = text[d, i]
        if (line !~ /^[ \t]*#[ \t]*include[ \t]*"/) continue
        name = line
        sub(/^[^"]*"/, "", name)
        sub(/".*/, "", name)
        p = resolve(f, name)
        if (p == "" || (p in seen)) continue
        seen[p] = 1
        if (p == target || visit(p, d + 1)) {
          chain[d] = f; upto[d] = i
          if (d > nchain) nchain = d
          return 1
        }
      }
      return 0
    }
    function directive(s) {
      if (s !~ /^[ \t]*#/) return ""
      sub(/^[ \t]*#[ \t]*/, "", s)
      return s
    }
    { dep[$0] = 1 }
    END {
      if (!visit(start, 1)) exit 1
      for (d = 1; d <= nchain; d++) {
        guard = 0; first = ""
        for (i = 1; i < upto[d]; i++) {
          s = directive(text[d, i])
          if (s == "") continue
          if (first == "") { first = s; fi = i; continue }
          if (first ~ /^ifndef[ \t]/ && s ~ /^define[ \t]/) {
            split(first, a, /[ \t]+/); split(s, b, /[ \t]+/)
            if (a[2] == b[2]) guard = fi
          }
          break
        }
        cut = upto[d]; depth = 0
        for (i = 1; i < upto[d]; i++) {
          if (i == guard) continue
          s = directive(text[d, i])
          if (s ~ /^if/) { if (depth++ == 0) open = i }
          else if (s ~ /^endif/) depth--
        }
        if (depth > 0) cut = open
        printf "#line 1 \"%s/%s\"\n", root, chain[d]
        for (i = 1; i < cut; i++) print (i == guard ? "" : text[d, i])
      }
    }' "$cache/deps_$1" > "$3"
}

out="$root/compile_commands.json"
list="$cache/entries"
: > "$list"

# entry <file> <flags>: one "file<TAB>flags" line, turned into JSON below.
entry() {
  printf '%s\t%s\n' "$1" "$2" >> "$list"
}

for f in src/*.c; do
  case "$f" in src/tier_*.c) continue ;; esac
  entry "$f" "$host_flags"
done
for t in $tiers; do
  entry "src/tier_$t.c" "$(tier_flags "$t") $pkg_flags"
done

for f in src/kernels/*.inc.c src/kernels/*.h src/*.h; do
  tier=
  for t in $lsp_tier $tiers; do
    if grep -qx "$f" "$cache/deps_$t"; then
      tier=$t
      break
    fi
  done
  case "$f" in *.h) lang="-x c-header" ;; *) lang= ;; esac
  if [ -z "$tier" ]; then
    entry "$f" "$lang $host_flags"
    continue
  fi
  prelude="$cache/$(echo "${f#src/}" | tr / _).prelude.h"
  write_prelude "$tier" "$f" "$prelude"
  entry "$f" "$lang $(tier_flags "$tier") $pkg_flags -I$src/kernels -include $prelude"
done

entry tools/vector_layer_test.c "$(tier_flags "$lsp_tier") $pkg_flags -DRSIMD_TIER=$lsp_tier"
entry tools/simde_probe.c "-std=c11 --target=x86_64-linux-gnu $(echo "$tier_rows" |
  awk -F'|' '$1 == "avx512" { sub(/ -DRSIMD_HAVE_SLEEF.*/, "", $3); print $3 }') -I$src/vendor/simde"

awk -F'\t' -v dir="$root" '
  function q(s) { gsub(/\\/, "\\\\", s); gsub(/"/, "\\\"", s); return "\"" s "\"" }
  BEGIN { print "[" }
  {
    n = split("clang " $2 " -c " dir "/" $1, args, / +/)
    line = "  {\"directory\": " q(dir) ", \"file\": " q(dir "/" $1) ", \"arguments\": ["
    for (i = 1; i <= n; i++) line = line (i > 1 ? ", " : "") q(args[i])
    printf "%s%s]}", (NR > 1 ? ",\n" : ""), line
  }
  END { print "\n]" }' "$list" > "$out"

echo "compile_commands.sh: wrote $(grep -c '"file"' "$out") entries to compile_commands.json" \
  "(kernel sources as the $lsp_tier tier)"
