# Contributing to rsimd

## Development dependencies

You need R (>= 4.3) and a C11 compiler. Install the packages used by the tests,
benchmarks and vignettes into your user library:

```sh
mkdir -p "$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"
Rscript -e 'install.packages(c("testthat", "bit64", "bench", "knitr", "rmarkdown"),
  lib = Sys.getenv("R_LIBS_USER"), repos = "https://cloud.r-project.org")'
```

Optional: `lintr` (configured by `.lintr`) and `styler` (run `Rscript tools/style.R`).

## Documentation and NAMESPACE

`NAMESPACE` and `man/*.Rd` are hand-written; roxygen2 is not used. When you add an
exported function:

- add an `export()` line to `NAMESPACE`, kept sorted within its `# <family>` section;
- add an `\alias` and `\usage` entry to the Rd file for its function family;
- use the shared macros in `man/macros/rsimd.Rd` instead of copying their wording;
- add a `NEWS.md` bullet and tests.

`tests/testthat/test-namespace-docs.R` checks that `NAMESPACE` and the Rd files agree.

## Checking the package

```sh
R CMD build .
_R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --as-cran rsimd_*.tar.gz
```

The check must finish without errors or warnings, and with no notes other than those
CRAN raises for every new submission.

To check that the C code compiles without warnings, add the following to
`~/.R/Makevars` before installing:

```make
CFLAGS = -g -O2 -Wall -Wextra -pedantic
```

## C language server

The C code is set up for [clangd](https://clangd.llvm.org/), which Claude Code's
LSP tool and the VS Code clangd extension use. `sh tools/compile_commands.sh`
writes `compile_commands.json` (not committed) without building anything; it runs
`./configure` first if `src/rsimd_config.h` is missing, for example after `cleanup`.
Every tier file is listed with its own flags and target, so on an arm64 machine the
x86 tiers are still indexed on their native SIMDe and SLEEF paths (and the other way
round); that needs the other architecture's cross compiler for its libc headers
(`gcc-x86-64-linux-gnu` or `gcc-aarch64-linux-gnu`). The kernel sources and headers
only compile inside a tier file, so the script gives each one a tier file's command
and a generated prelude (under `.cache/rsimd-lsp/`) holding what that tier file
defines and includes before it. The kernel sources are seen as the host's baseline
tier; set `RSIMD_LSP_TIER=avx2` (or another tier) to see them as that one. Re-run
the script after changing `tools/tiers.txt` or an `#include` line. `.clangd`
silences diagnostics in the vendored headers.

The devcontainer installs clangd and the cross compiler, runs the script, and
installs the `clangd-lsp` Claude Code plugin.

## Instruction-set tiers and the build

Kernels are compiled once per instruction-set tier. `src/tier_<tier>.c` is a
short file that defines `RSIMD_TIER`, includes `src/kernels/common.inc.h`,
then the kernel sources, and finally `src/kernels/table.inc.c`, which defines
the tier's dispatch table; it is compiled with that tier's flags (`-mavx2 -mfma`
for `avx2`, `-march=armv8-a+sve` for `sve` ...). Code in a tier file runs only
after a run-time check that the CPU supports the tier, and everything outside
the tier files is compiled with R's flags alone. `-march=native` is never used.

`./configure` (also run by `configure.win`) reads the tier table in
`tools/tiers.txt`, asks R for its compiler and flags, and test-compiles a probe
for each tier of the target. It writes `src/Makevars` (or `src/Makevars.win`)
from the `.in` template, with an explicit rule for every tier that compiles,
and `src/rsimd_config.h` with `RSIMD_HAVE_<TIER>` for each of them. For a tier
whose row names a SLEEF header (last column), it also test-compiles a call into
that header and, if it compiles, adds `-DRSIMD_HAVE_SLEEF_<TIER>=1` to the
tier's rule; otherwise the tier is built without SLEEF. Every tier file is
compiled with `-ffp-contract=off`. Tiers that fail are dropped with a message; the `none` tier and the platform baseline
(`sse2` on x86-64, `neon` on arm64) must compile. SVE tiers are built whenever
the compiler supports them, whatever the build machine's CPU. `cleanup`
removes the generated files. Useful environment variables:

- `RSIMD_DISABLE_TIERS=avx512,sve2` skips tiers;
- `RSIMD_DISABLE_SLEEF=1` stops tiers using SLEEF elementary functions (they
  then come from the `none` tier's libm);
- `RSIMD_CONFIGURE_VERBOSE=1` shows the compiler output of failed probes;
- `RSIMD_DEBUG=1` makes a debug build (`-DRSIMD_DEBUG`, extra checks at load);
- `RSIMD_TEST_HOLE=<slot>` leaves one dispatch slot empty in every tier but
  `none`, for testing the fill-down described below.

## Kernels and dispatch

Every operation has one slot per element type in `struct rsimd_kernels`,
generated from the one-line-per-slot X-macro list in `src/kernel_list.h`. To
add an operation, add its `RSIMD_OP(slot, return_type, (arguments))` line
there and write the kernel, named `RSIMD_KERNEL(slot)`, once in a
`src/kernels/<family>.inc.c` file; every tier then compiles it with its own
vector layer. A tier with no useful version of a kernel writes
`#define RSIMD_SKIP_<slot> 1` in the kernel source instead. The `none` tier
cannot skip: a missing `none` kernel is a compile error.

`src/dispatch.c` builds, at load, a resolved copy of each available tier's
table in which every empty slot is taken from the next lower available tier
(`avx512 > avx2 > sse2 > none`, `sve2 > sve > neon > none`). `simd_use()`
makes one of them active, and the `.Call` entry points call kernels through
`rsimd_active`. The internal `simd_kernel_tiers()` shows which tier each slot
of a table really runs. Every exported compute function starts with
`.sync_impl()`, so that a change to `options(rsimd.impl =)` takes effect.

`src/kernels/common.inc.h` documents the vector layer that kernels are written
against: types such as `rsimd_vf64`, operations such as `rsimd_vf64_add`, and
predicated loads and stores for loop tails. It maps to SIMDe on the x86 and NEON
tiers, to Arm SVE intrinsics on the SVE tiers and to plain C on `none`. On tiers
built with SLEEF it also wraps SLEEF's elementary functions as
`rsimd_sleef_<f>()` and states their accuracy policy; math kernels must fall
back to C99 libm where `RSIMD_HAVE_SLEEF` is not defined. Maintainer checks:

```sh
sh tools/check_vector_layer.sh   # run the layer, the kernels and the SLEEF wrappers on every tier
sh tools/check_tier_symbols.sh   # after building in place: tier symbols end in _<tier>
sh tools/check_build.sh install.log  # tier flags in an R CMD INSTALL log, no exported Sleef_*
```

`tools/check_vector_layer.sh` uses cross compilers and `qemu-user` when present
(Debian: `gcc-x86-64-linux-gnu`, `gcc-i686-linux-gnu`,
`gcc-arm-linux-gnueabihf`, `qemu-user`). It runs the SVE tiers at vector lengths
128 to 2048 bits; AVX-512 is compile-only because qemu cannot emulate it.

To add a tier, add its rows to `tools/tiers.txt` and a probe to `configure`,
the tier id to `src/tiers.h` and `src/tiers.c`, its layer to
`src/kernels/common.inc.h`, its table and place in the preference order to
`src/dispatch.c`, and copy an existing `src/tier_<tier>.c`.

## Entry points and R vectors

`.Call` entry points (`src/api_*.c`) never touch R vector data directly. They
go through the access layer in `src/rvec.h`, which classifies inputs, reads
ALTREP objects such as `1:n` in chunks without expanding them, applies the
length-1 broadcast rule, runs the chunk loop (`RSIMD_FOREACH_CHUNK`) with
interrupt checks every 2^20 elements, allocates results and applies the
attribute policy. Kernels see only the plain C types of `src/kernel_types.h`.
A grep, also run in CI, keeps raw accessors (`REAL()`, `DATAPTR_OR_NULL()`,
`XLENGTH()`, `R_CheckUserInterrupt()` ...) inside `src/rvec.c`:

```sh
sh tools/lint_c.sh
```

Setting `RSIMD_DEBUG_STRIDE=<n>` before loading the package shrinks the chunk
and interrupt stride to `n` elements, to test chunking on small inputs.

## Vendored code

`src/vendor/simde` is a pruned copy of the SIMDe headers, generated by
`tools/vendor_simde.sh` from the pinned commit and the header list in
`tools/simde_probe.c`; do not edit it by hand. After changing either, run:

```sh
sh tools/vendor_simde.sh        # regenerate src/vendor/simde and inst/COPYRIGHTS
sh tools/check_simde_subset.sh  # compile each SIMDe tier with -Werror
sh tools/check_size.sh          # vendored and tarball size budgets
```

`tools/check_simde_subset.sh` compiles tiers for another architecture only when a
cross compiler is available (for example `x86_64-linux-gnu-gcc` on arm64); other
tiers are reported as skipped. See `src/vendor/simde/README-rsimd.md` for details.

`src/vendor/sleef` holds SLEEF inline headers generated by `tools/vendor_sleef.sh`
(needs CMake and network access) from the pinned SLEEF tag; do not edit them by
hand. The script generates the x86-64 and arm64 headers on Linux, natively for the
machine's architecture and with a cross compiler (`gcc-x86-64-linux-gnu` or
`gcc-aarch64-linux-gnu`) for the other; both routes give identical files. The
manually triggered `vendor-sleef` GitHub Actions workflow runs it on both
architectures and checks the result against the vendored files. After updating:

```sh
sh tools/vendor_sleef.sh         # regenerate src/vendor/sleef and inst/COPYRIGHTS
sh tools/check_vector_layer.sh   # SLEEF wrappers against long double libm, every tier
sh tools/check_size.sh           # vendored, tarball and library size budgets
```

See `src/vendor/sleef/README-rsimd.md` for what is changed from SLEEF's output and
the known limits of its functions.
