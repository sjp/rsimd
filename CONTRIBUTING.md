# Contributing to rsimd

## Development dependencies

You need R (>= 4.3) and a C11 compiler. Install the packages used by the tests and
benchmarks into your user library:

```sh
mkdir -p "$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"
Rscript -e 'install.packages(c("testthat", "bit64", "bench", "callr"),
  lib = Sys.getenv("R_LIBS_USER"), repos = "https://cloud.r-project.org")'
```

Optional: `lintr` (configured by `.lintr`) and `styler` (run `Rscript tools/style.R`).

The vignettes in `vignettes/` are Sweave (`.Rnw`) files built with R's own
`utils::Sweave` engine, so they need no packages, but building them (`R CMD build`)
needs a LaTeX installation with `pdflatex`. Use `R CMD build --no-build-vignettes` without one.

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
R CMD build --compact-vignettes=both .
_R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --as-cran rsimd_*.tar.gz
```

`--compact-vignettes=both` shrinks the vignette PDFs with Ghostscript and qpdf; without it
`--as-cran` warns that they could be made much smaller.

The check must finish without errors or warnings, and with no notes other than those
explained in `cran-comments.md`.

To check that the C code compiles without warnings, add the following to
`~/.R/Makevars` before installing (or put it in a scratch file named by
`R_MAKEVARS_USER`):

```make
CFLAGS = -g -O2 -Wall -Wextra -pedantic -Wstrict-prototypes -Wconversion
```

The package's own files must compile without warnings with GCC and clang. The
vendored SIMDe and SLEEF headers are kept identical to upstream: they warn under
`-Wextra` and `-Wconversion`, though not under CRAN's `-Wall -pedantic`. R's own
`R_ext/Boolean.h` also warns under `-pedantic` before C23.

Other checks to run before a release:

```sh
sh ./configure && sh tools/check_makevars.sh   # no forbidden flags or variable overrides
Rscript tools/check_spelling.R                 # spelling, against inst/WORDLIST
sh tools/check_size.sh                         # tarball, installed and library sizes
shellcheck -s sh configure configure.win cleanup cleanup.win
checkbashisms configure configure.win cleanup cleanup.win
```

## Continuous integration

The GitHub Actions workflows in `.github/workflows/` call POSIX shell scripts in
`tools/ci/`, so each job can be run locally. Every script reports one row per run
(label, platform, result, duration) to the job summary, or to the terminal outside
Actions.

On every pull request and every push to `main` (except changes only to `issues/`,
`bench/`, `.devcontainer/`, `CONTRIBUTING.md`, `cran-comments.md` and
`LICENSE.md`, which do not ship):

| Workflow | Pull requests | Pushes to `main` add |
|----------|---------------|----------------------|
| `fast` | The first signal, meant to take under 5 minutes: on Linux x86-64 and arm64, install, then the whole suite at quick lengths (`RSIMD_TEST_SUBSET=quick`) with the kernel tests looping over every available tier, then the dispatch tests once per tier. | |
| `R-CMD-check` | `R CMD check --as-cran` (warnings fail; the tests loop over every tier) on Linux x86-64 (R release, devel, oldrel-1), Linux arm64, macOS arm64 and Windows x86-64. Each job then runs the dispatch tests once per tier with that tier as the default (`RSIMD_IMPL`), fails if `none`, or `sse2` and `avx2` on x86-64, or `neon` on arm64, is missing, and prints the runner's CPU features and tiers. Only Linux x86-64 R release builds the vignettes and the PDF manual; the others neither build nor check the vignettes. The C lint, `tools/check_makevars.sh` and `tools/check_build.sh` on Linux. | macOS x86-64; the vignettes on every platform; builds without SLEEF and with a compiler that rejects AVX-512. |
| `size-check` | `tools/check_size.sh`: tarball and installed package (library stripped) under 5 MB, no installed directory but `libs` over 1 MB, vendored code and library budgets. | |
| `coverage` | | covr, uploaded to Codecov (needs the `CODECOV_TOKEN` secret); informational. |

The per-platform check job is `check-os.yaml`, which `R-CMD-check` calls. It caches
TinyTeX with the packages listed in `tools/ci/latex-packages.txt`.

On a schedule, and on demand with "Run workflow":

| Workflow | When | What it does |
|----------|------|--------------|
| `sanitizers` | nightly | R CMD check, the suite (every tier) and the dispatch tests per tier on R-hub's `gcc-asan`, `clang-asan` and `clang-ubsan` R-devel containers; any sanitizer report fails. |
| `valgrind` | weekly | `R CMD check --use-valgrind` (examples), then the quick subset on the `avx2` and `none` tiers under memcheck with leak checking. |
| `noLD` | weekly | R CMD check, the suite (every tier) and the dispatch tests per tier on R-hub's `nold` container (R without long double). |
| `emulation-sde` | weekly | The `avx512` tier under Intel SDE (Sapphire Rapids and Ice Lake models, `tier_emulation` subset), and `auto` choosing `sse2` on a Merom model (quick subset). |
| `emulation-qemu` | weekly | The `sve` and `sve2` tiers under `qemu-aarch64` at 256- and 512-bit vectors (`tier_emulation` subset), and `auto` choosing `neon` on a Cortex-A72 model (quick subset). |
| `cran-incoming` | release branches (`release/**`), `v*` tags | `R CMD check --as-cran` with CRAN's remote incoming checks (URLs, maintainer), the PDF manual and the spelling check; and a check without the suggested packages bit64 and bench. Errors and warnings fail. |
| `benchmarks` | weekly | `bench/run.R` on Linux x86-64, Linux arm64 and macOS arm64; results uploaded as an artifact and shown in the job summary (see `bench/README.md`). Informational; never fails on timings. |
| `math-paths` | on demand only | `bench/math_paths.R` on Linux arm64 (Neoverse) and macOS arm64, for the package as built and for a build with `RSIMD_NO_MATH_LIBM=1` (SLEEF for every function): the functions the arm64 tiers hand to the C math library and the complex functions on their branch cuts, per tier and accuracy mode. Checks the libm list in `src/kernels/math.inc.c` on those cores; informational. |

The emulated runs use the reduced test subsets described in `tests/README.md`
(`RSIMD_TEST_SUBSET`, `RSIMD_TEST_TIERS`). To run the jobs locally, install the
package first (`R CMD INSTALL .`), then from the package root:

```sh
# The suite once (every tier), then the dispatch tests once per tier, checking the tiers
# every machine of this architecture must have. With RSIMD_TEST_SUBSET=quick, as the fast job.
sh tools/ci/test_tiers.sh --logs logs

# One run, optionally under an emulator; --expect-impl and --expect-tiers fail early
# if the session did not select or detect what it should.
sh tools/ci/run_tests.sh --impl none

# SVE at a 256-bit vector length on an arm64 Linux machine (Debian: apt install qemu-user).
RSIMD_TEST_SUBSET=tier_emulation RSIMD_TEST_TIERS=sve \
  sh tools/ci/run_tests.sh --impl sve --expect-impl sve --expect-tiers "sve sve2" \
    --runner "qemu-aarch64 -cpu max,sve-default-vector-length=32"
# No SVE: auto must choose neon.
RSIMD_TEST_SUBSET=quick sh tools/ci/run_tests.sh --expect-impl neon --runner "qemu-aarch64 -cpu cortex-a72"

# AVX-512 under Intel SDE on an x86-64 Linux machine. SDE may not be redistributed:
# get_sde.sh downloads it from Intel (SDE_URL and SDE_SHA256 pick another release).
sde=$(sh tools/ci/get_sde.sh .sde)
RSIMD_TEST_SUBSET=tier_emulation RSIMD_TEST_TIERS=avx512 \
  sh tools/ci/run_tests.sh --impl avx512 --expect-impl avx512 --runner "$sde -spr --"
# An SSE2-class CPU: auto must choose sse2, and SDE aborts on any AVX instruction.
RSIMD_TEST_SUBSET=quick sh tools/ci/run_tests.sh --expect-impl sse2 --runner "$sde -mrm --"

# Quick subset under valgrind.
RSIMD_TEST_SUBSET=quick sh tools/ci/run_tests.sh \
  --runner "valgrind --tool=memcheck --leak-check=full --track-origins=yes --error-exitcode=1"
```

With `--runner`, R's own executable is started directly under the emulator, with the
environment R's front-end script would set; tests that start a subprocess then run it
natively. For sanitizers without a sanitizer build of R, see `tests/README.md`; in a
sanitizer container (`docker run ghcr.io/r-hub/containers/gcc-asan`), run
`test_tiers.sh --logs logs` and then `sh tools/ci/scan_sanitizer_logs.sh logs/*.log`.

Not covered by CI: Windows on arm64 (no hosted runners), and 32-bit Linux (i686,
armv7), which could be added as Docker jobs later.

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
