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
