\name{simd_use}
\alias{simd_use}
\alias{simd_tiers}
\alias{simd_available}
\alias{simd_current}
\alias{simd_complex_variants}
\alias{simd_with_impl}
\title{Select the SIMD Implementation}
\usage{
simd_tiers()
simd_available()
simd_current()
simd_complex_variants()
simd_use(impl)
simd_with_impl(impl, expr)
}
\arguments{
\item{impl}{\code{"auto"} or one of the ids returned by \code{simd_tiers()}.}
\item{expr}{An expression to evaluate with \code{impl} selected.}
}
\description{
\pkg{rsimd} contains several implementations of its functions, one per
instruction set tier, and runs one of them at a time. These functions list
the tiers, report which one is in use and change it.
}
\details{
The tiers are \code{"none"} (plain C, available everywhere and the reference
the others are tested against), \code{"sse2"}, \code{"avx2"} (with FMA) and
\code{"avx512"} (F, BW, DQ and VL) on x86, and \code{"neon"}, \code{"sve"} and
\code{"sve2"} on arm. \code{"rvv"} and \code{"wasm128"} are reserved names
that no build provides yet.

A tier is available when the package was built with it and the CPU and
operating system support it. \code{"auto"} selects the best available tier:
\code{avx512}, \code{avx2}, \code{sse2}, \code{none} on x86 and \code{sve2},
\code{sve}, \code{neon}, \code{none} on arm.

Not every tier has its own version of every operation. When the selected
tier lacks one, the version of the next lower available tier is used, ending
with \code{"none"}. \code{simd_current()} therefore reports the selected
tier, while an individual operation may run at a lower one.

The selection is global. \code{simd_use()} sets it, together with the option
\code{rsimd.impl}. Setting that option directly also works: it takes effect
the next time an \pkg{rsimd} function is called. When the package is loaded,
the initial value comes from the option \code{rsimd.impl} if it is already
set, else from the environment variable \env{RSIMD_IMPL}, else it is
\code{"auto"}. An unknown or unavailable value there gives a warning and
\code{"auto"} is used instead.

\code{simd_with_impl()} is meant for tests and benchmarks: it selects
\code{impl}, evaluates \code{expr} and then restores the previous selection,
also when \code{expr} signals an error.

The implementation asked for, \code{"auto"} or a tier id, is the value of
the option \code{rsimd.impl}, which \code{simd_use()} keeps up to date.

\code{simd_complex_variants()} reports how complex multiplication, division
and \code{cumprod()} were matched to base R's rounding when the package was
loaded (see \code{\link{simd_mul}()}); it does not depend on the selected
tier.
}
\value{
\code{simd_tiers()}: a character vector of every tier id, including reserved
ones.

\code{simd_available()}: a character vector of the available tiers, best
first, always ending with \code{"none"}.

\code{simd_current()}: a character string, the selected tier.

\code{simd_complex_variants()}: a named character vector (\code{mul},
\code{div}, \code{cumprod}) describing how complex multiplication,
division and \code{cumprod()} reproduce the rounding of this R build: the
formula the kernels use, or \code{"scalar"} when no kernel variant matches and base R's operator is
used for every element, with \code{" (NA/NaN may differ)"} appended if
only the choice between \code{NA} and \code{NaN} could not be
reproduced.

\code{simd_use()}: the previously requested value, invisibly. It signals an
error, leaving the selection unchanged, if \code{impl} is unknown or not
available.

\code{simd_with_impl()}: the value of \code{expr}.
}
\seealso{
\code{\link{simd_cpu_features}} for the CPU features behind availability.
}
\examples{
simd_available()
simd_current()
getOption("rsimd.impl")
simd_complex_variants()

old <- simd_use("none")
simd_current()
simd_use(old)

simd_with_impl("none", simd_current())
}
