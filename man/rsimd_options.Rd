\name{rsimd_options}
\alias{rsimd_options}
\alias{simd_options}
\alias{print.simd_options}
\alias{simd_na_check}
\title{Options and Environment Variables of rsimd}
\usage{
simd_options()

simd_na_check(value)
}
\arguments{
\item{value}{\code{TRUE} or \code{FALSE}.}
}
\description{
The options and environment variables that control \pkg{rsimd}.
\code{simd_options()} reports all the options at once;
\code{simd_na_check()} gets or sets \code{rsimd.na_check}.
}
\value{
\code{simd_options()}: a list of class \code{"simd_options"} with elements
\code{impl} (the requested implementation, as returned by
\code{\link{simd_use}()}), \code{current} (the tier in use, as returned by
\code{\link{simd_current}()}), \code{precision}, \code{math_accuracy} and
\code{na_check}. Its print method shows one line per option. An option
holding an invalid value makes it an error, as for the function that reads
that option.

\code{simd_na_check()}: without \code{value}, the current value of
\code{rsimd.na_check}. With \code{value}, the previous value, invisibly
(\code{TRUE}, the default, when the option held an invalid value); an
invalid \code{value} is an error and leaves the option unchanged.
}
\section{Options}{
\describe{
\item{\code{rsimd.impl}}{The implementation in use: \code{"auto"} (the
default) or a tier id such as \code{"avx2"} or \code{"none"}. Set it with
\code{\link{simd_use}()}, which checks the value; setting the option
directly takes effect on the next call of an \pkg{rsimd} function.}
\item{\code{rsimd.precision}}{How floating-point sums are accumulated:
\code{"fast"} (the default), \code{"pairwise"} or \code{"compensated"}. See
\code{\link{simd_precision}()}.}
\item{\code{rsimd.math_accuracy}}{Which variants of the SLEEF library the
SIMD implementations of the elementary functions use: \code{"accurate"}
(the default) or \code{"fast"}, which is faster and less accurate in the
last bits. See \code{\link{simd_math_accuracy}()}.}
\item{\code{rsimd.na_check}}{\code{TRUE} (the default) or \code{FALSE}. Used
when the \code{na_check} argument is \code{NULL} (its default):
\code{FALSE} promises that inputs contain no missing values, so functions
skip looking for them. See \code{simd_na_check()} and \sQuote{Missing
values and integer overflow} below.}
}
The options are set to their defaults when the package is loaded, unless
they are already set. An option already set to a value that is not valid
gives a warning naming it, and is set to its default. A value that is not
valid set after loading makes the next call of a function that reads the
option an error naming it and how to reset it: with \code{\link{simd_use}()},
\code{\link{simd_precision}()}, \code{\link{simd_math_accuracy}()} or
\code{simd_na_check()}, which replace an invalid value.
}
\section{Environment variables}{
\describe{
\item{\env{RSIMD_IMPL}}{Read when the package is loaded, if the option
\code{rsimd.impl} is not set: the initial implementation. An unknown or
unavailable value gives a warning and \code{"auto"} is used.}
\item{\env{RSIMD_CPU_FEATURES_MASK}}{For testing: a comma-separated list of
CPU feature names (see \code{\link{simd_cpu_features}()}) to treat as
absent, read when the package is loaded.}
\item{\env{RSIMD_DEBUG_STRIDE}}{For testing: the number of elements
processed between checks for user interrupts (rounded up to a multiple of
128), read when the package is loaded. Long inputs are processed in pieces
of this size, so a small value exercises that splitting. The default is
\eqn{2^{20}}{2^20}.}
\item{\env{RSIMD_DISABLE_TIERS}}{Read when the package is installed: a
comma-separated list of tiers not to compile, such as
\code{"avx512,sve2"}. The \code{"none"} tier cannot be disabled: naming it,
or a tier that does not exist, gives a warning and is ignored.}
\item{\env{RSIMD_DISABLE_SLEEF}}{Read when the package is installed: set to
\code{1} to build no tier with the SLEEF elementary functions; they then
come from the C library through \code{"none"}.}
\item{\env{RSIMD_CONFIGURE_VERBOSE}}{Read when the package is installed: set
to \code{1} to show the compiler output of the tests that fail, to see why
a tier was not built.}
}
}
\section{Missing values and integer overflow}{
\rsimdSemantics
}
\seealso{
\code{\link{simd_use}()}, \code{\link{simd_precision}()},
\code{\link{simd_math_accuracy}()},
\code{\link{simd_cpu_features}()}.
}
\examples{
simd_options()

simd_na_check()
old <- simd_na_check(FALSE)
simd_sum(c(1, 2, 3))
simd_na_check(old)
}
