\name{simd_min}
\alias{simd_min}
\alias{simd_max}
\alias{simd_range}
\title{Minimum, Maximum and Range of Vector Elements}
\usage{
simd_min(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_max(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_range(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x}{A double, integer, logical or \pkg{bit64} \code{integer64} vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The smallest or largest element of a vector, or both, like
\code{\link{min}()}, \code{\link{max}()} and \code{\link{range}()} for a
single vector, computed with SIMD instructions.
}
\details{
The results are exact and the same for every implementation, including the
sign of a zero result: like base R, \code{simd_min(c(0, -0))} is \code{0}
and \code{simd_min(c(-0, 0))} is \code{-0} (the first zero wins).

Without \code{na.rm}, a missing value makes the result \code{NA} (if there
is an \code{NA}) or \code{NaN}, wherever it is; \code{simd_range()} then
returns it twice. When nothing is left after removing missing values (or
\code{x} is empty), \code{simd_min()} returns \code{Inf} and
\code{simd_max()} \code{-Inf}, as doubles, with base R's warning;
\code{simd_range()} returns both, with both warnings.

\rsimdSemantics

The precision mode does not apply. All implementations have their own
kernels for these functions, except that on 32-bit ARM the double versions
use the scalar \code{"none"} code. Raw and complex vectors are rejected, as
by base R. Integer64 vectors give integer64 results; for no values these
are the largest and smallest integer64 values, with \pkg{bit64}'s warnings
(see \link{rsimd-integer64}).
}
\value{
\code{simd_min()} and \code{simd_max()}: a double for double \code{x}, an
integer for integer and logical \code{x} (but see above for empty input).
\code{simd_range()}: a vector of length 2 of the same type,
\code{c(simd_min(x), simd_max(x))}.
}
\seealso{
\code{\link{min}}, \code{\link{range}}, \code{\link{simd_which_min}},
\link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_min(c(3, 1, 2))
simd_max(c(3L, NA, 2L), na.rm = TRUE)
simd_range(c(TRUE, FALSE))
simd_min(c(1, NaN, NA)) # NA: NA wins over NaN
}
