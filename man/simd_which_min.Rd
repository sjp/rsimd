\name{simd_which_min}
\alias{simd_which_min}
\alias{simd_which_max}
\alias{simd_which_min_abs}
\alias{simd_which_max_abs}
\title{Where the Minimum or Maximum Is}
\usage{
simd_which_min(x)

simd_which_max(x)

simd_which_min_abs(x)

simd_which_max_abs(x)
}
\arguments{
\item{x}{A double, integer, logical, raw or \pkg{bit64} \code{integer64}
vector; for \code{simd_which_min_abs()} and \code{simd_which_max_abs()} a
double, integer, logical, complex or integer64 one.}
}
\description{
The index of the first minimum or maximum of a vector, like
\code{\link{which.min}()} and \code{\link{which.max}()}, or of its smallest
or largest magnitude, like \code{which.min(abs(x))} and
\code{which.max(abs(x))} without the vector \code{abs(x)}, computed with SIMD
instructions.
}
\details{
As in base R, missing values (\code{NA} and \code{NaN}) are ignored, the
first of several equal extremes is returned, \code{0} and \code{-0} are
equal (\code{simd_which_max(c(0, -0))} is 1), and logical values count as 0
and 1. The results are exact and the same for every implementation.

The minimum or maximum is found first and then the first element equal to
it, so the data are read up to twice. The precision mode does not apply.
All implementations have their own kernels, except for raw vectors, and on
32-bit ARM for double vectors, which use the scalar \code{"none"} code.
Integer64 vectors compare as 64-bit integers (\pkg{bit64} has no
\code{which.min} method; see \link{rsimd-integer64}).

\code{simd_which_max_abs(x)} and \code{simd_which_min_abs(x)} compare
magnitudes with the same rules: \code{simd_which_max_abs(c(-3, 3))} is 1,
and \code{x[simd_which_max_abs(x)]} is the element of largest magnitude
with its sign. For complex \code{x} they compare the moduli that
\code{\link{simd_abs}()} computes, as \code{\link{simd_max_abs}()} does,
which can differ from \code{\link{Mod}(z)} in the last bit: two elements
whose moduli are that close may be ordered differently from
\code{which.max(Mod(z))}. Complex moduli are computed in blocks and
compared once each, so the data are read once.
}
\value{
An integer of length 1, or \code{integer(0)} if \code{x} is empty or has
only missing values. For vectors longer than \code{.Machine$integer.max}
the index is a double.
}
\seealso{
\code{\link{which.min}}, \code{\link{simd_min}}, \code{\link{simd_max_abs}},
\link{rsimd-integer64}.
}
\examples{
simd_which_max(c(1, 3, 3))
simd_which_min(c(NA, 2, NaN, 1))
simd_which_max(numeric(0))
simd_which_max_abs(c(2, -5, 5))
simd_which_min_abs(c(-1L, NA, 4L))
}
