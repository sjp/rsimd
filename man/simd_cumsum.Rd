\name{simd_cumsum}
\alias{simd_cumsum}
\alias{simd_cumprod}
\alias{simd_cummin}
\alias{simd_cummax}
\title{Cumulative Sums, Products, Minima and Maxima}
\usage{
simd_cumsum(x)

simd_cumprod(x)

simd_cummin(x)

simd_cummax(x)
}
\arguments{
\item{x}{A double, integer or logical vector; for \code{simd_cumsum()}
and \code{simd_cumprod()} also a complex vector, and except for
\code{simd_cumprod()} a \pkg{bit64} \code{integer64} vector.}
}
\description{
Running sums, products, minima and maxima of the elements of a vector, like
\code{\link{cumsum}()}, \code{\link{cumprod}()}, \code{\link{cummin}()} and
\code{\link{cummax}()}, computed with SIMD instructions.
}
\details{
The result types are base R's: \code{simd_cumprod()} always returns a
double; the others return a double for double \code{x} and an integer for
integer and logical \code{x}.

As in base R, everything from the first missing value onward is missing:
from the first \code{NaN} the results are \code{NaN}, and from the first
\code{NA} they are \code{NA} (\code{simd_cumsum(c(1, NaN, NA))} is
\code{1 NaN NA}). An integer \code{simd_cumsum()} that leaves the integer
range is \code{NA} from there on, with base R's warning
\code{"integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'"}.

Integer results, and all results of \code{simd_cummin()} and
\code{simd_cummax()}, are identical to base R's on every implementation,
including the sign of a zero (of equal values the later one is kept, so
\code{simd_cummax(c(0, -0))} ends in \code{-0}). Double sums and products
are computed a vector register at a time: the elements of a register are
combined in a tree before the running value is applied, so the results can
differ from base R's sequential ones in the last bits and between
implementations. In \code{"compensated"} precision mode
\code{simd_cumsum()} instead computes a sequential compensated running sum,
identical on every implementation and usually closer to the exact one, at
the speed of a scalar loop. \code{simd_cumprod()} ignores the precision
mode.

\code{simd_cumsum()} and \code{simd_cumprod()} of a complex vector are
complex and identical to base R's on every implementation, missing values
included: they are computed in order as base R computes them (each step of
\code{cumprod()} with the rounding of base R's build, see
\code{\link{simd_mul}()}), then base R's rule for missing values is
applied: from the first element with a missing part on, each part that is
missing at the end is \code{NA} where an \code{NA} has been seen and
\code{NaN} before that. Each step depends on the one before, so they are no
faster than base R.

The \code{"sve"} and \code{"sve2"} implementations use the \code{"neon"}
kernels for these functions, and on 32-bit ARM they use the scalar
\code{"none"} code.

\rsimdPrecision

\rsimdAttributes

Integer64 input gives integer64 results, computed sequentially on every
implementation; a \code{simd_cumsum()} that leaves the integer64 range is
\code{NA} from there on, with \pkg{bit64}'s warning (see \link{rsimd-integer64}).
}
\value{
A vector of the same length as \code{x}: double, integer or complex as
described in Details.
}
\seealso{
\code{\link{cumsum}}, \code{\link{simd_sum}}, \code{\link{simd_precision}},
\link{rsimd-integer64}.
}
\examples{
simd_cumsum(1:10)
simd_cumprod(c(1, 2, 3, 4))
simd_cummax(c(1, 3, 2, 5, 4))
simd_cumsum(c(1, NaN, 3, NA, 5))
}
