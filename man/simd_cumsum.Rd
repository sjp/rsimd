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
\item{x}{A double, integer or logical vector.}
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

The \code{"sve"} and \code{"sve2"} implementations use the \code{"neon"}
kernels for these functions, and on 32-bit ARM they use the scalar
\code{"none"} code.

\rsimdPrecision

\rsimdAttributes
}
\value{
A vector of the same length as \code{x}: double, or integer as described in
Details.
}
\seealso{
\code{\link{cumsum}}, \code{\link{simd_sum}}, \code{\link{simd_precision}}.
}
\examples{
simd_cumsum(1:10)
simd_cumprod(c(1, 2, 3, 4))
simd_cummax(c(1, 3, 2, 5, 4))
simd_cumsum(c(1, NaN, 3, NA, 5))
}
