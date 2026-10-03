\name{simd_sum}
\alias{simd_sum}
\alias{simd_prod}
\alias{simd_mean}
\title{Sums, Products and Means of Vector Elements}
\usage{
simd_sum(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_prod(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_mean(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The sum, product or arithmetic mean of the elements of a vector, like
\code{\link{sum}()}, \code{\link{prod}()} and \code{\link{mean}()} for a
single vector, computed with SIMD instructions.
}
\details{
Integer and logical vectors are summed exactly. Double vectors are summed in
the mode set with \code{\link{simd_precision}()}, so the result can differ
from base R's in the last bits.

\code{simd_prod()} multiplies in double precision in every mode, keeping
several partial products. Base R multiplies in \code{long double} where the
platform has it, whose wider exponent range can avoid an intermediate
overflow or underflow: \code{prod(c(1e300, 1e300, 1e-300, 1e-300))} is 1
in base R, while \code{simd_prod()} gives \code{Inf} or 1 depending on the
implementation in use, which forms its partial products in its own order.

\code{simd_mean()} of a double vector is the sum divided by the number of
elements; in \code{"pairwise"} and \code{"compensated"} modes it then adds
the mean of the deviations from that first estimate (a second pass over the
data), as base R does. The mean of an integer or logical vector is its exact
sum divided once, in \code{long double} where available, as in base R.

\rsimdSemantics

\rsimdPrecision

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double versions use the scalar \code{"none"} code.
}
\value{
\code{simd_sum()}: a double for double \code{x}. For integer and logical
\code{x}, an integer, or a double if the sum is outside the integer range,
as in base R. The sum of an empty vector is \code{0} (or \code{0L}).

\code{simd_prod()} and \code{simd_mean()}: a double. The product of an
empty vector is 1 and its mean \code{NaN}.
}
\seealso{
\code{\link{sum}}, \code{\link{prod}}, \code{\link{mean}},
\code{\link{simd_precision}}, \link{rsimd_options}.
}
\examples{
simd_sum(c(1.5, 2.5, NA), na.rm = TRUE)
simd_sum(1:10)
simd_sum(c(TRUE, FALSE, TRUE))
simd_sum(rep(.Machine$integer.max, 3L)) # a double, as base R

simd_prod(1:10)
simd_mean(c(1, 2, NaN, 4), na.rm = TRUE)

old <- simd_precision("compensated")
simd_sum(c(1e16, 1, -1e16))
simd_precision(old)
}
