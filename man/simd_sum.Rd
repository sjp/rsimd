\name{simd_sum}
\alias{simd_sum}
\title{Sum of Vector Elements}
\usage{
simd_sum(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) before summing?}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The sum of the elements of a vector, like \code{\link{sum}()} for a single
vector, computed with SIMD instructions.
}
\details{
Integer and logical vectors are summed exactly. Double vectors are summed in
the mode set with \code{\link{simd_precision}()}, so the result can differ
from base R's in the last bits.

\rsimdSemantics

\rsimdPrecision
}
\value{
A double for double \code{x}. For integer and logical \code{x}, an integer,
or a double if the sum is outside the integer range, as in base R. The sum of
an empty vector is \code{0} (or \code{0L}).
}
\seealso{
\code{\link{sum}}, \code{\link{simd_precision}}, \link{rsimd_options}.
}
\examples{
simd_sum(c(1.5, 2.5, NA), na.rm = TRUE)
simd_sum(1:10)
simd_sum(c(TRUE, FALSE, TRUE))
simd_sum(rep(.Machine$integer.max, 3L)) # a double, as base R

old <- simd_precision("compensated")
simd_sum(c(1e16, 1, -1e16))
simd_precision(old)
}
