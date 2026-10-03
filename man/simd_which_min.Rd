\name{simd_which_min}
\alias{simd_which_min}
\alias{simd_which_max}
\title{Where the Minimum or Maximum Is}
\usage{
simd_which_min(x)

simd_which_max(x)
}
\arguments{
\item{x}{A double, integer, logical or raw vector.}
}
\description{
The index of the first minimum or maximum of a vector, like
\code{\link{which.min}()} and \code{\link{which.max}()}, computed with SIMD
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
}
\value{
An integer of length 1, or \code{integer(0)} if \code{x} is empty or has
only missing values. For vectors longer than \code{.Machine$integer.max}
the index is a double.
}
\seealso{
\code{\link{which.min}}, \code{\link{simd_min}}.
}
\examples{
simd_which_max(c(1, 3, 3))
simd_which_min(c(NA, 2, NaN, 1))
simd_which_max(numeric(0))
}
