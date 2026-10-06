\name{simd_sum_sq}
\alias{simd_sum_sq}
\alias{simd_sum_abs}
\alias{simd_norm}
\title{Sums of Squares, Sums of Absolute Values and Euclidean Norms}
\usage{
simd_sum_sq(x, na.rm = FALSE, na_check = NULL)

simd_sum_abs(x, na.rm = FALSE, na_check = NULL)

simd_norm(x, na.rm = FALSE, na_check = NULL)
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
\code{simd_sum_sq(x)} is \code{sum(x^2)}, \code{simd_sum_abs(x)} is
\code{sum(abs(x))} and \code{simd_norm(x)} is the Euclidean (L2) norm
\code{sqrt(sum(x^2))}, computed in one pass with SIMD instructions and
without forming \code{x^2} or \code{abs(x)}. They correspond to
\code{SumOfSquares}, \code{SumOfMagnitudes} and \code{Norm} of .NET's
\code{TensorPrimitives}.
}
\details{
Integer and logical elements are converted to double as they are read,
except in \code{simd_sum_abs()}, which sums them exactly as integers.

\code{simd_norm()} takes the square root of the sum of squares, so it
overflows to \code{Inf} (or underflows to 0) when the squares do, even
where the norm itself is representable: \code{simd_norm(c(1e200, 1e200))}
is \code{Inf}. Rescale such data first.

\rsimdSemantics

\rsimdPrecision

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double versions use the scalar \code{"none"} code.

\rsimdAttributes
}
\value{
\code{simd_sum_sq()} and \code{simd_norm()}: a double; 0 for an empty
vector.

\code{simd_sum_abs()}: a double for double \code{x}. For integer and
logical \code{x}, an integer, or a double if the sum is outside the integer
range, as \code{sum(abs(x))}. The sum of an empty vector is \code{0} (or
\code{0L}).
}
\seealso{
\code{\link{simd_dot}}, \code{\link{simd_sum}}, \code{\link{simd_precision}},
\link{rsimd_options}.
}
\examples{
simd_sum_sq(1:4)
simd_sum_abs(c(-2L, 3L, NA), na.rm = TRUE)
simd_norm(c(3, 4))
}
