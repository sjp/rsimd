\name{simd_eq}
\alias{simd_eq}
\alias{simd_ne}
\alias{simd_lt}
\alias{simd_le}
\alias{simd_gt}
\alias{simd_ge}
\title{Elementwise Comparisons}
\usage{
simd_eq(x, y)

simd_ne(x, y)

simd_lt(x, y)

simd_le(x, y)

simd_gt(x, y)

simd_ge(x, y)
}
\arguments{
\item{x, y}{Double, integer, logical or raw vectors.}
}
\description{
\code{x == y}, \code{x != y}, \code{x < y}, \code{x <= y}, \code{x > y} and
\code{x >= y} computed with SIMD instructions. They correspond to
\code{Equals}, \code{LessThan}, \code{GreaterThan} and their variants in
.NET's \code{TensorPrimitives}, but return \code{NA} where base R does.
}
\details{
The result is the same as base R's operator: \code{NA} wherever either
operand is \code{NA} or \code{NaN}, and the IEEE comparison elsewhere, so
\code{-0 == 0} is \code{TRUE}. Integer and logical operands compare as
integers; an operand that is double makes the comparison double (integers
are converted exactly). Raw vectors compare as unsigned bytes with each
other; compared with another type they are converted as base R does: to
logical (non-zero is \code{TRUE}) when the other operand is logical, else
to integer.

\rsimdLengths

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double comparisons use the scalar \code{"none"} code.
Complex and \pkg{bit64} \code{integer64} operands are not supported.
}
\value{
A logical vector of the length of the longer operand.
}
\seealso{
\code{\link{Comparison}}, \code{\link{simd_is_na}}, \code{\link{simd_and}}.
}
\examples{
simd_eq(c(1, NA, NaN, -0), c(1, 1, 1, 0))
simd_lt(1:5, 3L)
simd_ge(c(TRUE, FALSE), 0.5)
simd_gt(as.raw(c(1, 200)), as.raw(100))
}
