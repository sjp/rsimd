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
\item{x, y}{Double, integer, logical, raw or \pkg{bit64} \code{integer64}
vectors; complex vectors for \code{simd_eq()} and \code{simd_ne()}.}
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
\pkg{bit64} \code{integer64} operands compare as 64-bit integers with
integer64, integer and logical operands; compared with a double they are
converted to double, with a warning (see \link{rsimd-integer64}).

\code{simd_eq()} and \code{simd_ne()} take complex operands: the other
operand is converted with \code{as.complex()} (integer64 through double,
with the same warning), both parts are compared, and the result is
\code{NA} wherever either operand has a missing part, as in base R. The
other comparisons give base R's error \code{"invalid comparison with
complex values"}.
}
\value{
A logical vector of the length of the longer operand.
}
\seealso{
\code{\link{Comparison}}, \code{\link{simd_is_na}}, \code{\link{simd_and}},
\link{rsimd-integer64}.
}
\examples{
simd_eq(c(1+1i, 2i, NA), 1+1i)
simd_eq(c(1, NA, NaN, -0), c(1, 1, 1, 0))
simd_lt(1:5, 3L)
simd_ge(c(TRUE, FALSE), 0.5)
simd_gt(as.raw(c(1, 200)), as.raw(100))
}
