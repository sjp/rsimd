\name{simd_and}
\alias{simd_and}
\alias{simd_or}
\alias{simd_xor}
\alias{simd_not}
\title{Logical Operators}
\usage{
simd_and(x, y)

simd_or(x, y)

simd_xor(x, y)

simd_not(x)
}
\arguments{
\item{x, y}{Logical, integer or double vectors, or two raw vectors.}
}
\description{
\code{x & y}, \code{x | y}, \code{xor(x, y)} and \code{!x} computed with
SIMD instructions, with base R's three-valued logic.
}
\details{
Each element is read as a logical value: non-zero is \code{TRUE}, zero
\code{FALSE}, and \code{NA} or \code{NaN} is \code{NA}. Then, as in base R,
\code{x & y} is \code{FALSE} if either is \code{FALSE}, else \code{NA} if
either is \code{NA}, else \code{TRUE}; \code{x | y} is \code{TRUE} if either
is \code{TRUE}, else \code{NA} if either is \code{NA}, else \code{FALSE};
\code{xor} is \code{NA} if either is \code{NA}; \code{!NA} is \code{NA}.
Integer and double operands give a logical result, like base R's operators.

Two raw vectors are combined bit by bit and give a raw vector, as base R
does (\code{simd_not()} of raw flips every bit). Mixing raw with another
type is an error with base R's message. For bitwise operations on integers
see \code{\link{simd_bit_and}}.

\rsimdLengths

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM double operands use the scalar \code{"none"} code. Complex and
\pkg{bit64} \code{integer64} operands are not supported.
}
\value{
A logical vector of the length of the longer operand, or a raw vector for
raw operands.
}
\seealso{
\code{\link{Logic}}, \code{\link{simd_bit_and}}, \code{\link{simd_eq}}.
}
\examples{
simd_and(c(TRUE, NA, NA), c(FALSE, TRUE, FALSE))
simd_or(c(TRUE, NA, NA), c(NA, TRUE, FALSE))
simd_xor(c(TRUE, FALSE), TRUE)
simd_not(c(0, 2.5, NaN))
simd_and(as.raw(0xF0), as.raw(0x3C))
}
