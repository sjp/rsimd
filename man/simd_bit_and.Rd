\name{simd_bit_and}
\alias{simd_bit_and}
\alias{simd_bit_or}
\alias{simd_bit_xor}
\alias{simd_bit_not}
\alias{simd_shl}
\alias{simd_shr}
\alias{simd_sar}
\alias{simd_rotl}
\alias{simd_rotr}
\title{Bitwise Operations, Shifts and Rotates}
\usage{
simd_bit_and(x, y, na_check = NULL)

simd_bit_or(x, y, na_check = NULL)

simd_bit_xor(x, y, na_check = NULL)

simd_bit_not(x, na_check = NULL)

simd_shl(x, n, na_check = NULL)

simd_shr(x, n, na_check = NULL)

simd_sar(x, n, na_check = NULL)

simd_rotl(x, n, na_check = NULL)

simd_rotr(x, n, na_check = NULL)
}
\arguments{
\item{x, y}{Integer, logical or \pkg{bit64} \code{integer64} vectors, or
raw vectors (not for \code{simd_sar}). \code{x} and \code{y} must both be
raw or both not.}
\item{n}{The number of bit positions, a single number (truncated toward
zero); see Details for its range.}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that the
integer operands have no \code{NA}, so they are not looked for; an
\code{NA} is then treated as its bit pattern \code{0x80000000}.}
}
\description{
Bit-by-bit and, or, exclusive or and complement, shifts and rotates of
integers (as 32-bit patterns) and of raw bytes, computed with SIMD
instructions. They correspond to base R's \code{\link{bitwAnd}} family and
\code{\link{rawShift}}.
}
\details{
\strong{Integers.} The operations act on the 32-bit two's complement
pattern of each element, and the result is read back as an integer, as
\code{bitwAnd()}, \code{bitwNot()}, \code{bitwShiftL()} and
\code{bitwShiftR()} do: \code{simd_shl()} and \code{simd_shr()} shift
logically (zeros come in), so \code{simd_shr(-1L, 1L)} is
\code{2147483647L}. \code{simd_sar()} shifts arithmetically (copies of the
sign bit come in), so \code{simd_sar(-8L, 1L)} is \code{-4L}; base R has no
equivalent. \code{simd_rotl()} and \code{simd_rotr()} rotate the 32 bits.
Logical operands count as integers.

An \code{NA} operand gives \code{NA}. Beware that the bit pattern
\code{0x80000000} is \code{NA_integer_}, so any result with that pattern
reads as \code{NA}: \code{simd_shl(1L, 31L)}, \code{simd_rotl(1L, 31L)}
and \code{simd_bit_not(.Machine$integer.max)} are all \code{NA}, as
\code{bitwShiftL(1L, 31L)} and \code{bitwNot(.Machine$integer.max)} are in
base R.

\strong{Raw.} Raw vectors are combined byte by byte and give raw vectors.
Shifts lose the bits shifted out, as \code{rawShift()} does
(\code{simd_shr(x, n)} is \code{rawShift(x, -n)}), and rotates are 8-bit
rotates. \code{simd_sar()} does not take raw vectors.

\strong{Counts.} For integers, a shift count \code{n} outside 0 to 31 (or
\code{NA}) gives \code{NA} for every element, like \code{bitwShiftL()};
for raw, \code{n} must be in 0 to 8 (8 gives \code{00}), otherwise it is an
error with \code{rawShift()}'s message. Rotate counts are taken modulo the
width (32, or 8 for raw), so a negative \code{n} rotates the other way;
\code{n} must not be \code{NA}. Only a single count is supported.

\rsimdLengths

\rsimdAttributes

All implementations have their own kernels for these functions. Double
operands are rejected (convert them with \code{\link{simd_as_integer}()}
first), and complex operands are not supported. \pkg{bit64}
\code{integer64} operands (with integer64, integer or logical operands)
give integer64 results; shift counts are then 0 to 63 and rotate counts are
taken modulo 64 (see \link{rsimd-integer64}).
}
\value{
An integer vector, or a raw vector for raw operands, of the length of the
longer operand (of \code{x} for the functions of one vector).
}
\seealso{
\code{\link{bitwAnd}}, \code{\link{rawShift}}, \code{\link{simd_popcount}},
\code{\link{simd_and}}, \link{rsimd-integer64}.
}
\examples{
simd_bit_and(5L, 3L)
simd_bit_xor(c(1L, 2L, 3L), 1L)
simd_bit_not(0L)
simd_shl(3L, 30L)
simd_shr(-1L, 1L)
simd_sar(-8L, 1L)
simd_rotl(1L, 31L) # 0x80000000 is NA
simd_rotr(as.raw(0x81), 1)
simd_shl(as.raw(0x81), 1)
}
