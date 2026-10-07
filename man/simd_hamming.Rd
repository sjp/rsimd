\name{simd_hamming}
\alias{simd_hamming}
\alias{simd_hamming_bits}
\title{Hamming Distances}
\usage{
simd_hamming(x, y, na.rm = FALSE)

simd_hamming_bits(x, y)
}
\arguments{
\item{x, y}{For \code{simd_hamming()}, double, integer, logical, raw,
complex or \pkg{bit64} \code{integer64} vectors. For
\code{simd_hamming_bits()}, two integer or logical vectors, two raw vectors
or two \code{integer64} vectors.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: should pairs with a missing
element be skipped?}
}
\description{
The number of positions where two vectors differ, and the number of bits in
which they differ, computed with SIMD instructions and without creating a
logical vector.
}
\details{
\code{simd_hamming(x, y)} is \code{sum(x != y)}, as a double. A pair in
which either element is \code{NA} or \code{NaN} makes the result
\code{NA}, unless \code{na.rm = TRUE}, which skips such pairs, as
\code{sum(x != y, na.rm = TRUE)} does. Doubles are compared as numbers, so
\code{0} and \code{-0} are equal. Complex numbers differ when either part
does, and are missing when either part is.

The operands are converted as for \code{\link{simd_eq}()} and base R's
\code{!=}: integers and logicals mixed with doubles compare as doubles, a
raw vector with a logical one as logical and with any other type as integer,
an \code{integer64} vector with integers as 64-bit integers, and with
doubles as doubles (with a warning, as in \code{\link{simd_add}()}). A
complex operand converts the other operand to complex.

\code{simd_hamming_bits(x, y)} is the total number of set bits of the
bitwise exclusive or of each pair: integers are compared as 32-bit patterns
(\code{NA} as its pattern \code{0x80000000}, without special treatment),
raw bytes as 8 bits and \code{integer64} values as 64 bits. Both operands
must have the same width: an integer cannot be compared with a raw or
\code{integer64} vector. The result is a double.

\rsimdLengths

Both functions return plain values, also for \code{\link{simd_vec}}
operands, which select their pinned implementation. Every implementation
has its own kernels for these functions, except that on 32-bit ARM double
and complex input uses the scalar \code{"none"} code.
}
\value{
A single double: the number of differing pairs (or \code{NA}), or the
number of differing bits.
}
\seealso{
\code{\link{simd_ne}}, \code{\link{simd_popcount_total}},
\code{\link{simd_bit_xor}}, \link{rsimd-integer64}.
}
\examples{
x <- c(1, 2, 3, NA)
y <- c(1, 5, 3, 4)
simd_hamming(x, y)
simd_hamming(x, y, na.rm = TRUE)
simd_hamming(c(1 + 1i, 2i), c(1 + 1i, 3i))
simd_hamming_bits(c(0L, 7L), c(1L, 0L))
simd_hamming_bits(as.raw(0xFF), as.raw(0x0F))
}
