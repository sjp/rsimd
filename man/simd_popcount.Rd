\name{simd_popcount}
\alias{simd_popcount}
\alias{simd_popcount_total}
\alias{simd_lzcnt}
\alias{simd_tzcnt}
\title{Bit Counts}
\usage{
simd_popcount(x, na_check = NULL)

simd_popcount_total(x, ..., na.rm = FALSE, na_check = NULL)

simd_lzcnt(x, na_check = NULL)

simd_tzcnt(x, na_check = NULL)
}
\arguments{
\item{x}{An integer, logical, raw or \pkg{bit64} \code{integer64} vector.}
\item{\dots}{Must be empty. \code{simd_popcount_total()} takes one vector,
and \code{na.rm} and \code{na_check} must be named, so that an extra vector
is not read as \code{na.rm}.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: should \code{NA} elements be
skipped?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} has no \code{NA}, so it is not looked for; an \code{NA} is then
counted as its bit pattern \code{0x80000000}. \code{na.rm = TRUE} always
looks for \code{NA}.}
}
\description{
The number of set bits (population count), of leading zero bits and of
trailing zero bits of each element, and the total number of set bits of a
vector, computed with SIMD instructions. Base R has no equivalent.
}
\details{
Integers are counted as 32-bit two's complement patterns, so
\code{simd_popcount(-1L)} is 32, and \code{simd_lzcnt(0L)} and
\code{simd_tzcnt(0L)} are 32. Raw bytes are counted as 8 bits, so the counts
are between 0 and 8. Logical elements count as integers. An \code{NA}
element gives \code{NA}.

\code{simd_popcount_total()} is the sum of \code{simd_popcount(x)}, as a
double (it can exceed the integer range). It is \code{NA} if \code{x}
contains an \code{NA}, unless \code{na.rm = TRUE}.

\rsimdAttributes

All implementations have their own kernels for these functions; on x86 the
counts use shift-and-add sequences, as the \code{"avx512"} tier does not
include the AVX-512 extensions with lane-wise bit counts. Double input is
rejected. \pkg{bit64} \code{integer64} input counts the bits of the 64-bit
patterns (0 to 64; see \link{rsimd-integer64}).
}
\value{
An integer vector of the length of \code{x}, or for
\code{simd_popcount_total()} a single double.
}
\seealso{
\code{\link{simd_bit_and}}, \code{\link{bitwAnd}}, \link{rsimd-integer64}.
}
\examples{
simd_popcount(c(0L, 7L, -1L, NA))
simd_popcount(as.raw(c(0x0F, 0xFF)))
simd_lzcnt(c(0L, 1L, 65536L, -1L))
simd_tzcnt(c(0L, 1L, 8L, 65536L))
simd_popcount_total(c(1L, 3L, NA), na.rm = TRUE)
}
