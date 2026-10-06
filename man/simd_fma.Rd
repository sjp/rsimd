\name{simd_fma}
\alias{simd_fma}
\alias{simd_mul_add}
\alias{simd_add_mul}
\alias{simd_lerp}
\title{Fused and Unfused Multiply-Add, Linear Interpolation}
\usage{
simd_fma(x, y, z, na_check = NULL)

simd_mul_add(x, y, z, na_check = NULL)

simd_add_mul(x, y, z, na_check = NULL)

simd_lerp(x, y, t, na_check = NULL)
}
\arguments{
\item{x, y, z}{Double, integer or logical vectors, or for
\code{simd_mul_add()} and \code{simd_add_mul()} \pkg{bit64}
\code{integer64} vectors.}
\item{t}{A double, integer or logical vector of interpolation weights.}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that the
operands have no missing values, so they are not looked for; see Details.}
}
\description{
\code{simd_fma(x, y, z)} is \code{x * y + z} computed with a single rounding
(a fused multiply-add), \code{simd_mul_add(x, y, z)} is \code{x * y + z} and
\code{simd_add_mul(x, y, z)} is \code{(x + y) * z}, each step rounded, and
\code{simd_lerp(x, y, t)} interpolates linearly from \code{x} (at
\code{t = 0}) to \code{y} (at \code{t = 1}). They correspond to
\code{FusedMultiplyAdd}, \code{MultiplyAdd}, \code{AddMultiply} and
\code{Lerp} of .NET's \code{TensorPrimitives}.
}
\details{
\rsimdBroadcast
Every operand may be broadcast, in any position.

\code{simd_fma()} and \code{simd_lerp()} always return doubles. They round
once on every implementation, so their results are identical everywhere: the
implementations without a fused multiply-add instruction (\code{"none"} and
\code{"sse2"}, and 32-bit ARM) call the C library's \code{fma()}, which is
exact but can be much slower. \code{simd_fma(0.1, 10, -1)} is
\eqn{2^{-54}} where \code{0.1 * 10 - 1} is 0.

\code{simd_lerp(x, y, t)} is \code{fma(t, y, (1 - t) * x)}, the formula
\code{x * (1 - t) + y * t} of \code{TensorPrimitives} with the last step
fused. It returns \code{x} exactly at \code{t = 0} and \code{y} exactly at
\code{t = 1} (for finite operands), which \code{x + t * (y - x)} does not;
it is not guaranteed to be monotonic in \code{t}.

\code{simd_mul_add()} and \code{simd_add_mul()} of integer and logical
operands return integers and follow base R's \code{x * y + z} and
\code{(x + y) * z} exactly: an intermediate result out of range is
\code{NA}, with the overflow warning, even if the final value would fit.

\rsimdSemantics

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double operations use the scalar \code{"none"} code.
\code{simd_mul_add()} and \code{simd_add_mul()} of integer64 operands (with
integer64, integer or logical ones) follow the integer rules in 64 bits and
return integer64 (see \link{rsimd-integer64}); \code{simd_fma()} and \code{simd_lerp()}
do not take integer64. None of them takes complex operands.
}
\value{
A vector of the length of the longest operand: integer for
\code{simd_mul_add()} and \code{simd_add_mul()} when every operand is
integer or logical, otherwise double.
}
\seealso{
\code{\link{simd_add}}, \code{\link{simd_mul}},
\code{\link{simd_mul_add_approx}} (fused where that is fast),
\link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_fma(0.1, 10, -1)
simd_mul_add(0.1, 10, -1)
simd_add_mul(1:3, 1L, 2L)
simd_lerp(0, 10, c(0, 0.25, 1))
}
