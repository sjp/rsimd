\name{simd_recip_approx}
\alias{simd_recip_approx}
\alias{simd_rsqrt_approx}
\alias{simd_mul_add_approx}
\title{Approximate Reciprocals, Reciprocal Square Roots and Multiply-Add}
\usage{
simd_recip_approx(x)

simd_rsqrt_approx(x)

simd_mul_add_approx(x, y, z, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x}{A double, integer or logical vector; for
\code{simd_mul_add_approx()} also a \pkg{bit64} \code{integer64} vector.}
\item{y, z}{As \code{x} for \code{simd_mul_add_approx()}.}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that the
operands have no missing values; see \code{\link{simd_mul_add}}.}
}
\description{
Faster, less accurate versions of \code{1 / x}, \code{1 / sqrt(x)} and
\code{x * y + z}, using whatever is fastest on each implementation. They
correspond to \code{ReciprocalEstimate}, \code{ReciprocalSqrtEstimate} and
\code{MultiplyAddEstimate} of .NET's \code{TensorPrimitives}. Results are
not reproducible across implementations: use \code{\link{simd_recip}},
\code{\link{simd_rsqrt}} and \code{\link{simd_fma}} or
\code{\link{simd_mul_add}} when they must be.
}
\details{
\code{simd_recip_approx(x)} and \code{simd_rsqrt_approx(x)} are within a
relative error of \eqn{2^{-22}}{2^-22} (about \code{2.4e-7}, roughly single
precision) of \code{1 / x} and \code{1 / sqrt(x)} for every finite double,
on every implementation. Each implementation starts from a hardware
estimate and refines it with Newton steps, or computes the exact value
where that is faster:
\tabular{llll}{
\strong{Implementation} \tab \code{simd_recip_approx} \tab \code{simd_rsqrt_approx} \tab \code{simd_mul_add_approx} \cr
\code{"none"} \tab exact \tab exact (two roundings) \tab rounded twice \cr
\code{"sse2"} \tab RCPPS + 1 step, \eqn{2^{-22.8}}{2^-22.8} \tab RSQRTPS + 1 step, \eqn{2^{-22.2}}{2^-22.2} \tab rounded twice \cr
\code{"avx2"} \tab RCPPS + 1 step, \eqn{2^{-22.8}}{2^-22.8} \tab RSQRTPS + 1 step, \eqn{2^{-22.2}}{2^-22.2} \tab fused \cr
\code{"avx512"} \tab VRCP14PD + 1 step, \eqn{2^{-28}}{2^-28} \tab VRSQRT14PD + 1 step, \eqn{2^{-27.4}}{2^-27.4} \tab fused \cr
\code{"neon"} (64-bit) \tab exact \tab FRSQRTE + 2 steps, \eqn{2^{-31.3}}{2^-31.3} \tab fused \cr
\code{"sve"}, \code{"sve2"} \tab FRECPE + 2 steps, \eqn{2^{-33.8}}{2^-33.8} \tab FRSQRTE + 2 steps, \eqn{2^{-31.3}}{2^-31.3} \tab fused \cr
}
The errors are the largest relative errors. For SSE2 and AVX2 they follow
from the largest error Intel documents for its float estimates
(\eqn{1.5 \times 2^{-12}}{1.5 * 2^-12}) after one step; the others were
measured. SSE2
and AVX2 have estimate instructions for floats only, so they estimate the
reciprocal of the significand of \code{x} and put the exponent back
exactly, which works over the whole range of doubles. On 64-bit ARM an
exact division was measured to be faster than the estimate and the two
steps it needs (on the Apple M-series), so \code{simd_recip_approx} divides
there. 32-bit ARM, and an implementation built without SLEEF (see
\code{\link{simd_available}}), use the next lower implementation's
functions, ultimately \code{"none"}'s.

The values where an estimate is least reliable are computed exactly on
every implementation: zeros (\code{1 / 0} is \code{Inf}, \code{1 / -0} is
\code{-Inf}), infinities, subnormal numbers, \code{|x| >= 2^1022} (whose
reciprocal is subnormal) and, for \code{simd_rsqrt_approx()}, negative
numbers, which give \code{NaN} with the warning \code{"NaNs produced"}, as
\code{1 / sqrt(x)} does. \code{NA} and \code{NaN} are returned as they are.
Integer and logical input is read as double.

\code{simd_mul_add_approx(x, y, z)} is \code{x * y + z} rounded once (a
fused multiply-add) on the implementations with a fused multiply-add
instruction, and rounded twice, as \code{simd_mul_add()}, on the others
(\code{"none"}, \code{"sse2"} and 32-bit ARM), where a fused multiply-add
would be slow; either way it is at least as accurate as
\code{simd_mul_add()}. For
integer, logical and integer64 operands there is nothing to approximate: it
is \code{simd_mul_add()} exactly, with its integer result and overflow
checks. Every operand may be broadcast, and missing values follow
\code{simd_mul_add()}.

None of these functions depends on option \code{rsimd.math_accuracy}.

\rsimdAttributes
}
\value{
A double vector, of the length of the longest operand; for
\code{simd_mul_add_approx()} of integer or integer64 operands, as
\code{\link{simd_mul_add}()}.
}
\seealso{
\code{\link{simd_recip}}, \code{\link{simd_rsqrt}}, \code{\link{simd_fma}},
\code{\link{simd_mul_add}}, \link{rsimd_options}.
}
\examples{
x <- c(3, 0.1, 1e300, 7e-300)
(simd_recip_approx(x) - 1 / x) * x
(simd_rsqrt_approx(x) - 1 / sqrt(x)) * sqrt(x)
simd_recip_approx(c(0, -0, Inf, NA))
simd_mul_add_approx(0.1, 10, -1)
simd_mul_add_approx(2L, 3L, 1L)
}
