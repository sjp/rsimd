\name{simd_sin}
\alias{simd_sin}
\alias{simd_cos}
\alias{simd_tan}
\alias{simd_sincos}
\alias{simd_asin}
\alias{simd_acos}
\alias{simd_atan}
\alias{simd_atan2}
\title{Trigonometric Functions}
\usage{
simd_sin(x)

simd_cos(x)

simd_tan(x)

simd_sincos(x)

simd_asin(x)

simd_acos(x)

simd_atan(x)

simd_atan2(y, x)
}
\arguments{
\item{x, y}{Double, integer, logical or complex vectors; angles in radians.
\code{simd_sincos()} takes no complex input.}
}
\description{
\code{sin(x)}, \code{cos(x)}, \code{tan(x)}, both \code{sin(x)} and
\code{cos(x)} in one pass, \code{asin(x)}, \code{acos(x)}, \code{atan(x)}
and \code{atan2(y, x)} computed with SIMD instructions. They correspond to
\code{Sin}, \code{Cos}, \code{Tan}, \code{SinCos}, \code{Asin},
\code{Acos}, \code{Atan} and \code{Atan2} of .NET's
\code{TensorPrimitives}.
}
\details{
Large arguments are reduced exactly (SLEEF uses Payne-Hanek reduction), so
\code{simd_sin(1e22)} is \code{-0.8522008}, as in base R.
\code{simd_sin(Inf)} is \code{NaN} with a warning, as are \code{simd_asin()}
and \code{simd_acos()} outside \eqn{[-1, 1]}.

\code{simd_sincos(x)} returns \code{list(sin = , cos = )} and warns at most
once. SLEEF computes both together, so its results may differ from
\code{simd_sin(x)} and \code{simd_cos(x)} in the last bit.

\code{simd_atan2(y, x)} is the angle of the point \code{(x, y)}, with C's
signed-zero rules (\code{simd_atan2(0, -0)} is \code{pi},
\code{simd_atan2(-0, -0)} is \code{-pi}). Its operands follow the broadcast
rule:

\rsimdBroadcast

See \code{\link{simd_sinpi}()} for \code{sin(pi * x)} and friends, which are
exact at multiples of 1/2.

\strong{Complex numbers.} These functions take complex input, as base R's
do: see the Complex numbers section of \code{\link{simd_exp}} for the
\code{"none"} implementation (base R exactly), the accuracy of the others and
the special values. \code{simd_atan2(y, x)} of complex operands (one may be
double, integer or logical) is base R's: \code{atan(y / x)}, plus \eqn{\pi}
for \code{Re(x) < 0} (minus \eqn{2\pi} above \eqn{\pi}); \code{atan2(0i, 0i)}
is \code{NA} with a warning.

\rsimdMath

\rsimdAttributes
}
\value{
A double vector of the length of \code{x} (for \code{simd_atan2()}, of the
longer operand), complex for complex input; for \code{simd_sincos()}, a list of two such vectors named
\code{sin} and \code{cos}.
}
\seealso{
\code{\link{Trig}}, \code{\link{simd_sinpi}}, \code{\link{simd_sinh}},
\link{rsimd_options}.
}
\examples{
simd_sin(c(0, pi / 2, 1e22))
simd_cos(pi)
simd_sincos(c(0, pi / 4))
simd_asin(c(-1, 0.5, 2))
simd_atan2(c(1, -1, 0), c(0, 0, -1))
simd_asin(2 + 0i)
simd_tan(1 + 30i)
}
