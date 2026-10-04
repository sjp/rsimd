\name{simd_pow}
\alias{simd_pow}
\alias{simd_cbrt}
\alias{simd_hypot}
\title{Powers, Cube Roots and Hypotenuses}
\usage{
simd_pow(x, y)

simd_cbrt(x)

simd_hypot(x, y)
}
\arguments{
\item{x, y}{Double, integer or logical vectors.}
}
\description{
\code{x^y}, the real cube root of \code{x} and \code{sqrt(x^2 + y^2)}
computed with SIMD instructions. They correspond to \code{Pow},
\code{Cbrt} and \code{Hypot} of .NET's \code{TensorPrimitives}.
}
\details{
\code{simd_pow(x, y)} follows base R's \code{x^y} exactly: \code{x^2} is
\code{x * x}, \code{x^3} and \code{x^4} are products for \code{|x| <= 11},
\code{1^y} and \code{x^0} are 1 for every \code{x} and \code{y}, even
\code{NA} and \code{NaN}, \code{0^y} is \code{0} for \code{y > 0} and
\code{Inf} for \code{y < 0} (also for \code{-0}), and infinite operands
follow R's rules, which differ from C's: \code{(-Inf)^0.5} and
\code{(-2)^Inf} are \code{NaN}. Everything else is the SLEEF (or C) power
function, within 1 ULP. Apart from \code{x^0} and \code{1^y}, an \code{NA}
operand gives \code{NA}, as in \code{\link{simd_add}()}. Like \code{^},
\code{simd_pow()} never warns: \code{simd_pow(-8, 1/3)} is a silent
\code{NaN}.

\code{simd_cbrt()} is C's \code{cbrt()}: the cube root of a negative number
is negative (\code{simd_cbrt(-8)} is \code{-2}), unlike
\code{(-8)^(1/3)}, which is \code{NaN}. Base R has no \code{cbrt()}; the
\code{"none"} implementation calls the C library's, which may be less
accurate than SLEEF's (up to 3 ULP with glibc on arm64, where
\code{cbrt(27)} is \code{3 + 4.4e-16}).

\code{simd_hypot()} is C's \code{hypot()}, which does not overflow or
underflow in the squares: \code{simd_hypot(1e200, 1e200)} is
\code{1.414214e200}. Base R has no \code{hypot()}. A missing operand gives a
missing result (\code{simd_hypot(Inf, NaN)} is \code{NaN}, unlike C's
\code{Inf}), as for base R's \code{atan2()}.

The two-operand functions follow the broadcast rule:

\rsimdBroadcast

\rsimdMath

\rsimdAttributes
}
\value{
A double vector, of the length of the longer operand for \code{simd_pow()}
and \code{simd_hypot()}.
}
\seealso{
\code{\link{Arithmetic}}, \code{\link{simd_sqrt}}, \code{\link{simd_exp}},
\link{rsimd_options}.
}
\examples{
simd_pow(c(2, -2, 10), 3)
simd_pow(c(NA, NaN), 0)
simd_pow(-8, 1/3)
simd_cbrt(c(-8, 27, 2))
simd_hypot(3, 4)
simd_hypot(1e200, 1e200)
}
