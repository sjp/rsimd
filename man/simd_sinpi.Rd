\name{simd_sinpi}
\alias{simd_sinpi}
\alias{simd_cospi}
\alias{simd_tanpi}
\alias{simd_sincospi}
\title{Trigonometric Functions of Multiples of Pi}
\usage{
simd_sinpi(x)

simd_cospi(x)

simd_tanpi(x)

simd_sincospi(x)
}
\arguments{
\item{x}{A double, integer or logical vector.}
}
\description{
\code{sin(pi * x)}, \code{cos(pi * x)} and \code{tan(pi * x)} computed with
SIMD instructions, without the rounding error of \code{pi * x}, and the sine
and cosine together. They correspond to \code{SinPi}, \code{CosPi},
\code{TanPi} and \code{SinCosPi} of .NET's \code{TensorPrimitives}.
}
\details{
\code{x} is first reduced exactly, as base R does (modulo 2 for
\code{sinpi} and \code{cospi}, modulo 1 for \code{tanpi}), so the results
are exact at every integer and half-integer, however large:
\itemize{
\item \code{simd_sinpi(n)} is a zero with the sign of \code{n}
(\code{1 / simd_sinpi(-1)} is \code{-Inf}), and \code{simd_sinpi(n + 1/2)}
is 1 or -1.
\item \code{simd_cospi(n)} is 1 or -1, and \code{simd_cospi(n + 1/2)} is
\code{0}, never \code{-0}.
\item \code{simd_tanpi(n)} is \code{0} for every integer, including
\code{-0} and the negative ones (as base R's \code{tanpi()}),
\code{simd_tanpi(n + 1/4)} is 1 and \code{simd_tanpi(n - 1/4)} is -1, and
\code{simd_tanpi(n + 1/2)} is \code{NaN} with a warning.
}
Infinite arguments give \code{NaN} with a warning.

Elsewhere \code{simd_sinpi()} and \code{simd_cospi()} are within 0.5 ULP of
the exact value (SLEEF's \code{_u05} variants) and \code{simd_tanpi()},
their ratio, within about 2 ULP. The \code{"none"} implementation reduces
the same way and then calls the C library's \code{sin()}, \code{cos()} and
\code{tan()} on a small argument, within about 1 ULP. Base R calls the
platform's \code{sinpi()} and \code{cospi()} where there are any (and its
own, less accurate code elsewhere), and computes \code{tanpi(x)} as
\code{tan(pi * x)} after the reduction, which loses accuracy near the
poles at the half-integers; results may therefore differ from base R's by
more than 1 ULP there.

\code{simd_sincospi(x)} is \code{list(sin = simd_sinpi(x), cos =
simd_cospi(x))}, the same values bit for bit, computed in one pass (one call
of SLEEF's \code{sincospi} per vector), with at most one \code{"NaNs
produced"} warning; it is to \code{simd_sinpi()} and \code{simd_cospi()}
what \code{\link{simd_sincos}()} is to \code{simd_sin()} and
\code{simd_cos()}.

None of these functions takes complex input.

\rsimdMath

\rsimdAttributes
}
\value{
A double vector of the length of \code{x}; for \code{simd_sincospi()} a
list of two, \code{sin} and \code{cos}.
}
\seealso{
\code{\link{Trig}}, \code{\link{simd_sin}}, \link{rsimd_options}.
}
\examples{
simd_sinpi(c(0.5, 1, 1e300, 1/6))
1 / simd_sinpi(-1)
simd_cospi(c(0.5, 1, 2^53 + 2))
simd_tanpi(c(0.25, -0.25, 1, 0.3))
simd_sincospi(c(0, 1/6, 0.5))
}
