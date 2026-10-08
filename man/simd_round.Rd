\name{simd_round}
\alias{simd_round}
\alias{simd_floor}
\alias{simd_ceiling}
\alias{simd_trunc}
\title{Rounding}
\usage{
simd_floor(x)

simd_ceiling(x)

simd_trunc(x)

simd_round(x, digits = 0)
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{digits}{A single number: the number of decimal places to round to
(negative to round to tens, hundreds ...). A logical counts as a number
(\code{TRUE} is 1), as in base R; an integer64 is not accepted.}
}
\description{
\code{floor(x)}, \code{ceiling(x)}, \code{trunc(x)} and
\code{round(x, digits)} computed with SIMD instructions.
}
\details{
The results are doubles for every input, and identical to base R's:
exact for every double, \code{-0} kept (\code{simd_ceiling(-0.5)} and
\code{simd_round(-0.4)} are \code{-0}), \code{NA}, \code{NaN} and
infinities returned as they are. \code{simd_round(x)} rounds halves to even,
as base R does: \code{simd_round(c(0.5, 1.5, 2.5))} is \code{0 2 2}.

\code{simd_round(x, digits)} with \code{digits} other than 0 uses base R's
own algorithm for \code{round(x, digits)}, which is more careful than
scaling by a power of ten: of the two numbers with \code{digits} decimal
places either side of \code{x}, the closer one, with ties to the one
whose last digit is even. It is vectorised, with the same results as base
R bit for bit, when \code{digits} rounds to a whole number from -308 to
308; other values of \code{digits} go through base R's own code one
element at a time.
Unlike base R's \code{round()}, \code{digits} must be a single number
(\code{NA} gives \code{NA}).

\rsimdAttributes

All implementations have their own kernels for these functions. On
\code{"sse2"}, which has no rounding instruction, they add and subtract
\eqn{2^{52}}, which is exact. On 32-bit ARM they use the scalar
\code{"none"} code.
}
\value{
A double vector of the length of \code{x}.
}
\seealso{
\code{\link{round}}, \link{rsimd_options}.
}
\examples{
simd_floor(c(-1.5, 1.5))
simd_ceiling(c(-1.5, 1.5))
simd_trunc(c(-1.5, 1.5))
simd_round(c(0.5, 1.5, 2.5, -0.5))
simd_round(2.675, 2)
}
