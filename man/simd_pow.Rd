\name{simd_pow}
\alias{simd_pow}
\alias{simd_cbrt}
\alias{simd_hypot}
\alias{simd_rootn}
\alias{simd_rsqrt}
\title{Powers, Roots and Hypotenuses}
\usage{
simd_pow(x, y)

simd_cbrt(x)

simd_hypot(x, y)

simd_rootn(x, n)

simd_rsqrt(x)
}
\arguments{
\item{x, y}{Double, integer or logical vectors; for \code{simd_pow()}, also
complex.}
\item{n}{An integer or logical vector, or a double vector of whole numbers
in the integer range (\code{NA} allowed).}
}
\description{
\code{x^y}, the real cube root of \code{x}, \code{sqrt(x^2 + y^2)}, the
real \code{n}-th root of \code{x} and \code{1 / sqrt(x)}, computed with
SIMD instructions.
}
\details{
\code{simd_pow(x, y)} follows base R's \code{x^y}: \code{x^2} is
\code{x * x}, \code{x^3} and \code{x^4} are products for \code{|x| <= 11}
(as in R >= 4.6.0; older R calls \code{pow()} there, which can differ
by 1 ULP),
\code{1^y} and \code{x^0} are 1 for every \code{x} and \code{y}, even
\code{NA} and \code{NaN}, \code{0^y} is \code{0} for \code{y > 0} and
\code{Inf} for \code{y < 0} (also for \code{-0}), and infinite operands
follow R's rules, which differ from C's: \code{(-Inf)^0.5} and
\code{(-2)^Inf} are \code{NaN}. A single exponent of \code{0.5} or
\code{-1} gives \code{sqrt(x)} or \code{1 / x}, with base R's signed
zeros (\code{(-0)^0.5} is \code{0}, \code{(-0)^-1} is \code{Inf},
\code{(-Inf)^-1} is \code{0}): these are correctly rounded, so they can
differ by 1 ULP from base R, which calls \code{pow()}. Everything else is
the SLEEF (or C) power function, within 1 ULP; the \code{"none"}
implementation and those of arm64 (\code{"neon"}, \code{"sve"},
\code{"sve2"}) call the C library's, so they are identical to base R (see
\link{simd_math_accuracy}). On arm64 that makes \code{simd_pow()} about
as fast as \code{^} (between 0.85 and 1.25 times its speed, depending on
the data), and much faster only for the exponents 2, 0.5 and -1, which
use vector instructions on every implementation. Apart from \code{x^0} and \code{1^y}, an \code{NA}
operand gives \code{NA}, as in \code{\link{simd_add}()}. Like \code{^},
\code{simd_pow()} never warns: \code{simd_pow(-8, 1/3)} is a silent
\code{NaN}.

\code{simd_cbrt()} is C's \code{cbrt()}: the cube root of a negative number
is negative (\code{simd_cbrt(-8)} is \code{-2}), unlike
\code{(-8)^(1/3)}, which is \code{NaN}. Base R has no \code{cbrt()}; the
\code{"none"} implementation calls the C library's, which may be less
accurate than SLEEF's (up to 3 ULP with glibc on arm64, where
\code{cbrt(27)} is \code{3 + 4.4e-16}). On arm64 that accuracy costs
speed: in accurate mode \code{"neon"} is about 15\% slower than the C
library's \code{cbrt()} and 25\% slower than \code{abs(x)^(1/3)}; fast
mode (\link{simd_math_accuracy}) is about 1.4 times as fast as the
latter.

\code{simd_hypot()} is C's \code{hypot()}, which does not overflow or
underflow in the squares: \code{simd_hypot(1e200, 1e200)} is
\code{1.414214e200}. Base R has no \code{hypot()}. A missing operand gives a
missing result (\code{simd_hypot(Inf, NaN)} is \code{NaN}, unlike C's
\code{Inf}), as for base R's \code{atan2()}.

\code{simd_rootn(x, n)} is C23's \code{rootn()}: the real \code{n}-th root
of \code{x} for a whole number \code{n}, positive or negative. For an odd
\code{n} the root of a negative number is negative (\code{simd_rootn(-8, 3)}
is \code{-2}, where \code{(-8)^(1/3)} is \code{NaN}); for an even \code{n}
it is \code{NaN}. It is within 1 ULP of the exact root on every
implementation, and exact whenever the root is: \code{simd_rootn(1000, 3)}
is exactly \code{10} and \code{simd_rootn(2^60, 4)} exactly \code{2^15},
where \code{1000^(1/3)} is \code{9.999999999999998}: \code{x^(1/n)} rounds
\code{1/n} first, and can be more than 100 ULP off. The method is a first
approximation from \code{exp2(log2(m) / n)}, where \code{m} is \code{|x|}
scaled by a power of two, refined by one Newton step whose \code{n}-th power
is computed in double-double arithmetic, which makes the result correctly
rounded in nearly every case (measured: never more than 0.502 ULP from the
exact root); \code{simd_rootn(x, 2)} is \code{sqrt(x)} exactly,
\code{simd_rootn(x, 1)} is \code{x} and \code{simd_rootn(x, -1)} is \code{1
/ x}. The special values are C23's: \code{n = 0} gives \code{NaN}; a zero
gives a zero for a positive \code{n} and \code{Inf} for a negative one, and
an infinity the reverse, each with the sign of \code{x} for an odd \code{n}
and positive for an even one (so \code{simd_rootn(-0, 2)} is \code{0} and
\code{simd_rootn(-0, -3)} is \code{-Inf}); a negative \code{x}, including
\code{-Inf}, with an even \code{n} gives \code{NaN}. Every \code{NaN} that
does not come from a missing operand warns \code{"NaNs produced"}. An
\code{NA} (or \code{NaN}) in \code{n} gives \code{NA}; a double \code{n}
that is not a whole number in the integer range is an error. Base R has no
\code{rootn()}.

\code{simd_rsqrt(x)} is base R's \code{1 / sqrt(x)}, two correctly
rounded operations, so identical to it on every implementation, including
\code{-Inf} for \code{-0} and \code{"NaNs produced"} for a negative
\code{x}. \code{\link{simd_rsqrt_approx}()} is a faster approximation.

Neither \code{simd_rootn()} nor \code{simd_rsqrt()} depends on option
\code{rsimd.math_accuracy}.

The two-operand functions follow the broadcast rule:

\rsimdBroadcast

\code{simd_pow(x, y)} with a complex operand converts the other with
\code{as.complex()} and is base R's \code{x^y}: whole-number real powers up
to 65536 are identical to base R on every implementation, other powers are
\eqn{\exp(y \log x)}; see the Complex numbers section of
\code{\link{simd_exp}}.

\rsimdMath

\rsimdAttributes
}
\value{
A double vector, of the length of the longer operand for \code{simd_pow()},
\code{simd_hypot()} and \code{simd_rootn()}; complex for \code{simd_pow()}
with a complex operand.
}
\seealso{
\code{\link{Arithmetic}}, \code{\link{simd_sqrt}}, \code{\link{simd_exp}},
\code{\link{simd_recip_approx}}, \link{rsimd_options}.
}
\examples{
simd_pow(1 + 1i, c(2, -3, 0.5))
simd_pow(c(2, -2, 10), 3)
simd_pow(c(NA, NaN), 0)
simd_pow(-8, 1/3)
simd_cbrt(c(-8, 27, 2))
simd_hypot(3, 4)
simd_hypot(1e200, 1e200)
simd_rootn(c(27, -8, 2^60, 16), c(3L, 3L, 4L, -2L))
1000^(1/3) == 10
simd_rootn(1000, 3) == 10
simd_rsqrt(c(4, 0, -0))
}
