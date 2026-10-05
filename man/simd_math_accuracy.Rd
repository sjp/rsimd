\name{simd_math_accuracy}
\alias{simd_math_accuracy}
\title{Accuracy of the Elementary Functions}
\usage{
simd_math_accuracy(mode)
}
\arguments{
\item{mode}{\code{"accurate"} or \code{"fast"}.}
}
\description{
Gets or sets which variants of the bundled SLEEF library the SIMD
implementations of the elementary functions use: the accurate ones (the
default), or faster ones that are less accurate in the last bits.
}
\details{
In \code{"accurate"} mode the SIMD implementations use SLEEF's 1-ULP
(\code{_u10}) functions, and its 0.5-ULP (\code{_u05}) ones for
\code{sinpi}, \code{cospi} and \code{hypot}. In \code{"fast"} mode they use
SLEEF's 3.5-ULP (\code{_u35}) functions wherever SLEEF has one that is
faster, which is typically 1.5 to 3 times as fast (6 times for
\code{hypot}). The bounds are in units in the last place (ULP) of the exact
result:
\tabular{lll}{
\strong{Function} \tab \strong{accurate} \tab \strong{fast} \cr
\code{sin}, \code{cos}, \code{tan}, \code{simd_sincos} \tab 1 \tab 3.5 \cr
\code{asin}, \code{acos}, \code{atan}, \code{atan2} \tab 1 \tab 3.5 \cr
\code{sinh}, \code{cosh}, \code{tanh} \tab 1 \tab 3.5 \cr
\code{log}, \code{log2}, \code{cbrt} \tab 1 \tab 3.5 \cr
\code{log(x, base)} for a base other than 2 and 10 \tab \code{log}'s, then a division \tab \code{log}'s, then a division \cr
\code{sinpi}, \code{cospi} \tab 0.5 \tab 3.5 \cr
\code{tanpi} \tab about 2 \tab about 4 \cr
\code{hypot} \tab 0.5 \tab 3.5 \cr
\code{exp}, \code{exp2}, \code{exp10}, \code{expm1} \tab 1 \tab 1 \cr
\code{log10}, \code{log1p}, \code{pow} \tab 1 \tab 1 \cr
\code{asinh}, \code{acosh}, \code{atanh} \tab 1 \tab 1 \cr
\code{simd_sigmoid}, \code{simd_softmax}, \code{simd_log_softmax} \tab as \code{exp} \tab as \code{exp} \cr
}
The functions whose bound is the same in both modes give identical results
in both: SLEEF has no 3.5-ULP variant of them, or (for \code{exp2} and
\code{exp10}) none that is faster. Only the last bits of finite results
change between the modes. Missing values (\code{NA} and \code{NaN} kept
apart), infinities, the signs of zero results, the exact cases (such as
\code{simd_sinpi(1)} or \code{simd_log2(8)}), the range reductions and the
\code{"NaNs produced"} warnings are the same in both.

The \code{"none"} implementation always uses the C math library, so the
mode makes no difference to it, and neither does it to an implementation
built without SLEEF (see \code{\link{simd_available}}). The mode applies to
every call, whichever implementation a \code{\link{simd_vec}} is pinned to.

The mode is stored in the option \code{rsimd.math_accuracy} (see
\link{rsimd_options}). Setting that option directly works too; an invalid
value is reported by the next function that uses it.
}
\value{
Without \code{mode}, the current mode. With \code{mode}, the previous mode,
invisibly; an invalid \code{mode} is an error and leaves the mode unchanged.
}
\seealso{
\link{simd_exp}, \link{simd_sin}, \link{simd_sinpi}, \link{simd_sinh} and
\link{simd_pow} for the functions; \link{rsimd_options} for all options.
}
\examples{
simd_math_accuracy()
x <- seq(-3, 3, length.out = 7)
old <- simd_math_accuracy("fast")
simd_sin(x)
simd_math_accuracy(old)
simd_sin(x)
}
