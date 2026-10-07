\name{simd_exp}
\alias{simd_exp}
\alias{simd_exp2}
\alias{simd_exp10}
\alias{simd_expm1}
\alias{simd_log}
\alias{simd_log2}
\alias{simd_log10}
\alias{simd_log1p}
\alias{simd_exp2m1}
\alias{simd_exp10m1}
\alias{simd_log2p1}
\alias{simd_log10p1}
\title{Exponentials and Logarithms}
\usage{
simd_exp(x)

simd_exp2(x)

simd_exp10(x)

simd_expm1(x)

simd_log(x, base = exp(1))

simd_log2(x)

simd_log10(x)

simd_log1p(x)

simd_exp2m1(x)

simd_exp10m1(x)

simd_log2p1(x)

simd_log10p1(x)
}
\arguments{
\item{x}{A double, integer or logical vector; complex for \code{simd_exp()},
\code{simd_log()}, \code{simd_log2()} and \code{simd_log10()}.}
\item{base}{A single number, the base of the logarithm; with complex
\code{x} or \code{base}, a vector following the broadcast rule.}
}
\description{
\code{exp(x)}, \code{2^x}, \code{10^x}, \code{expm1(x)} (\code{exp(x) - 1}),
\code{log(x, base)}, \code{log2(x)}, \code{log10(x)} and \code{log1p(x)}
(\code{log(1 + x)}), and \code{2^x - 1}, \code{10^x - 1},
\code{log2(1 + x)} and \code{log10(1 + x)} without their loss of accuracy
near 0, computed with SIMD instructions.
}
\details{
\code{simd_exp2(x)} and \code{simd_exp10(x)} correspond to base R's
\code{2^x} and \code{10^x}; the \code{"none"} implementation computes them
that way. They are exact where the result is a representable power
(\code{simd_exp10(0:22)}, \code{simd_exp2(-1074:1023)}), and so are
\code{simd_log2()} of the powers of 2 and \code{simd_log10()} of
\code{10^(0:22)}.

\code{simd_log(x, base)} takes a single \code{base} (base R's \code{log()}
also takes a vector of bases): base 2 and 10 are \code{simd_log2()} and
\code{simd_log10()}, any other base \code{simd_log(x) / log(base)}. Logarithms
of negative numbers are \code{NaN} with a warning; \code{simd_log(0)} is
\code{-Inf}.

Results overflow to \code{Inf} and underflow through the subnormal numbers
to zero where base R's do: \code{simd_exp(709.78)} is finite,
\code{simd_exp(709.79)} is \code{Inf}, \code{simd_exp(-745)} is
\code{5e-324} and \code{simd_exp(-746)} is 0.

\code{simd_exp2m1(x)}, \code{simd_exp10m1(x)}, \code{simd_log2p1(x)} and
\code{simd_log10p1(x)} are \code{2^x - 1}, \code{10^x - 1},
\code{log2(1 + x)} and \code{log10(1 + x)} (C23's \code{exp2m1},
\code{exp10m1}, \code{log2p1} and \code{log10p1}), without the
cancellation of those expressions near 0: \code{simd_exp2m1(1e-20)} is
\code{6.93e-21}, where \code{2^1e-20 - 1} is 0. They are computed from
\code{expm1} and \code{log1p}, with the scaling by \code{log(2)} or
\code{log(10)} carried in extra precision, and are within 2 ULP of the
exact value on every implementation (at most 1.6 ULP measured), with the
same results in both accuracy modes (SLEEF has no faster \code{expm1} or
\code{log1p}). They are exact where the base R expression is:
\code{simd_exp2m1()} at the whole numbers from -53 to 53,
\code{simd_exp10m1()} at 1 to 15, and \code{simd_log2p1()} and
\code{simd_log10p1()} wherever \code{1 + x} is exact and a power of 2 or
10 (\code{simd_log2p1(7)} is 3). Zeros keep their sign,
\code{simd_log2p1(-1)} is \code{-Inf}, and arguments below -1 give
\code{NaN} with a warning. They do not take complex input.

\rsimdMath

\rsimdAttributes
}
\section{Complex numbers}{
\code{simd_exp()}, \code{simd_log()} (also with a \code{base}),
\code{simd_log2()}, \code{simd_log10()}, \code{\link{simd_sqrt}()},
\code{\link{simd_pow}()}, the trigonometric and hyperbolic functions and
their inverses (\code{\link{simd_sin}()} to \code{\link{simd_atan}()},
\code{\link{simd_sinh}()} to \code{\link{simd_atanh}()}) and
\code{\link{simd_atan2}()} take complex input, as their base R counterparts
do. In the two-operand ones a double, integer or logical operand is
converted with \code{as.complex()}, and the operands follow the broadcast
rule (base R recycles any lengths). \code{simd_expm1()}, \code{simd_log1p()},
\code{simd_exp2()}, \code{simd_exp10()}, \code{\link{simd_cbrt}()},
\code{\link{simd_hypot}()}, \code{\link{simd_sinpi}()} and the other
functions base R does not define on complex numbers reject complex input.

The \code{"none"} implementation computes every element as base R does: the
C library's \code{csqrt}, \code{cexp}, \code{clog}, \code{csin} and so on,
and base R's own code for \code{tan}, for \code{asin}, \code{acos},
\code{atan}, \code{asinh} and \code{atanh} on their branch cuts, for powers,
for logarithms with a base and for \code{atan2}. It is identical to base R
on every platform, \code{NA} and \code{NaN} results and base R's warning
(such as \samp{NaNs produced in function "sqrt"}) included, except for
\code{simd_acosh()} (below).

The other implementations compute each element with both parts finite from
SLEEF's real functions (for example \eqn{e^x(\cos y + i \sin y)} for
\code{exp}, \eqn{\log|z| + i\,\mathrm{atan2}(y, x)} for \code{log} with
\eqn{\log|z|} as \code{log1p(x^2 + y^2 - 1) / 2} near \eqn{|z| = 1},
Kahan's formula for \code{tan} and \code{tanh}, and Hull, Fairgrieve and
Tang's algorithm for \code{asin} and \code{acos}). On the lines where base
R has its own branch-cut code for \code{asin}, \code{acos}, \code{asinh},
\code{atan} and \code{atanh} (such as \code{asin(2+0i)} or
\code{atan(2i)}) they compute that code's formula with base R's roundings,
only its \code{asin} and \code{log} being SLEEF's, so the results are
within 2 ULP of base R's there (bit for bit in most elements). The C
library's (and base R's) code computes every other element: those with an
infinite, \code{NaN} or \code{NA} part, a result with a \code{NaN} part,
and the regions the formulas do not handle accurately (parts beyond
\eqn{2^{500}} or near the subnormal range, \eqn{|x| > 708} in \code{exp},
a hyperbolic part beyond 20 in \code{tan} and \code{tanh}). So special
values, the choices on the branch cuts and the signs of zero parts are base
R's.
Each part of a result is within the bound below of the exact value, in
units in the last place (ULP), or within that many ULP of the result's
modulus where the part is much smaller than the modulus. The C library's
complex functions, and so base R and the \code{"none"} implementation, are
within about 5 ULP.
\tabular{lrr}{
  function \tab accurate \tab fast \cr
  \code{sqrt} \tab 2 \tab 3 \cr
  \code{exp} \tab 2.5 \tab 3 \cr
  \code{log} \tab 3 \tab 5.5 \cr
  \code{log} with a base, \code{log2}, \code{log10} \tab 4.5 \tab 6.5 \cr
  \code{sin}, \code{cos}, \code{sinh}, \code{cosh} \tab 3 \tab 4 \cr
  \code{tan} \tab 4.5 \tab 6 \cr
  \code{tanh} \tab 5 \tab 6 \cr
  \code{asin}, \code{acos}, \code{asinh}, \code{acosh} \tab 3.5 \tab 4.5 \cr
  \code{atan}, \code{atanh} \tab 3 \tab 3.5 \cr
  \code{atan2} \tab 3.5 \tab 4 \cr
  \code{pow}, per \eqn{1 + |w \log z|} \tab 3 \tab 5 \cr
}
In \link[=simd_math_accuracy]{fast mode} the formulas use SLEEF's 3.5-ULP
functions where SLEEF has them.

\strong{Powers.} \code{simd_pow(z, w)} for a whole-number real \code{w}
with \eqn{|w| \le 65536} multiplies by repeated squaring as base R does,
with base R's complex multiplication and division, so it is identical to
base R on every implementation. Other powers are \eqn{\exp(w \log z)}, as
in the C library's \code{cpow}, and their error grows with
\eqn{|w \log z|}; the sign of a zero part may differ from base R's.
\code{0^w} is base R's: \code{0^w} for a real \code{w}, \code{NaN} in both
parts otherwise.

\strong{acosh.} \code{simd_acosh()} is C's \code{cacosh}, the principal
value: its real part is never negative, and its imaginary part has the sign
of \code{Im(z)}, a zero's sign included (\code{simd_acosh(-2-0i)} is
\code{1.317-3.142i}). Base R computes \code{acosh(z)} as
\code{acos(z) * 1i}, which is the negative of the principal value wherever
its real part comes out negative, that is for the whole lower half-plane and
for real \eqn{z > 1}: \code{acosh(2-1i)} is \code{-1.469+0.507i} in base R
and \code{1.469-0.507i} in \pkg{rsimd}, \code{acosh(2+0i)} is
\code{-1.317+0i} in base R and \code{1.317+0i} in \pkg{rsimd}. Base R also
ignores the sign of a zero imaginary part there.
}

\value{
A double vector of the length of \code{x}; complex for complex input (with
\code{simd_log()}'s complex \code{base}, of the longer operand).
}
\seealso{
\code{\link{exp}}, \code{\link{log}}, \code{\link{simd_pow}},
\link{rsimd_options}.
}
\examples{
simd_exp(c(0, 1, -Inf, NA))
simd_exp2(c(10, 0.5))
simd_expm1(1e-10)
simd_log(c(1, exp(2), 0))
simd_log(c(8, 81), base = 3)
simd_log10(c(1000, 0.01))
simd_log1p(-0.5)
simd_exp(c(1i * pi, 1 + 1i))
simd_log(c(-1+0i, -1-0i, 1i))
simd_log(1i, base = 2i)
simd_exp2m1(c(1e-20, 3, -1))
2^1e-20 - 1 # 0: the difference is lost
simd_log10p1(c(1e-15, 99))
}
