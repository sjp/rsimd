\name{simd_exp}
\alias{simd_exp}
\alias{simd_exp2}
\alias{simd_exp10}
\alias{simd_expm1}
\alias{simd_log}
\alias{simd_log2}
\alias{simd_log10}
\alias{simd_log1p}
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
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{base}{A single number, the base of the logarithm.}
}
\description{
\code{exp(x)}, \code{2^x}, \code{10^x}, \code{expm1(x)} (\code{exp(x) - 1}),
\code{log(x, base)}, \code{log2(x)}, \code{log10(x)} and \code{log1p(x)}
(\code{log(1 + x)}) computed with SIMD instructions. They correspond to
\code{Exp}, \code{Exp2}, \code{Exp10}, \code{ExpM1}, \code{Log},
\code{Log2}, \code{Log10} and \code{LogP1} of .NET's
\code{TensorPrimitives}.
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

\rsimdMath

\rsimdAttributes
}
\value{
A double vector of the length of \code{x}.
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
}
