\name{simd_sinh}
\alias{simd_sinh}
\alias{simd_cosh}
\alias{simd_tanh}
\alias{simd_asinh}
\alias{simd_acosh}
\alias{simd_atanh}
\title{Hyperbolic Functions}
\usage{
simd_sinh(x)

simd_cosh(x)

simd_tanh(x)

simd_asinh(x)

simd_acosh(x)

simd_atanh(x)
}
\arguments{
\item{x}{A double, integer or logical vector.}
}
\description{
\code{sinh(x)}, \code{cosh(x)}, \code{tanh(x)}, \code{asinh(x)},
\code{acosh(x)} and \code{atanh(x)} computed with SIMD instructions. They
correspond to \code{Sinh}, \code{Cosh}, \code{Tanh}, \code{Asinh},
\code{Acosh} and \code{Atanh} of .NET's \code{TensorPrimitives}.
}
\details{
\code{simd_sinh()} and \code{simd_cosh()} overflow to \code{Inf} where base
R's do (just above \eqn{|x| = 710.47}); \code{simd_asinh()} and
\code{simd_acosh()} are finite up to the largest double. SLEEF is less
accurate in those ranges (\eqn{|x| > 709} and \eqn{|x| > 10^{154}}), so
there the C library's functions are used instead. \code{simd_acosh(x)} for
\code{x < 1} and \code{simd_atanh(x)} for \eqn{|x| > 1} are \code{NaN} with a
warning; \code{simd_atanh(1)} is \code{Inf}.

\rsimdMath

\rsimdAttributes
}
\value{
A double vector of the length of \code{x}.
}
\seealso{
\code{\link{sinh}}, \code{\link{simd_exp}}, \link{rsimd_options}.
}
\examples{
simd_sinh(c(0, 1, 710))
simd_cosh(c(0, 1, 711))
simd_tanh(c(-20, 0.5))
simd_asinh(1e300)
simd_acosh(c(1, 10))
simd_atanh(c(0.5, 1))
}
