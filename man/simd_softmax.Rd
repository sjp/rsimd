\name{simd_softmax}
\alias{simd_softmax}
\alias{simd_log_softmax}
\alias{simd_sigmoid}
\title{Sigmoid, Softmax and Log-Softmax}
\usage{
simd_sigmoid(x)

simd_softmax(x)

simd_log_softmax(x)
}
\arguments{
\item{x}{A double, integer or logical vector.}
}
\description{
The logistic sigmoid \eqn{1 / (1 + e^{-x})} of each element, and the
softmax \eqn{e^{x_i} / \sum_j e^{x_j}} and its logarithm over the whole of
\code{x}, computed with SIMD instructions.
}
\details{
\strong{Sigmoid.} \code{simd_sigmoid(x)} is computed as
\eqn{1 / (1 + e)} for \eqn{x \ge 0} and \eqn{e / (1 + e)} for \eqn{x < 0},
with \eqn{e = \exp(-|x|)}{e = exp(-|x|)}, which is accurate in both tails and never
overflows. It equals \code{\link{plogis}(x)} within a few units in the last
place, except below \eqn{x = -709.78}, where \code{plogis()}
(\eqn{1 / (1 + \exp(-x))}{1 / (1 + exp(-x))}) underflows to 0 while \code{simd_sigmoid()}
gives the subnormal value, about \eqn{e^x}, down to about \eqn{x = -745}.
\code{simd_sigmoid(0)} is exactly 0.5, \code{simd_sigmoid(-Inf)} is 0 and
\code{simd_sigmoid(Inf)} is 1. It follows the missing-value rules below
and never warns.

\strong{Softmax.} \code{simd_softmax(x)} and \code{simd_log_softmax(x)}
normalise over the whole vector: a matrix is treated as one vector (there
is no per-row or per-column version). They use the stable forms
\deqn{\mathrm{softmax}(x)_i = \frac{e^{x_i - m}}{\sum_j e^{x_j - m}},
\qquad \mathrm{log\_softmax}(x)_i = (x_i - m) - \log \sum_j e^{x_j - m}}{
softmax(x)[i] = exp(x[i] - m) / sum(exp(x - m)),
log_softmax(x)[i] = (x[i] - m) - log(sum(exp(x - m)))}
with \eqn{m = \max(x)}{m = max(x)}, the same formulas as base R's
\code{exp(x - max(x)) / sum(exp(x - max(x)))} and
\code{x - max(x) - log(sum(exp(x - max(x))))}. They do not overflow, and
\code{simd_log_softmax()} stays finite where \code{log(simd_softmax(x))}
underflows to \code{-Inf}: \code{simd_log_softmax(c(0, 800, 1600))} is
\code{c(-1600, -800, 0)}. The sum of the exponentials is accumulated in the
precision mode (see below), and the division is a true division, not a
multiplication by the reciprocal. Each element is then within a few units
in the last place of base R's formula plus the relative error of the sum:
at most about \eqn{8n} units in \code{"fast"} mode (much less in practice),
\eqn{128 + \log_2 n}{128 + log2(n)} in \code{"pairwise"} mode and a few in
\code{"compensated"} mode.

The result has the length of \code{x}; a zero-length \code{x} gives
\code{numeric(0)} and a single element gives 1 (\code{simd_softmax}) or 0
(\code{simd_log_softmax}).

\strong{Missing and infinite values.} If \code{x} has an \code{NA}, every
element of the result is \code{NA}; otherwise if it has a \code{NaN},
every element is \code{NaN}. \code{-Inf} elements give 0 (\code{-Inf} for
\code{simd_log_softmax}) next to finite ones. If the maximum is
\code{Inf}, or every element is \code{-Inf}, every element of the result is
\code{NaN}, as base R's formula gives (\code{Inf - Inf}), although the
limit would be a one-hot or uniform vector. None of these warn. Base R's
formula gives \code{NA} or \code{NaN} depending on the order of the
missing values; these functions always prefer \code{NA}.

\strong{Implementations.} The SIMD implementations use the 1-ULP
exponential of the bundled SLEEF library, and \code{"none"} the C math
library's \code{exp()}, as base R does; it is the reference the others are
tested against. An implementation built without SLEEF (see
\code{\link{simd_available}}) uses the next lower one that has it,
ultimately \code{"none"}. SLEEF has no faster, less accurate exponential,
so option \code{rsimd.math_accuracy} (see
\code{\link{simd_math_accuracy}()}) makes no difference to these
functions. Integer and logical input is converted to
double, an integer \code{NA} becoming \code{NA_real_}; complex and
\code{integer64} input is not supported.

\rsimdPrecision

\rsimdAttributes
}
\value{
A double vector of the length of \code{x}.
}
\seealso{
\code{\link{plogis}}, \code{\link{simd_exp}}, \code{\link{simd_tanh}},
\code{\link{simd_precision}}, \link{rsimd_options}.
}
\examples{
simd_sigmoid(c(-800, -1, 0, 1, 800))
plogis(-720)
simd_sigmoid(-720)

x <- c(1, 2, 3)
simd_softmax(x)
exp(x - max(x)) / sum(exp(x - max(x)))
simd_log_softmax(c(1000, 1001))
simd_softmax(c(-Inf, 1, 2))
simd_softmax(c(1, Inf))

old <- simd_precision("compensated")
sum(simd_softmax(seq(-10, 10, length.out = 1e5)))
simd_precision(old)
}
