\name{simd_dot}
\alias{simd_dot}
\alias{simd_dist}
\alias{simd_cosine}
\title{Dot Products, Euclidean Distances and Cosine Similarity}
\usage{
simd_dot(x, y, na.rm = FALSE, na_check = NULL)

simd_dist(x, y, na.rm = FALSE, na_check = NULL)

simd_cosine(x, y, na_check = NULL)
}
\arguments{
\item{x, y}{Double, integer or logical vectors of the same length. The
types may differ.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove the pairs
\code{(x[i], y[i])} in which either element is missing (\code{NA} or
\code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} and \code{y} have no missing values, so they are not looked for; see
Details.}
}
\description{
\code{simd_dot(x, y)} is the dot product \code{sum(x * y)},
\code{simd_dist(x, y)} the Euclidean distance \code{sqrt(sum((x - y)^2))}
and \code{simd_cosine(x, y)} the cosine similarity
\code{sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))}, computed in one pass
with SIMD instructions.
}
\details{
\code{x} and \code{y} must have the same length: unlike the elementwise
functions, these never recycle a length-one argument. Integer and logical
elements are converted to double as they are read.

Missing values are handled per pair: \code{(x[i], y[i])} is missing when
either element is, and \code{na.rm = TRUE} drops the pair from both sums
(base R has no equivalent; \code{sum(x * y, na.rm = TRUE)} does the same for
\code{simd_dot()}). \code{simd_cosine()} has no \code{na.rm}.

The distance of two empty vectors is 0. The cosine similarity is \code{NaN} for empty input, when
either vector is all zero, and when any element is infinite. It is not
clamped to \eqn{[-1, 1]}, so rounding can take it just past either end.
Like \code{\link{simd_norm}()}, these do not guard against overflow of the
sums of squares.

\rsimdSemantics

\rsimdPrecision

All implementations have their own kernels for these functions, except that
on 32-bit ARM they use the scalar \code{"none"} code.
}
\value{
A double. \code{simd_dot()} and \code{simd_dist()} of empty vectors are 0.
}
\seealso{
\code{\link{simd_sum_sq}}, \code{\link{crossprod}},
\code{\link{simd_precision}}, \link{rsimd_options}.
}
\examples{
simd_dot(1:3, c(0.5, 1, 2))
simd_dist(c(0, 0), c(3, 4))
simd_cosine(c(1, 0), c(1, 1))
simd_dot(c(1, NA, 3), c(2, 5, 1), na.rm = TRUE)
}
