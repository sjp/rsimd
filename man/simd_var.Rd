\name{simd_var}
\alias{simd_var}
\alias{simd_sd}
\title{Variance and Standard Deviation}
\usage{
simd_var(x, ..., na.rm = FALSE, na_check = NULL)

simd_sd(x, ..., na.rm = FALSE, na_check = NULL)
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{\dots}{Must be empty. Each function takes one vector, and \code{na.rm}
and \code{na_check} must be named, so that an extra vector is not read as
\code{na.rm}. Unlike \code{sd(x, TRUE)}, \code{simd_sd(x, TRUE)} is an
error.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The sample variance and standard deviation of the elements of a vector, as
\code{\link[stats]{var}()} and \code{\link[stats]{sd}()} compute them for a
single vector (denominator \eqn{n - 1}), using SIMD instructions.
}
\details{
Like base R, these compute the mean first and then the sum of squared
deviations from it. For doubles the mean is the sum (accumulated in the
precision mode) divided by the count, refined in every mode by adding the
mean of the deviations from it, as base R and the \code{"pairwise"} and
\code{"compensated"} modes of \code{\link{simd_mean}()} do; so there are
three passes over the data, and a small spread is measured accurately even
next to a large mean. They do not use the one-pass formula
\eqn{E[x^2] - E[x]^2}, which loses all precision when the mean is large
compared with the spread. Integer and logical elements are converted to
double as they are read.

With fewer than two elements (after \code{na.rm}) the result is \code{NA}.
As in base R, a missing value gives \code{NA} without \code{na.rm}, whether
it is \code{NA} or \code{NaN}.

When finite values are so large that their sum overflows, the mean is
recomputed as in \code{\link{simd_mean}()}, so \code{simd_var(rep(1e308, 3))}
is 0 as in base R. A variance beyond the double range (say of
\code{c(1.7e308, -1.7e308)}) is \code{Inf}, and an infinite element
gives \code{NaN}.

\rsimdSemantics

\rsimdPrecision

All implementations have their own kernels for these functions, except that
on 32-bit ARM they use the scalar \code{"none"} code.

\code{stats::var()} and \code{stats::sd()} are not generic, so these are
separate functions rather than methods.
}
\value{
A double.
}
\seealso{
\code{\link[stats]{var}}, \code{\link[stats]{sd}}, \code{\link{simd_mean}},
\code{\link{simd_precision}}, \link{rsimd_options}.
}
\examples{
simd_var(1:10)
simd_sd(c(1, 2, 3, 4))
simd_var(c(1, NA, 3), na.rm = TRUE)
simd_var(1) # NA, as var(1)
}
