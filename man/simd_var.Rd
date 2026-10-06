\name{simd_var}
\alias{simd_var}
\alias{simd_sd}
\title{Variance and Standard Deviation}
\usage{
simd_var(x, na.rm = FALSE, na_check = NULL)

simd_sd(x, na.rm = FALSE, na_check = NULL)
}
\arguments{
\item{x}{A double, integer or logical vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The sample variance and standard deviation of the elements of a vector, as
\code{\link[stats]{var}()} and \code{\link[stats]{sd}()} compute them for a
single vector (denominator \eqn{n - 1}), using SIMD instructions. They
correspond to \code{StdDev} of .NET's \code{TensorPrimitives}, which leaves
the denominator unspecified.
}
\details{
Like base R, these take two passes over the data: the mean (computed as
\code{\link{simd_mean}()} does, in the same precision mode), then the sum of
squared deviations from it. They do not use the one-pass formula
\eqn{E[x^2] - E[x]^2}, which loses all precision when the mean is large
compared with the spread. Integer and logical elements are converted to
double as they are read.

With fewer than two elements (after \code{na.rm}) the result is \code{NA}.
As in base R, a missing value gives \code{NA} without \code{na.rm}, whether
it is \code{NA} or \code{NaN}.

When finite values are so large that their mean overflows to \code{Inf} or
\code{-Inf}, the result is \code{Inf} on every implementation, and an
infinite element gives \code{NaN}. Base R computes the mean in long
double, which has a wider range than double on most platforms, so it can
give a finite result there (\code{var(rep(1e308, 3))} is 0 on x86-64).

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
