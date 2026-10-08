\name{simd_min}
\alias{simd_min}
\alias{simd_max}
\alias{simd_range}
\alias{simd_max_abs}
\alias{simd_min_abs}
\title{Minimum, Maximum and Range of Vector Elements and Magnitudes}
\usage{
simd_min(x, ..., na.rm = FALSE, na_check = NULL)

simd_max(x, ..., na.rm = FALSE, na_check = NULL)

simd_range(x, ..., na.rm = FALSE, na_check = NULL)

simd_max_abs(x, ..., na.rm = FALSE, na_check = NULL)

simd_min_abs(x, ..., na.rm = FALSE, na_check = NULL)
}
\arguments{
\item{x}{A double, integer, logical or \pkg{bit64} \code{integer64} vector;
also complex for \code{simd_max_abs()} and \code{simd_min_abs()}.}
\item{\dots}{Must be empty. Each function takes one vector, and \code{na.rm}
and \code{na_check} must be named, so that an extra vector is not read as
\code{na.rm}. Combine several vectors with \code{c()}: unlike
\code{max(x, y)}, \code{simd_max(x, y)} is an error.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} has no missing values, so they are not looked for; see Details.}
}
\description{
The smallest or largest element of a vector, or both, like
\code{\link{min}()}, \code{\link{max}()} and \code{\link{range}()} for a
single vector, and the largest or smallest magnitude, like
\code{max(abs(x))} and \code{min(abs(x))} without the vector
\code{abs(x)}, computed with SIMD instructions.
}
\details{
The results are exact and the same for every implementation, including the
sign of a zero result: like base R, \code{simd_min(c(0, -0))} is \code{0}
and \code{simd_min(c(-0, 0))} is \code{-0} (the first zero wins).

Without \code{na.rm}, a missing value makes the result \code{NA} (if there
is an \code{NA}) or \code{NaN}, wherever it is; \code{simd_range()} then
returns it twice. When nothing is left after removing missing values (or
\code{x} is empty), \code{simd_min()} returns \code{Inf} and
\code{simd_max()} \code{-Inf}, as doubles, with base R's warning;
\code{simd_range()} returns both, with both warnings.

\rsimdSemantics

The precision mode does not apply. All implementations have their own
kernels for these functions, except that on 32-bit ARM the double versions
use the scalar \code{"none"} code. Raw and complex vectors are rejected, as
by base R. Integer64 vectors give integer64 results; for no values these
are the largest and smallest integer64 values, with \pkg{bit64}'s warnings
(see \link{rsimd-integer64}).

\strong{Magnitudes.} \code{simd_max_abs(x)} is \code{max(abs(x))}, the
magnitude, which is never negative (a zero one is \code{0}, never
\code{-0}); the element itself is \code{x[\link{simd_which_max_abs}(x)]}.
\code{simd_min_abs(x)} is \code{min(abs(x))}. They follow the rules above:
\code{NA} or \code{NaN} without \code{na.rm}, \code{-Inf} and \code{Inf}
with base R's warnings for no values, integer results for integer and
logical \code{x} (\code{abs(NA_integer_)} stays \code{NA}; every other
integer has a magnitude that fits), and integer64 results for integer64
\code{x}. For complex \code{x} the magnitudes are the moduli that
\code{\link{simd_abs}()} computes (a double result), so
\code{simd_max_abs(z)} is \code{simd_max(simd_abs(z))}. These can differ
from \code{\link{Mod}(z)} in the last bit, and in fast mode by more (option
\code{rsimd.math_accuracy}, see \code{\link{simd_math_accuracy}}). Raw
vectors are rejected.
}
\value{
\code{simd_min()}, \code{simd_max()}, \code{simd_max_abs()} and
\code{simd_min_abs()}: a double for double (and, for the last two,
complex) \code{x}, an integer for integer and logical \code{x} (but see
above for empty input), an integer64 for integer64 \code{x}.
\code{simd_range()}: a vector of length 2 of the same type,
\code{c(simd_min(x), simd_max(x))}.
}
\seealso{
\code{\link{min}}, \code{\link{range}}, \code{\link{simd_which_min}},
\code{\link{simd_which_max_abs}}, \link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_min(c(3, 1, 2))
simd_max(c(3L, NA, 2L), na.rm = TRUE)
simd_range(c(TRUE, FALSE))
simd_min(c(1, NaN, NA)) # NA: NA wins over NaN

x <- c(-7, 2, NA, 5)
simd_max_abs(x, na.rm = TRUE)
x[simd_which_max_abs(x)]
simd_min_abs(complex(real = 3, imaginary = c(4, 0)))
}
