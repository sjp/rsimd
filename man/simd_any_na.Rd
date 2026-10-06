\name{simd_any_na}
\alias{simd_any_na}
\alias{simd_count_na}
\alias{simd_which_na}
\alias{simd_which}
\alias{simd_count}
\title{Find Missing Values and TRUE Elements}
\usage{
simd_any_na(x)

simd_count_na(x)

simd_which_na(x)

simd_which(x)

simd_count(x, na.rm = FALSE)
}
\arguments{
\item{x}{A double, integer, logical, complex, raw or \pkg{bit64}
\code{integer64} vector; a logical vector for \code{simd_which()} and
\code{simd_count()}.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: ignore \code{NA} elements?}
}
\description{
Whether a vector has missing values, how many, and where, like
\code{\link{anyNA}(x)}, \code{sum(is.na(x))} and \code{which(is.na(x))},
computed with SIMD instructions and without creating the logical vector
\code{is.na(x)}; and how many elements of a logical vector are
\code{TRUE}, and where, like \code{\link{sum}(x)} and
\code{\link{which}(x)}.
}
\details{
As for \code{\link{is.na}()}, \code{NaN} counts as missing, and a complex
value is missing when either its real or its imaginary part is. Raw vectors
have no missing values, and an integer64 value is missing when it is
\eqn{-2^{63}} (see \link{rsimd-integer64}). \code{simd_any_na()} stops at the first missing
value, and none of the functions reads a vector that R already knows to have
no missing values, such as \code{1:n}.

\code{simd_which()} and \code{simd_count()} take only logical vectors; any
other type is an error, where base R's \code{sum()} would add numbers and
\code{which()} would also error. As in base R, \code{simd_which()} ignores
\code{NA} elements, and \code{simd_count()} returns \code{NA} if there is
one, unless \code{na.rm = TRUE}, which counts the \code{TRUE} elements
among the others. \code{simd_count()} does not look for \code{NA} in a
vector known to have none, such as a \code{\link{simd_vec}} marked
NA-free. Names are dropped:
\code{simd_which(c(a = TRUE))} is \code{1L}, not \code{c(a = 1L)}.

The results are exact and the same for every implementation. All
implementations have their own kernels, except that on 32-bit ARM double
and complex vectors use the scalar \code{"none"} code.
}
\value{
\code{simd_any_na()}: \code{TRUE} or \code{FALSE}.

\code{simd_count_na()}: the number of missing values, as a double (so that
it can exceed \code{.Machine$integer.max}).

\code{simd_which_na()} and \code{simd_which()}: the indices of the missing
values or of the \code{TRUE} elements, in increasing order: an integer
vector, or a double vector when an index exceeds
\code{.Machine$integer.max}.

\code{simd_count()}: the number of \code{TRUE} elements, an integer, or a
double when it exceeds \code{.Machine$integer.max}; \code{NA_integer_}
when \code{x} has an \code{NA} and \code{na.rm} is \code{FALSE}.
}
\seealso{
\code{\link{anyNA}}, \code{\link{is.na}}, \code{\link{which}},
\code{\link{simd_any}}, \link{rsimd-integer64}.
}
\examples{
x <- c(1, NA, 3, NaN)
simd_any_na(x)
simd_count_na(x)
simd_which_na(x)

b <- c(TRUE, FALSE, NA, TRUE)
simd_which(b)
simd_count(b)
simd_count(b, na.rm = TRUE)
}
