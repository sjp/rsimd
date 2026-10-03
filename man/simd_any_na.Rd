\name{simd_any_na}
\alias{simd_any_na}
\alias{simd_count_na}
\alias{simd_which_na}
\title{Find Missing Values}
\usage{
simd_any_na(x)

simd_count_na(x)

simd_which_na(x)
}
\arguments{
\item{x}{A double, integer, logical, complex or raw vector.}
}
\description{
Whether a vector has missing values, how many, and where, like
\code{\link{anyNA}(x)}, \code{sum(is.na(x))} and \code{which(is.na(x))},
computed with SIMD instructions and without creating the logical vector
\code{is.na(x)}.
}
\details{
As for \code{\link{is.na}()}, \code{NaN} counts as missing, and a complex
value is missing when either its real or its imaginary part is. Raw vectors
have no missing values. \code{simd_any_na()} stops at the first missing
value, and none of the functions reads a vector that R already knows to have
no missing values, such as \code{1:n}.

The results are exact and the same for every implementation. All
implementations have their own kernels, except for complex vectors, and on
32-bit ARM for double vectors, which use the scalar \code{"none"} code.
}
\value{
\code{simd_any_na()}: \code{TRUE} or \code{FALSE}.

\code{simd_count_na()}: the number of missing values, as a double (so that
it can exceed \code{.Machine$integer.max}).

\code{simd_which_na()}: the indices of the missing values, in increasing
order: an integer vector, or a double vector when an index exceeds
\code{.Machine$integer.max}.
}
\seealso{
\code{\link{anyNA}}, \code{\link{is.na}}.
}
\examples{
x <- c(1, NA, 3, NaN)
simd_any_na(x)
simd_count_na(x)
simd_which_na(x)
}
