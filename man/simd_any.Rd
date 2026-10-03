\name{simd_any}
\alias{simd_any}
\alias{simd_all}
\title{Are Some or All Values True?}
\usage{
simd_any(x, na.rm = FALSE)

simd_all(x, na.rm = FALSE)
}
\arguments{
\item{x}{A logical, integer, double or raw vector.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: ignore missing values?}
}
\description{
Whether at least one, or every, element of a vector is \code{TRUE}, like
\code{\link{any}()} and \code{\link{all}()} for a single vector, computed
with SIMD instructions.
}
\details{
The result uses R's three-valued logic: \code{simd_any()} is \code{TRUE} if
an element is \code{TRUE}, otherwise \code{NA} if an element is \code{NA}
(and \code{na.rm} is \code{FALSE}), otherwise \code{FALSE};
\code{simd_all()} is \code{FALSE} if an element is \code{FALSE}, otherwise
\code{NA} if an element is \code{NA}, otherwise \code{TRUE}. Reading stops
once the answer is known. The empty vector gives \code{FALSE} and
\code{TRUE}.

Integer, double and raw values are read as logical values without being
converted: non-zero is \code{TRUE} and zero \code{FALSE}; an integer
\code{NA} and a double \code{NA} or \code{NaN} are \code{NA}. As in base R,
a non-empty double or raw \code{x} gives the warning \code{"coercing
argument of type 'double' to logical"} (or \code{'raw'}).

The results are the same for every implementation; all have their own
kernels, except that on 32-bit ARM the double version uses the scalar
\code{"none"} code. The precision mode does not apply.
}
\value{
\code{TRUE}, \code{FALSE} or \code{NA}.
}
\seealso{
\code{\link{any}}, \code{\link{all}}, \code{\link{simd_any_na}}.
}
\examples{
simd_any(c(FALSE, NA, TRUE))
simd_all(c(TRUE, NA))
simd_all(c(TRUE, NA), na.rm = TRUE)
simd_any(1:3 > 2)
}
