\name{simd_precision}
\alias{simd_precision}
\title{Accuracy of Floating-Point Sums}
\usage{
simd_precision(mode)
}
\arguments{
\item{mode}{One of \code{"fast"}, \code{"pairwise"} or \code{"compensated"}.}
}
\description{
Gets or sets how \pkg{rsimd} accumulates floating-point sums: as fast as
possible, or more accurately at some cost in speed.
}
\details{
\rsimdPrecision

The mode is stored in the option \code{rsimd.precision} (see
\link{rsimd_options}). Setting that option directly works too; an invalid
value is reported by the next function that uses it.
}
\value{
Without \code{mode}, the current mode. With \code{mode}, the previous mode,
invisibly; an invalid \code{mode} is an error and leaves the mode unchanged.
}
\seealso{
\link{rsimd_options} for all options.
}
\examples{
simd_precision()
old <- simd_precision("compensated")
simd_precision()
simd_precision(old)
}
