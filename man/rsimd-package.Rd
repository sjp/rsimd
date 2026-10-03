\docType{package}
\name{rsimd-package}
\alias{rsimd}
\alias{rsimd-package}
\title{rsimd: SIMD-Accelerated Vector Primitives}
\description{
Fast vectorised primitives for R atomic vectors, implemented with single
instruction, multiple data (SIMD) instructions on x86-64 and arm64. Kernels
are written with the portable \sQuote{SIMDe} library and use \sQuote{SLEEF}
for vectorised elementary functions. CPU features are detected at run time
and the best supported implementation is chosen automatically, while a
specific implementation, or a plain scalar one, can be selected by name.
}
\section{Options}{
\describe{
  \item{\code{rsimd.impl}}{The implementation to use. Defaults to
    \code{"auto"}, or to the value of the \env{RSIMD_IMPL} environment
    variable when it is set before the package is loaded.}
  \item{\code{rsimd.precision}}{Accumulation mode for floating-point
    reductions. Defaults to \code{"fast"}.}
  \item{\code{rsimd.na_check}}{Whether integer and logical kernels check for
    \code{NA}. Defaults to \code{TRUE}.}
}
}
\seealso{
\code{\link{simd_version}}
}
\author{
Simon Potter \email{simon@sjp.co.nz}
}
\keyword{internal}
