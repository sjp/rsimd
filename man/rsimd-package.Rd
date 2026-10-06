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
The options that control the package (\code{rsimd.impl},
\code{rsimd.precision}, \code{rsimd.na_check} and
\code{rsimd.math_accuracy}) are described in \link{rsimd_options}.
}
\seealso{
\code{\link{simd_version}}
}
\author{
Simon Potter \email{simon@sjp.co.nz}
}
\keyword{internal}
