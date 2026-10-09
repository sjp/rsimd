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
\section{Main functions}{
\describe{
  \item{Reductions}{\code{\link{simd_sum}}, \code{\link{simd_min}},
    \code{\link{simd_which_min}}.}
  \item{Linear algebra and statistics}{\code{\link{simd_dot}},
    \code{\link{simd_var}}.}
  \item{Scans}{\code{\link{simd_cumsum}}.}
  \item{Elementwise arithmetic}{\code{\link{simd_add}},
    \code{\link{simd_fma}}, \code{\link{simd_round}}.}
  \item{Predicates and comparisons}{\code{\link{simd_is_na}},
    \code{\link{simd_eq}}, \code{\link{simd_hamming}}.}
  \item{Logic and bits}{\code{\link{simd_and}}, \code{\link{simd_bit_and}},
    \code{\link{simd_popcount}}.}
  \item{Conversions}{\code{\link{simd_as_integer}}.}
  \item{Elementary functions}{\code{\link{simd_exp}},
    \code{\link{simd_sin}}, \code{\link{simd_math_accuracy}}.}
  \item{Floating-point extras}{\code{\link{simd_ilogb}},
    \code{\link{simd_recip_approx}}.}
  \item{Machine learning}{\code{\link{simd_softmax}}.}
  \item{Complex vectors}{\code{\link{simd_conj}}.}
  \item{64-bit integers}{\link{rsimd-integer64}.}
  \item{The \code{simd_vec} class}{\code{\link{simd_vec}}.}
  \item{Implementation selection}{\code{\link{simd_use}},
    \code{\link{simd_precision}}, \code{\link{simd_cpu_features}}.}
}
}
\section{Vignettes}{
\describe{
  \item{\code{vignette("rsimd", package = "rsimd")}}{Getting started with
    rsimd.}
  \item{\code{vignette("implementation-selection", package = "rsimd")}}{Choosing
    an implementation.}
  \item{\code{vignette("numerical-semantics", package = "rsimd")}}{Numerical
    semantics and precision.}
  \item{\code{vignette("benchmarks", package = "rsimd")}}{Benchmarks.}
}
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
\concept{SIMD}
