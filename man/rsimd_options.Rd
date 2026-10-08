\name{rsimd_options}
\alias{rsimd_options}
\title{Options and Environment Variables of rsimd}
\description{
The options and environment variables that control \pkg{rsimd}.
}
\section{Options}{
\describe{
\item{\code{rsimd.impl}}{The implementation in use: \code{"auto"} (the
default) or a tier id such as \code{"avx2"} or \code{"none"}. Set it with
\code{\link{simd_use}()}, which checks the value; setting the option
directly takes effect on the next call of an \pkg{rsimd} function.}
\item{\code{rsimd.precision}}{How floating-point sums are accumulated:
\code{"fast"} (the default), \code{"pairwise"} or \code{"compensated"}. See
\code{\link{simd_precision}()}.}
\item{\code{rsimd.math_accuracy}}{Which variants of the SLEEF library the
SIMD implementations of the elementary functions use: \code{"accurate"}
(the default) or \code{"fast"}, which is faster and less accurate in the
last bits. See \code{\link{simd_math_accuracy}()}.}
\item{\code{rsimd.na_check}}{\code{TRUE} (the default) or \code{FALSE}. Used
when the \code{na_check} argument is \code{NULL} (its default):
\code{FALSE} promises that inputs contain no missing values, so functions
skip looking for them. See \sQuote{Missing values and integer overflow}
below.}
}
The options are set to their defaults when the package is loaded, unless
they are already set. An option already set to a value that is not valid
gives a warning naming it, and is set to its default. A value that is not
valid set after loading makes the next call of a function that reads the
option an error naming it and how to reset it: with \code{\link{simd_use}()},
\code{\link{simd_precision}()} or \code{\link{simd_math_accuracy}()}, which
replace an invalid value, or for \code{rsimd.na_check} with
\code{options(rsimd.na_check = TRUE)}.
}
\section{Environment variables}{
\describe{
\item{\env{RSIMD_IMPL}}{Read when the package is loaded, if the option
\code{rsimd.impl} is not set: the initial implementation. An unknown or
unavailable value gives a warning and \code{"auto"} is used.}
\item{\env{RSIMD_CPU_FEATURES_MASK}}{For testing: a comma-separated list of
CPU feature names (see \code{\link{simd_cpu_features}()}) to treat as
absent, read when the package is loaded.}
\item{\env{RSIMD_DISABLE_TIERS}}{Read when the package is installed: a
comma-separated list of tiers not to compile.}
\item{\env{RSIMD_DEBUG_STRIDE}}{For testing: the number of elements
processed between checks for user interrupts (rounded up to a multiple of
128), read when the package is loaded. Long inputs are processed in pieces
of this size, so a small value exercises that splitting. The default is
\eqn{2^{20}}{2^20}.}
}
}
\section{Missing values and integer overflow}{
\rsimdSemantics
}
\seealso{
\code{\link{simd_use}()}, \code{\link{simd_precision}()},
\code{\link{simd_math_accuracy}()},
\code{\link{simd_cpu_features}()}.
}
