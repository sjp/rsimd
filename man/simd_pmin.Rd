\name{simd_pmin}
\alias{simd_pmin}
\alias{simd_pmax}
\alias{simd_pmin_num}
\alias{simd_pmax_num}
\alias{simd_clamp}
\title{Parallel Minimum and Maximum, Clamping}
\usage{
simd_pmin(x, y, ..., na.rm = FALSE)

simd_pmax(x, y, ..., na.rm = FALSE)

simd_pmin_num(x, y)

simd_pmax_num(x, y)

simd_clamp(x, lo, hi)
}
\arguments{
\item{x, y}{Double, integer, logical or \pkg{bit64} \code{integer64}
vectors.}
\item{\dots}{Must be empty. \code{simd_pmin()} and \code{simd_pmax()} take
two vectors, and \code{na.rm} must be named, so that an extra vector is not
read as \code{na.rm}. Nest calls for more: unlike \code{pmin(x, y, z)},
\code{simd_pmin(x, y, z)} is an error.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: ignore a missing operand?}
\item{lo, hi}{Double, integer, logical or integer64 vectors of lower and
upper bounds.}
}
\description{
\code{simd_pmin(x, y)} and \code{simd_pmax(x, y)} are base R's
\code{pmin(x, y)} and \code{pmax(x, y)} for two vectors, computed with SIMD
instructions. \code{simd_pmin_num()} and \code{simd_pmax_num()} are the
same with \code{na.rm = TRUE}. \code{simd_clamp(x, lo, hi)} is
\code{pmin(pmax(x, lo), hi)}.
}
\details{
\rsimdBroadcast
Unlike base R's \code{pmin()}, these take exactly two operands (three for
\code{simd_clamp()}).

The results are identical to base R's, bit for bit. Each element of the
result is one of the operands: a tie returns \code{x}, so
\code{simd_pmin(0, -0)} is \code{0} and \code{simd_pmin(-0, 0)} is
\code{-0}. With \code{na.rm = FALSE} a missing operand is returned as it is
(\code{NA} or \code{NaN}), and when both are missing the second one is:
\code{simd_pmax(NaN, NA)} is \code{NA} and \code{simd_pmax(NA, NaN)} is
\code{NaN}. With \code{na.rm = TRUE} (and in the \code{_num} functions) a
missing operand is ignored, and the second one is returned when both are
missing.

\code{simd_clamp()} errors if \code{lo > hi} at any position (missing bounds
are not compared). A missing bound gives a missing result, following
\code{pmin(pmax(x, lo), hi)}.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double operations use the scalar \code{"none"} code.
Integer64 operands (with integer64, integer or logical ones) give integer64
results; with a double operand they are converted to double, with a warning
(see \link{rsimd-integer64}). Complex operands are not taken.
}
\value{
A vector of the length of the longest operand: integer when every operand is
integer or logical, otherwise double.
}
\seealso{
\code{\link{pmin}}, \code{\link{simd_min}} (reductions),
\link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_pmax(c(1, 5, NA), 3)
simd_pmax(c(1, 5, NA), 3, na.rm = TRUE)
simd_pmin(1:5, 3L)
simd_clamp(c(-2, 0.5, 7), 0, 1)
}
