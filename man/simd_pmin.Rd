\name{simd_pmin}
\alias{simd_pmin}
\alias{simd_pmax}
\alias{simd_pmin_num}
\alias{simd_pmax_num}
\alias{simd_clamp}
\title{Parallel Minimum and Maximum, Clamping}
\usage{
simd_pmin(x, y, na.rm = FALSE)

simd_pmax(x, y, na.rm = FALSE)

simd_pmin_num(x, y)

simd_pmax_num(x, y)

simd_clamp(x, lo, hi)
}
\arguments{
\item{x, y}{Double, integer or logical vectors.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: ignore a missing operand?}
\item{lo, hi}{Double, integer or logical vectors of lower and upper bounds.}
}
\description{
\code{simd_pmin(x, y)} and \code{simd_pmax(x, y)} are base R's
\code{pmin(x, y)} and \code{pmax(x, y)} for two vectors, computed with SIMD
instructions. \code{simd_pmin_num()} and \code{simd_pmax_num()} are the
same with \code{na.rm = TRUE}. \code{simd_clamp(x, lo, hi)} is
\code{pmin(pmax(x, lo), hi)}. They correspond to \code{Min}, \code{Max},
\code{MinNumber}, \code{MaxNumber} and \code{Clamp} of .NET's
\code{TensorPrimitives}.
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
missing. \code{TensorPrimitives} has no \code{NA}: there \code{Max}
propagates \code{NaN} and \code{MaxNumber} ignores it, which is what these
do with \code{NaN}, and \code{NA} is treated the same way.

\code{simd_clamp()} errors if \code{lo > hi} at any position (missing bounds
are not compared). A missing bound gives a missing result, following
\code{pmin(pmax(x, lo), hi)}.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double operations use the scalar \code{"none"} code.
}
\value{
A vector of the length of the longest operand: integer when every operand is
integer or logical, otherwise double.
}
\seealso{
\code{\link{pmin}}, \code{\link{simd_min}} (reductions),
\link{rsimd_options}.
}
\examples{
simd_pmax(c(1, 5, NA), 3)
simd_pmax(c(1, 5, NA), 3, na.rm = TRUE)
simd_pmin(1:5, 3L)
simd_clamp(c(-2, 0.5, 7), 0, 1)
}
