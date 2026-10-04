\name{simd_as_integer}
\alias{simd_as_integer}
\alias{simd_as_double}
\alias{simd_as_logical}
\alias{simd_as_raw}
\title{Type Conversions}
\usage{
simd_as_integer(x, mode = c("checked", "saturating", "truncating"))

simd_as_double(x)

simd_as_logical(x)

simd_as_raw(x, mode = c("checked", "saturating", "truncating"))
}
\arguments{
\item{x}{A double, integer, logical or raw vector.}
\item{mode}{How values that do not fit are converted: \code{"checked"}
(the default, as base R), \code{"saturating"} or \code{"truncating"}; see
Details. It matters only for doubles converted to integer and for
conversions to raw.}
}
\description{
\code{as.integer(x)}, \code{as.double(x)}, \code{as.logical(x)} and
\code{as.raw(x)} computed with SIMD instructions, with the conversion modes
of \code{ConvertChecked}, \code{ConvertSaturating} and
\code{ConvertTruncating} of .NET's \code{TensorPrimitives}.
}
\details{
\strong{To integer.} Doubles are truncated toward zero, then
\describe{
\item{\code{"checked"}}{gives \code{NA} for \code{NaN}, \code{NA} and
values outside the integer range, with base R's warning \code{"NAs
introduced by coercion to integer range"} once, exactly like
\code{as.integer()}. As \code{-2147483648} is \code{NA_integer_}, doubles
up to \code{-2147483648} are out of range, while \code{-2147483647.5}
becomes \code{-2147483647}.}
\item{\code{"saturating"}}{clamps to
\code{[-.Machine$integer.max, .Machine$integer.max]}, without a warning;
\code{NaN} and \code{NA} give \code{NA}.}
\item{\code{"truncating"}}{wraps around modulo \eqn{2^{32}} (two's
complement), without a warning; a result of \code{-2147483648} reads as
\code{NA}, and \code{NaN}, \code{NA} and infinities give \code{NA}.}
}
Logical and raw input converts exactly.

\strong{To raw.} Values are truncated toward zero, then \code{"checked"}
gives \code{00} for missing values and values outside 0 to 255, with base
R's warning \code{"out-of-range values treated as 0 in coercion to raw"}
once, exactly like \code{as.raw()} (which also warns \code{"NAs introduced
by coercion to integer range"} for doubles outside the integer range);
\code{"saturating"} clamps to 0 to 255; \code{"truncating"} keeps the low 8
bits (modulo 256). The last two never warn and give \code{00} for missing
values.

\strong{To double and logical.} \code{simd_as_double()} is exact for
integer, logical and raw input, with \code{NA} for \code{NA}.
\code{simd_as_logical()} gives \code{FALSE} for zero, \code{TRUE} for other
numbers and bytes, and \code{NA} for \code{NA} and \code{NaN}.

Input of the target type is returned as it is, without attributes.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM all conversions use the scalar \code{"none"} code. Complex input
is not supported; \pkg{bit64} \code{integer64} input is not supported yet,
except by \code{simd_as_raw}, which will not take it.
}
\value{
An integer, double, logical or raw vector of the length of \code{x}.
}
\seealso{
\code{\link{as.integer}}, \code{\link{as.raw}}, \code{\link{simd_round}}.
}
\examples{
simd_as_integer(c(2.9, -2.9, 2^31, NaN))
simd_as_integer(c(2^31, -2^31, 1e10), mode = "saturating")
simd_as_integer(c(2^31, 2^32 + 5), mode = "truncating")
simd_as_double(c(1L, NA))
simd_as_logical(c(0, 2.5, NaN))
simd_as_raw(c(255, 256, -1), mode = "saturating")
simd_as_raw(c(255, 256, -1), mode = "truncating")
}
