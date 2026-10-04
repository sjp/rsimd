\name{simd_as_integer}
\alias{simd_as_integer}
\alias{simd_as_double}
\alias{simd_as_logical}
\alias{simd_as_raw}
\alias{simd_as_integer64}
\title{Type Conversions}
\usage{
simd_as_integer(x, mode = c("checked", "saturating", "truncating"))

simd_as_double(x)

simd_as_logical(x)

simd_as_raw(x, mode = c("checked", "saturating", "truncating"))

simd_as_integer64(x, mode = c("checked", "saturating", "truncating"))
}
\arguments{
\item{x}{A double, integer, logical, raw or \pkg{bit64} \code{integer64}
vector.}
\item{mode}{How values that do not fit are converted: \code{"checked"}
(the default, as base R), \code{"saturating"} or \code{"truncating"}; see
Details. It matters only for doubles and integer64 converted to integer,
doubles converted to integer64 and conversions to raw.}
}
\description{
\code{as.integer(x)}, \code{as.double(x)}, \code{as.logical(x)},
\code{as.raw(x)} and \pkg{bit64}'s \code{as.integer64(x)} computed with SIMD
instructions, with the conversion modes
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

\strong{To integer64.} \code{simd_as_integer64()} gives an
\code{"integer64"} vector (see \link{rsimd-integer64}) without needing
\pkg{bit64}. Doubles are truncated toward zero, then
\describe{
\item{\code{"checked"}}{gives \code{NA} for \code{NaN}, \code{NA} and
values outside \eqn{(-2^{63}, 2^{63})}, with \pkg{bit64}'s warning
\code{"NAs produced by integer64 overflow"} once (but not for \code{NA} and
\code{NaN}), as \code{bit64::as.integer64()}; \eqn{-2^{63}} is the bit
pattern of \code{NA}, so it is out of range.}
\item{\code{"saturating"}}{clamps to \eqn{\pm(2^{63} - 1)}, without a
warning; \code{NaN} and \code{NA} give \code{NA}.}
\item{\code{"truncating"}}{wraps around modulo \eqn{2^{64}}, without a
warning; a result of \eqn{-2^{63}} reads as \code{NA}, and \code{NaN},
\code{NA} and infinities give \code{NA}.}
}
Integer, logical and raw input converts exactly.

From integer64, \code{simd_as_integer()} gives \code{NA} for values outside
\code{[-.Machine$integer.max, .Machine$integer.max]} in \code{"checked"}
mode, with \pkg{bit64}'s warning \code{"NAs produced by integer overflow"}
once; \code{"saturating"} clamps to that range and \code{"truncating"} keeps
the low 32 bits. \code{simd_as_raw()} treats integer64 as it treats
integers.

\strong{To double and logical.} \code{simd_as_double()} is exact for
integer, logical and raw input, with \code{NA} for \code{NA}; integer64 is
exact up to \eqn{2^{53}} in magnitude and rounded to nearest beyond, with
\pkg{bit64}'s warning \code{"integer precision lost while converting to
double"} once when any value has a magnitude of \eqn{2^{53}} or more.
\code{simd_as_logical()} gives \code{FALSE} for zero, \code{TRUE} for other
numbers and bytes, and \code{NA} for \code{NA} and \code{NaN}.

Input of the target type is returned as it is, without attributes.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM all conversions use the scalar \code{"none"} code (and see
\link{rsimd-integer64} for the integer64 conversions). Complex input is not
supported.
}
\value{
An integer, double, logical, raw or integer64 vector of the length of
\code{x}.
}
\seealso{
\code{\link{as.integer}}, \code{\link{as.raw}}, \code{\link{simd_round}},
\link{rsimd-integer64}.
}
\examples{
simd_as_integer(c(2.9, -2.9, 2^31, NaN))
simd_as_integer(c(2^31, -2^31, 1e10), mode = "saturating")
simd_as_integer(c(2^31, 2^32 + 5), mode = "truncating")
simd_as_double(c(1L, NA))
simd_as_logical(c(0, 2.5, NaN))
simd_as_raw(c(255, 256, -1), mode = "saturating")
simd_as_raw(c(255, 256, -1), mode = "truncating")
x <- simd_as_integer64(c(2^62, -7.9, NA))
simd_as_integer(x, mode = "saturating")
simd_as_integer64(c(Inf, -Inf), mode = "saturating")
}
