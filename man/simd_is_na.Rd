\name{simd_is_na}
\alias{simd_is_na}
\alias{simd_is_na_any}
\alias{simd_is_na_all}
\alias{simd_is_nan}
\alias{simd_is_nan_any}
\alias{simd_is_nan_all}
\alias{simd_is_finite}
\alias{simd_is_finite_any}
\alias{simd_is_finite_all}
\alias{simd_is_infinite}
\alias{simd_is_infinite_any}
\alias{simd_is_infinite_all}
\alias{simd_is_negative}
\alias{simd_is_negative_any}
\alias{simd_is_negative_all}
\alias{simd_is_zero}
\alias{simd_is_zero_any}
\alias{simd_is_zero_all}
\title{Missing, Finite, Negative and Zero Values}
\usage{
simd_is_na(x)
simd_is_na_any(x)
simd_is_na_all(x)

simd_is_nan(x)
simd_is_nan_any(x)
simd_is_nan_all(x)

simd_is_finite(x)
simd_is_finite_any(x)
simd_is_finite_all(x)

simd_is_infinite(x)
simd_is_infinite_any(x)
simd_is_infinite_all(x)

simd_is_negative(x)
simd_is_negative_any(x)
simd_is_negative_all(x)

simd_is_zero(x)
simd_is_zero_any(x)
simd_is_zero_all(x)
}
\arguments{
\item{x}{A double, integer, logical, raw or complex vector (neither raw nor
complex for \code{simd_is_negative} and \code{simd_is_zero}).}
}
\description{
\code{is.na(x)}, \code{is.nan(x)}, \code{is.finite(x)} and
\code{is.infinite(x)}, and tests for negative and zero values, computed with
SIMD instructions. Each test comes in three forms: elementwise, \code{_any}
(is it true for any element?) and \code{_all} (is it true for every
element?). They correspond to \code{IsNaN}, \code{IsFinite},
\code{IsInfinity}, \code{IsNegative} and \code{IsZero} and their
\code{All}/\code{Any} forms in .NET's \code{TensorPrimitives}.
}
\details{
The tests never return \code{NA}. For doubles they follow base R:
\code{simd_is_na()} is \code{TRUE} for \code{NA} and \code{NaN},
\code{simd_is_nan()} for \code{NaN} but not \code{NA} (whatever the
payload or sign of the \code{NaN}), \code{simd_is_finite()} for numbers
that are neither missing nor infinite. \code{simd_is_zero()} is
\code{x == 0}, so it is \code{TRUE} for \code{-0} too.
\code{simd_is_negative()} tests the sign bit of numbers, as
\code{TensorPrimitives.IsNegative} and IEEE 754's \code{isSignMinus} do: it
is \code{TRUE} for \code{-0} and \code{-Inf}, and \code{FALSE} for
\code{NA} and \code{NaN}, so it differs from \code{x < 0} for \code{-0}.

A complex number is \code{NA} (\code{simd_is_na()}) when either part is
\code{NA} or \code{NaN}, \code{NaN} (\code{simd_is_nan()}) when either
part is \code{NaN} but not \code{NA}, finite when both parts are, and
infinite when either part is, as in base R: \code{complex(real = NaN,
imaginary = NA)} is both \code{NA} and \code{NaN}.

Integers and logicals are \code{NA} where they are \code{NA} and finite
elsewhere; they are never \code{NaN} or infinite. Raw bytes are never
missing, and, as in base R's \code{is.finite()}, are not finite either.

The \code{_any} forms are \code{FALSE} for empty input and the \code{_all}
forms \code{TRUE}, as \code{any()} and \code{all()} are in R
(\code{TensorPrimitives} returns \code{false} for both). They may stop
reading at the first element that decides the answer.
\code{simd_is_na_any(x)} is the same as \code{\link{simd_any_na}(x)}.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double and complex tests use the scalar \code{"none"}
code. \pkg{bit64} \code{integer64} input is not supported yet for
\code{simd_is_na}; the other tests do not take it.
}
\value{
A logical vector of the length of \code{x}, without \code{NA}, or a single
\code{TRUE} or \code{FALSE} for the \code{_any} and \code{_all} forms.
}
\seealso{
\code{\link{is.na}}, \code{\link{is.finite}}, \code{\link{simd_any_na}},
\code{\link{simd_eq}}.
}
\examples{
x <- c(1, NA, NaN, -0, Inf, -2)
simd_is_na(x)
simd_is_nan(x)
simd_is_finite(x)
simd_is_negative(x)
simd_is_zero_any(x)
simd_is_finite_all(numeric(0))
}
