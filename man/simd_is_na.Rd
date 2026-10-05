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
\alias{simd_is_normal}
\alias{simd_is_normal_any}
\alias{simd_is_normal_all}
\alias{simd_is_subnormal}
\alias{simd_is_subnormal_any}
\alias{simd_is_subnormal_all}
\alias{simd_is_whole}
\alias{simd_is_whole_any}
\alias{simd_is_whole_all}
\alias{simd_is_even}
\alias{simd_is_even_any}
\alias{simd_is_even_all}
\alias{simd_is_odd}
\alias{simd_is_odd_any}
\alias{simd_is_odd_all}
\alias{simd_is_pow2}
\alias{simd_is_pow2_any}
\alias{simd_is_pow2_all}
\title{Missing, Finite, Negative and Zero Values, and Number Classes}
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

simd_is_normal(x)
simd_is_normal_any(x)
simd_is_normal_all(x)

simd_is_subnormal(x)
simd_is_subnormal_any(x)
simd_is_subnormal_all(x)

simd_is_whole(x)
simd_is_whole_any(x)
simd_is_whole_all(x)

simd_is_even(x)
simd_is_even_any(x)
simd_is_even_all(x)

simd_is_odd(x)
simd_is_odd_any(x)
simd_is_odd_all(x)

simd_is_pow2(x)
simd_is_pow2_any(x)
simd_is_pow2_all(x)
}
\arguments{
\item{x}{A double, integer, logical, raw, complex or \pkg{bit64}
\code{integer64} vector. Only \code{simd_is_na}, \code{simd_is_nan},
\code{simd_is_finite} and \code{simd_is_infinite} take raw and complex
vectors.}
}
\description{
\code{is.na(x)}, \code{is.nan(x)}, \code{is.finite(x)} and
\code{is.infinite(x)}, tests for negative and zero values, and tests of
number classes (normal, subnormal, whole, even, odd, power of two),
computed with SIMD instructions. Each test comes in three forms:
elementwise, \code{_any} (is it true for any element?) and \code{_all}
(is it true for every element?). They correspond to \code{IsNaN},
\code{IsFinite}, \code{IsInfinity}, \code{IsNegative}, \code{IsZero},
\code{IsNormal}, \code{IsSubnormal}, \code{IsInteger},
\code{IsEvenInteger}, \code{IsOddInteger} and \code{IsPow2} and their
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

The number classes are never \code{TRUE} for \code{NA}, \code{NaN} or
an infinity:
\itemize{
\item \code{simd_is_normal()}: a finite, non-zero double whose magnitude
is at least \code{.Machine$double.xmin}; an integer that is not 0;
\item \code{simd_is_subnormal()}: a non-zero double smaller in magnitude
than \code{.Machine$double.xmin}, such as \code{5e-324}; never an
integer;
\item \code{simd_is_whole()}: \code{x == trunc(x)}, without tolerance;
every integer;
\item \code{simd_is_even()} and \code{simd_is_odd()}: whole numbers
with \code{x \%\% 2} equal to 0 or 1 (so negative numbers too, and
\code{-0} is even). Every double of magnitude \eqn{2^{53}} or more is
even;
\item \code{simd_is_pow2()}: \eqn{2^k} for an integer \eqn{k}, which for
doubles includes 0.5 and the subnormal \eqn{2^{-1074}}; never zero or
negative.
}
\code{simd_is_whole()} is named so because base R's \code{is.integer()}
tests the storage type. Logical values are tested as the integers 0 and 1,
so \code{TRUE} is normal, odd and a power of two.

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
code. \pkg{bit64} \code{integer64} input is tested like integers: it is
\code{NA} when it is \eqn{-2^{63}}, never \code{NaN} or infinite (see
\link{rsimd-integer64}); it is never subnormal.
}
\value{
A logical vector of the length of \code{x}, without \code{NA}, or a single
\code{TRUE} or \code{FALSE} for the \code{_any} and \code{_all} forms.
}
\seealso{
\code{\link{is.na}}, \code{\link{is.finite}}, \code{\link{simd_any_na}},
\code{\link{simd_eq}}, \link{rsimd-integer64}.
}
\examples{
x <- c(1, NA, NaN, -0, Inf, -2)
simd_is_na(x)
simd_is_nan(x)
simd_is_finite(x)
simd_is_negative(x)
simd_is_zero_any(x)
simd_is_finite_all(numeric(0))
simd_is_subnormal(c(5e-324, 2.2e-308, 1))
simd_is_whole(c(2, 2.5, -0, Inf))
simd_is_even(c(-4L, 3L, NA))
simd_is_pow2(c(0.25, 1, 3, 1024))
}
