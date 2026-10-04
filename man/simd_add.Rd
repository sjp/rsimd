\name{simd_add}
\alias{simd_add}
\alias{simd_sub}
\alias{simd_mul}
\alias{simd_div}
\alias{simd_idiv}
\alias{simd_mod}
\alias{simd_add_wrap}
\alias{simd_sub_wrap}
\alias{simd_mul_wrap}
\title{Elementwise Arithmetic}
\usage{
simd_add(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_sub(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_mul(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_div(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_idiv(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_mod(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_add_wrap(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_sub_wrap(x, y, na_check = getOption("rsimd.na_check", TRUE))

simd_mul_wrap(x, y, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x, y}{Double, integer, logical or \pkg{bit64} \code{integer64}
vectors; not double for the \code{_wrap} functions. \code{simd_add()} and
\code{simd_sub()} also take complex vectors.}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} and \code{y} have no missing values, so they are not looked for; see
Details.}
}
\description{
\code{x + y}, \code{x - y}, \code{x * y}, \code{x / y}, \code{x \%/\% y} and
\code{x \%\% y} computed with SIMD instructions, and integer addition,
subtraction and multiplication that wrap around on overflow. They
correspond to \code{Add}, \code{Subtract}, \code{Multiply}, \code{Divide}
and \code{Remainder} of .NET's \code{TensorPrimitives} (whose
\code{Remainder} truncates; \code{simd_mod()} floors, as R's \code{\%\%}).
}
\details{
\rsimdBroadcast

\code{simd_add()} and \code{simd_sub()} take complex operands: if either
operand is complex, the other (double, integer or logical) is converted
with \code{as.complex()}, so an integer or logical \code{NA} becomes
\code{NA} with a zero imaginary part. The real and imaginary parts are then
added or subtracted separately, each following the rules for doubles
below, so a missing part makes only that part of the result missing.
Complex operands do not combine with \code{integer64} ones.

\pkg{bit64} \code{integer64} operands (with integer64, integer or logical
operands) give integer64 results computed exactly in 64 bits, with
\pkg{bit64}'s overflow warning; \code{simd_div()} converts them to double,
and a double operand makes them convert to double with a warning. See
\link{rsimd-integer64}.

\code{simd_div()} always returns a double, like \code{/}. On integers,
\code{simd_idiv()} and \code{simd_mod()} return \code{NA} (without a
warning) for a zero divisor, as base R does.

On doubles, \code{simd_idiv(x, y)} is the exact floor of the real quotient
\eqn{x / y} and \code{simd_mod(x, y)} the exact floored remainder
\code{x - simd_idiv(x, y) * y}, rounded once; it has the sign of \code{y}.
These are what base R computes in \code{long double}: \code{-2 \%/\% Inf} is
\code{-1}, \code{-2 \%\% Inf} is \code{Inf}, \code{x \%/\% 0} is
\code{x / 0} and \code{x \%\% 0} and \code{Inf \%\% y} are \code{NaN}. Where
base R's \code{long double} is not wide enough to be exact (on x86, or where
there is none, as on arm64 macOS) base R can differ in the last bits for
large quotients, and it warns \code{"probable complete loss of accuracy in
modulus"} above a platform-dependent quotient; \pkg{rsimd} stays exact and
never warns. A zero remainder is \code{0}, not \code{-0}, except that
\code{simd_mod(x, Inf)} returns \code{x} unchanged as base R does.

\rsimdSemantics

\rsimdAttributes

All implementations have their own kernels for these functions. Integer
\code{\%/\%} and \code{\%\%} divide in double precision, which is exact for
32-bit integers. Double \code{\%\%} with a quotient of \eqn{2^{52}} or more
falls back to scalar code for that element. On \code{"sse2"}, double
\code{\%/\%} and \code{\%\%} use the C library's \code{fma()} for the exact
remainder, which is slow on CPUs without a fused multiply-add. On 32-bit ARM
the double operations use the scalar \code{"none"} code.
}
\value{
A vector of the length of the longer operand: complex when an operand is
complex, integer when both operands are integer or logical (except for
\code{simd_div()}), otherwise double.
}
\seealso{
\link{Arithmetic}, \code{\link{simd_fma}}, \code{\link{simd_neg}},
\link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_add(1:5, 10L)
simd_mul(c(1.5, NA, 3), 2)
simd_idiv(c(-7, 7), 2)
simd_mod(c(-7, 5.5), 2)
simd_add(.Machine$integer.max, 1L)    # NA with a warning
simd_add_wrap(.Machine$integer.max, 2L)
simd_add(c(1, 2), 1i)
}
