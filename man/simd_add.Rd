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
simd_add(x, y, na_check = NULL)

simd_sub(x, y, na_check = NULL)

simd_mul(x, y, na_check = NULL)

simd_div(x, y, na_check = NULL)

simd_idiv(x, y, na_check = NULL)

simd_mod(x, y, na_check = NULL)

simd_add_wrap(x, y, na_check = NULL)

simd_sub_wrap(x, y, na_check = NULL)

simd_mul_wrap(x, y, na_check = NULL)
}
\arguments{
\item{x, y}{Double, integer, logical or \pkg{bit64} \code{integer64}
vectors; not double for the \code{_wrap} functions. \code{simd_add()},
\code{simd_sub()}, \code{simd_mul()} and \code{simd_div()} also take
complex vectors.}
\item{na_check}{\code{TRUE} or \code{FALSE}, or \code{NULL} (the default) for option \code{rsimd.na_check} (\code{\link{rsimd_options}}). \code{FALSE} promises that
\code{x} and \code{y} have no missing values, so they are not looked for; see
Details.}
}
\description{
\code{x + y}, \code{x - y}, \code{x * y}, \code{x / y}, \code{x \%/\% y} and
\code{x \%\% y} computed with SIMD instructions, and integer addition,
subtraction and multiplication that wrap around on overflow.
}
\details{
\rsimdBroadcast

\code{simd_add()}, \code{simd_sub()}, \code{simd_mul()} and
\code{simd_div()} take complex operands: if either operand is complex, the
other (double, integer or logical) is converted with \code{as.complex()},
so an integer or logical \code{NA} becomes \code{NA} with a zero imaginary
part. For \code{simd_add()} and \code{simd_sub()} the real and imaginary
parts are then added or subtracted separately, each following the rules for
doubles below, so a missing part makes only that part of the result
missing. Complex operands do not combine with \code{integer64} ones.

Complex \code{simd_mul()} and \code{simd_div()} give exactly base R's
results, special values and the choice between \code{NA} and \code{NaN}
included (\code{na_check} does not apply). Base R multiplies and divides
complex numbers with the C compiler's complex arithmetic, so its results
depend on the compiler that built R: whether it fuses a product into the
following addition (a fused multiply-add), in which order, and which
division algorithm its runtime library has. When \pkg{rsimd} loads, base R
computes products, quotients and a \code{cumprod()} of a fixed set of
inputs, and \pkg{rsimd} picks the variant of its kernels that reproduces
all of them; \code{\link{simd_complex_variants}()} names the
choice. A product's parts are \eqn{ac - bd} and \eqn{ad + bc} with one
product of each fused or not; division follows the scaled Smith's method of
GCC's runtime (GCC 12 and later). Elements whose result has a missing part
(infinite operands, a zero divisor, missing values) are recomputed with
base R's own operator. Where no variant matches (for example division on
macOS, whose runtime uses another algorithm), every element is computed with
base R's operator, which is exact but no faster than base R.

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
\code{x / 0} and \code{x \%\% 0} and \code{Inf \%\% y} are \code{NaN}. But
base R's \code{\%\%} is only exact for quotients up to the precision of its
\code{long double}: above about \eqn{2^{52}} where there is none (as on
arm64 macOS), \eqn{2^{63}} on x86 and \eqn{2^{112}} with a 128-bit
\code{long double}, it warns \code{"probable complete loss of accuracy in
modulus"} and its result can be wrong in every bit (\code{1e300 \%\% 7} is
often \code{0}, not \code{1}); and where \eqn{x / y} overflows, as with a
subnormal divisor, it returns \code{NaN} without a warning (\code{1 \%\%
5e-324} is \code{NaN}, not \code{0}). \pkg{rsimd} stays exact and never
warns. A zero remainder is \code{0}, not \code{-0}, except that
\code{simd_mod(x, Inf)} returns \code{x} unchanged as base R does.

\rsimdSemantics

\rsimdAttributes

All implementations have their own kernels for these functions. Integer
\code{\%/\%} and \code{\%\%} divide in double precision, which is exact for
32-bit integers. Double \code{\%\%} with a quotient of \eqn{2^{52}} or more
falls back to scalar code for that element. On \code{"sse2"}, double
\code{\%/\%} and \code{\%\%} use the C library's \code{fma()} for the exact
remainder, which is slow on CPUs without a fused multiply-add, and so are
complex products where R's compiler fuses them. On 32-bit ARM the double and
complex operations use the scalar \code{"none"} code.

A length-one double operand is fastest as \code{y}: \code{simd_sub(x, 1)}
keeps it in a register, while \code{simd_sub(1, x)} reads it from memory
for every vector and can take somewhat longer (about as long as
\code{1 - x}). An integer or logical vector combined with a double operand
is converted to double element by element, which is slower than a double
vector.
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
simd_mod(c(1e300, 1), c(7, 5e-324))   # 1 and 0, exactly; base R gives 0 or NaN
simd_add(.Machine$integer.max, 1L)    # NA with a warning
simd_add_wrap(.Machine$integer.max, 2L)
simd_add(c(1, 2), 1i)
z <- complex(real = c(1, 0, Inf), imaginary = c(2, 0, 0))
identical(simd_mul(z, 1 - 1i), z * (1 - 1i))
identical(simd_div(1, z), 1 / z)
simd_complex_variants()
}
