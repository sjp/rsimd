\name{simd_ilogb}
\alias{simd_ilogb}
\alias{simd_scaleb}
\alias{simd_nextafter}
\alias{simd_next_up}
\alias{simd_next_down}
\alias{simd_remainder}
\title{Binary Exponents, Scaling by Powers of Two, Neighbouring Doubles and
IEEE Remainders}
\usage{
simd_ilogb(x)

simd_scaleb(x, n)

simd_nextafter(x, y)

simd_next_up(x)

simd_next_down(x)

simd_remainder(x, y)
}
\arguments{
\item{x, y}{Double, integer or logical vectors.}
\item{n}{An integer or logical vector, or a double vector of whole numbers
in the integer range (\code{NA} allowed).}
}
\description{
Exact operations on the representation of doubles, computed with SIMD
instructions: the binary exponent of \code{x}, \code{x * 2^n}, the
neighbouring doubles of \code{x}, and the remainder of \code{x / y} with
the quotient rounded to the nearest integer. They correspond to C's
\code{ilogb()}, \code{scalbn()}, \code{nextafter()}, \code{nextup()},
\code{nextdown()} and \code{remainder()}.
}
\details{
Every implementation gives exactly the results of the C functions, so the
results are identical on every implementation.

\code{simd_ilogb(x)} is the integer \eqn{e} with
\eqn{2^e \le |x| < 2^{e+1}}{2^e <= |x| < 2^(e+1)}, also for subnormal
numbers (\code{simd_ilogb(5e-324)} is \code{-1074}). It is \code{NA} for
zero, infinities, \code{NaN} and \code{NA}, where C returns special
integers.

\code{simd_scaleb(x, n)} is \code{x * 2^n} computed exactly and rounded
once, so a result in the subnormal range is rounded correctly, a result
too large is \code{Inf} and one too small is a zero (with the sign of
\code{x}), even for \code{n} far outside the range of exponents
(\code{simd_scaleb(1, .Machine$integer.max)} is \code{Inf}). An \code{NA}
(or \code{NaN}) in \code{n} gives \code{NA}. A double \code{n} that is not
a whole number in the integer range is an error.

\code{simd_next_up(x)} is the smallest double greater than \code{x} and
\code{simd_next_down(x)} the largest double smaller than \code{x} (IEEE
754's \code{nextUp} and \code{nextDown}): \code{simd_next_up(0)} and
\code{simd_next_up(-0)} are the smallest subnormal \code{5e-324},
\code{simd_next_up(-5e-324)} is \code{-0}, \code{simd_next_up(Inf)} is
\code{Inf} and \code{simd_next_up(-Inf)} is
\code{-.Machine$double.xmax}. \code{simd_nextafter(x, y)} is the next double
after \code{x} in the direction of \code{y}, or \code{y} when the two are
equal (so \code{simd_nextafter(0, -0)} is \code{-0}).

\code{simd_remainder(x, y)} is \code{x - k * y} computed exactly, where
\code{k} is \code{x / y} rounded to the nearest integer, ties to even. The
result is at most \code{|y| / 2} in magnitude and can be negative where
\code{x} and \code{y} are positive: \code{simd_remainder(5, 3)} is \code{-1}
where \code{5 \%\% 3} is \code{2}, because \code{\%\%} rounds the quotient
down (see \code{\link{simd_mod}}). A zero result has the sign of \code{x}.
\code{simd_remainder(x, Inf)} is \code{x}; \code{simd_remainder(x, 0)} and
\code{simd_remainder(Inf, y)} are \code{NaN} and warn
\code{"NaNs produced"}, unlike \code{\%\%}, which is silent.

None of these functions depends on option \code{rsimd.math_accuracy}, and
none takes complex input.

\rsimdBroadcast

\rsimdMath

\rsimdAttributes
}
\value{
An integer vector for \code{simd_ilogb()}, otherwise a double vector, of
the length of the longer operand.
}
\seealso{
\code{\link{simd_mod}}, \code{\link{simd_rootn}}, \code{\link{simd_log2}},
\code{\link{simd_is_subnormal}}.
}
\examples{
simd_ilogb(c(1, 3, 1024, 0.1, 5e-324, 0))
simd_scaleb(3, c(-2L, 10L, 2000L))
simd_scaleb(.Machine$double.xmax, -2098)
simd_next_up(1) - 1
simd_next_down(c(0, Inf))
simd_nextafter(1, c(0, 2)) - 1
simd_remainder(c(5, 7, -5), 2)
c(5, 7, -5) \%\% 2
}
