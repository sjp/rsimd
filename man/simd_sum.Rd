\name{simd_sum}
\alias{simd_sum}
\alias{simd_prod}
\alias{simd_mean}
\alias{simd_prod_sums}
\alias{simd_prod_diffs}
\title{Sums, Products and Means of Vector Elements}
\usage{
simd_sum(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_prod(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_mean(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_prod_sums(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))

simd_prod_diffs(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE))
}
\arguments{
\item{x, y}{A double, integer, logical or complex vector, or for
\code{simd_sum()} a \pkg{bit64} \code{integer64} vector. \code{y} (for
\code{simd_prod_sums()} and \code{simd_prod_diffs()}) must have the length
of \code{x}, or either may have length 1.}
\item{na.rm}{\code{TRUE} or \code{FALSE}: remove missing values (\code{NA}
and \code{NaN}) first?}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} (and \code{y}) have no missing values, so they are not looked
for; see Details.}
}
\description{
The sum, product or arithmetic mean of the elements of a vector, like
\code{\link{sum}()}, \code{\link{prod}()} and \code{\link{mean}()} for a
single vector, and the product of the sums or differences of two vectors,
like \code{prod(x + y)} and \code{prod(x - y)} without the vector of sums,
computed with SIMD instructions.
}
\details{
Integer and logical vectors are summed exactly. Double vectors are summed in
the mode set with \code{\link{simd_precision}()}, so the result can differ
from base R's in the last bits.

\code{simd_sum()} of a complex vector sums the real and the imaginary
parts separately, each exactly as a double vector of that part would be
summed, and its missing values follow the rules below part by part:
\code{simd_sum(c(1+2i, NA))} has an \code{NA} real part and an imaginary
part of 2, as in base R. With \code{na.rm = TRUE} an element is removed
when either part is missing, as \code{is.na()} says.

\code{simd_prod()} of a double vector multiplies in double precision in
every mode, keeping several partial products. Base R multiplies in \code{long double} where the
platform has it, whose wider exponent range can avoid an intermediate
overflow or underflow: \code{prod(c(1e300, 1e300, 1e-300, 1e-300))} is 1
in base R, while \code{simd_prod()} gives \code{Inf} or 1 depending on the
implementation in use, which forms its partial products in its own order.

\code{simd_prod()} of a complex vector multiplies with base R's complex
multiplication \code{x * y} (see \code{\link{simd_mul}()}), in an order set
by the precision mode: \code{"fast"} keeps one partial product per vector
lane and multiplies them together at the end, so the result depends on the
implementation in the last bits; \code{"pairwise"} multiplies blocks of 128
elements in order and the blocks' products pairwise; \code{"compensated"}
multiplies in order in double-double arithmetic (about 106 bits), which is
usually as accurate as base R. The last two do not depend on the
implementation, and are not vectorised. Base R multiplies in \code{long
double} with the plain formula, without the recovery of infinities of
\code{x * y}, so results with infinite or overflowing parts can differ. A
missing element (either part \code{NA} or \code{NaN}) makes both parts of
the product \code{NA}, or \code{NaN} when there is no \code{NA}.

\code{simd_prod_sums(x, y)} is \code{simd_prod(x + y)}, and
\code{simd_prod_diffs(x, y)} is \code{simd_prod(x - y)}, with the sums or
differences added into the products as they are read. The rules of
\code{simd_prod()} apply, to the sums. Integer and logical operands are
read as doubles, so their sums never overflow: \code{simd_prod_sums(x, y)}
is \code{prod(as.double(x) + y)}, where \code{prod(x + y)} of integer
vectors would give \code{NA} with a warning when a sum leaves the integer
range. A pair whose sum is \code{NaN} is missing: when \code{x} or
\code{y} is \code{NA} or \code{NaN} there, and also for \code{Inf - Inf}.
Without \code{na.rm} the result is then \code{NA}, or \code{NaN} if no
operand was \code{NA}; with \code{na.rm = TRUE} the pair is left out, as
\code{prod(x + y, na.rm = TRUE)} leaves out its sum. A complex operand makes
the other complex (as \code{\link{as.complex}()}) and the result complex,
as \code{simd_prod(\link{simd_add}(x, y))}, in the precision mode;
integer64 operands are not taken. Operands must have the same length, or
one of them length 1 (see \code{\link{simd_add}}).

\code{simd_mean()} of a double vector is the sum divided by the number of
elements; in \code{"pairwise"} and \code{"compensated"} modes it then adds
the mean of the deviations from that first estimate (a second pass over the
data), as base R does. The mean of an integer or logical vector is its exact
sum divided once, in \code{long double} where available, as in base R.
The mean of a complex vector is \code{simd_sum()}'s result divided by the
number of elements, part by part (without base R's second pass in
\code{long double}), so a missing part makes that part of the mean
missing, as in base R.

\rsimdSemantics

\rsimdPrecision

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double versions use the scalar \code{"none"} code.

\code{simd_sum()} of an integer64 vector is an exact integer64, \code{NA}
with a warning when the total does not fit (see \link{rsimd-integer64}); \code{simd_prod()}
and \code{simd_mean()} do not take integer64 yet.
}
\value{
\code{simd_sum()}: a double for double \code{x}. For integer and logical
\code{x}, an integer, or a double if the sum is outside the integer range,
as in base R. For complex \code{x}, a complex number. The sum of an empty
vector is \code{0} (or \code{0L}, or \code{0+0i}).

\code{simd_prod()}, \code{simd_mean()}, \code{simd_prod_sums()} and
\code{simd_prod_diffs()}: a double, or for complex \code{x} a complex
number. The product of an empty vector is 1 (or \code{1+0i}) and its mean
\code{NaN} (or \code{NaN+NaNi}).
}
\seealso{
\code{\link{sum}}, \code{\link{prod}}, \code{\link{mean}},
\code{\link{simd_precision}}, \link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_sum(c(1.5, 2.5, NA), na.rm = TRUE)
simd_sum(1:10)
simd_sum(c(TRUE, FALSE, TRUE))
simd_sum(rep(.Machine$integer.max, 3L)) # a double, as base R
simd_sum(c(1+2i, 3-1i, NA), na.rm = TRUE)

simd_prod(1:10)
simd_prod(c(1i, 1i, 2))
simd_mean(c(1+2i, 3-1i))
simd_mean(c(1, 2, NaN, 4), na.rm = TRUE)

# The probability that none of independent events happens.
p <- c(0.1, 0.25, 0.05)
simd_prod_diffs(1, p)
simd_prod_sums(1:3, c(0.5, NA, 1), na.rm = TRUE)

old <- simd_precision("compensated")
simd_sum(c(1e16, 1, -1e16))
simd_precision(old)
}
