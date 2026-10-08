\name{rsimd-integer64}
\alias{rsimd-integer64}
\title{64-bit Integers (integer64)}
\description{
\pkg{rsimd} works on the 64-bit integers of the \pkg{bit64} package: a
double vector whose elements hold the bits of signed 64-bit integers and
which has the class \code{"integer64"}. \code{NA} is the smallest value,
\eqn{-2^{63}}, so values range over
\eqn{\pm(2^{63} - 1) = \pm 9223372036854775807}. Such vectors are recognised
by their class, and results get that class. \pkg{bit64} is only suggested,
but any function that returns an integer64 result, including
\code{\link{simd_as_integer64}()}, loads it, because without its
methods registered an integer64 vector sorts, prints and compares as the
doubles its bits make. If \pkg{bit64} is not installed, those functions
give the error \dQuote{integer64 results need package 'bit64'}. Functions
whose results are not integer64, such as \code{\link{simd_as_double}()},
comparisons and \code{\link{simd_which_max}()}, still work on integer64
data without it.
}
\details{
\strong{Functions that take integer64.} Results are integer64 unless noted.
\describe{
\item{Reductions}{\code{\link{simd_sum}}, \code{\link{simd_min}},
\code{\link{simd_max}}, \code{\link{simd_range}}, \code{\link{simd_cumsum}},
\code{\link{simd_cummin}}, \code{\link{simd_cummax}},
\code{\link{simd_max_abs}} and \code{\link{simd_min_abs}};
\code{\link{simd_which_min}}, \code{\link{simd_which_max}},
\code{\link{simd_which_min_abs}} and \code{\link{simd_which_max_abs}}
(integer, \code{NA} ignored); \code{\link{simd_any}} and \code{\link{simd_all}}
(logical; non-zero is \code{TRUE}); \code{\link{simd_any_na}},
\code{\link{simd_count_na}} (double) and \code{\link{simd_which_na}}.}
\item{Arithmetic}{\code{\link{simd_add}}, \code{\link{simd_sub}},
\code{\link{simd_mul}}, their \code{_wrap} variants, \code{\link{simd_idiv}},
\code{\link{simd_mod}}, \code{\link{simd_neg}}, \code{\link{simd_abs}},
\code{\link{simd_sign}} (\eqn{-1}, 0 or 1), \code{\link{simd_pmin}},
\code{\link{simd_pmax}}, \code{\link{simd_pmin_num}},
\code{\link{simd_pmax_num}}, \code{\link{simd_clamp}},
\code{\link{simd_mul_add}} and \code{\link{simd_add_mul}};
\code{\link{simd_div}} returns a double: it converts both operands to double
and then divides, so a quotient can differ by 1 ULP from \pkg{bit64}'s
\code{/}, which divides in \code{long double} and rounds once, and a zero
divisor gives \code{Inf}, \code{-Inf} or \code{NaN} without a warning, as
\code{/} does for doubles (\pkg{bit64}'s \code{/} gives \code{NA} with a
warning).}
\item{Predicates and comparisons}{\code{\link{simd_is_na}},
\code{\link{simd_is_nan}} (always \code{FALSE}), \code{\link{simd_is_finite}}
(not \code{NA}), \code{\link{simd_is_infinite}} (always \code{FALSE}),
\code{\link{simd_is_negative}}, \code{\link{simd_is_zero}},
\code{\link{simd_is_normal}} (not 0), \code{\link{simd_is_subnormal}}
(always \code{FALSE}), \code{\link{simd_is_whole}} (not \code{NA}),
\code{\link{simd_is_even}}, \code{\link{simd_is_odd}},
\code{\link{simd_is_pow2}}, each with its \code{_any} and \code{_all}
forms; \code{\link{simd_eq}}, \code{\link{simd_ne}}, \code{\link{simd_lt}},
\code{\link{simd_le}}, \code{\link{simd_gt}} and \code{\link{simd_ge}}
(logical); \code{\link{simd_hamming}} and \code{\link{simd_hamming_bits}}
(double).}
\item{Bit operations}{\code{\link{simd_bit_and}}, \code{\link{simd_bit_or}},
\code{\link{simd_bit_xor}}, \code{\link{simd_bit_not}}, \code{\link{simd_shl}},
\code{\link{simd_shr}}, \code{\link{simd_sar}}, \code{\link{simd_rotl}},
\code{\link{simd_rotr}}; \code{\link{simd_popcount}},
\code{\link{simd_lzcnt}} and \code{\link{simd_tzcnt}} (integer counts, 0 to
64); \code{\link{simd_popcount_total}} (double).}
\item{Conversions}{\code{\link{simd_as_double}}, \code{\link{simd_as_integer}},
\code{\link{simd_as_logical}}, \code{\link{simd_as_raw}} from integer64, and
\code{\link{simd_as_integer64}} from double, integer, logical, raw and
integer64.}
}
Other functions reject integer64 with an error. \code{simd_prod} and
\code{simd_mean} are not available for integer64 yet; transcendental and
statistical functions (\code{simd_exp}, \code{simd_sqrt}, \code{simd_var},
\code{simd_dot} \ldots), \code{simd_fma}, \code{simd_lerp},
\code{simd_copysign}, the rounding functions, \code{simd_cumprod}, the
three-valued logic of \code{simd_and} \ldots and the complex functions do
not take it.

\strong{Mixed operands.} An integer or logical operand combined with an
integer64 one is converted to integer64 (exactly; \code{NA} stays
\code{NA}). A double operand makes the integer64 operands convert to double:
the result is double, with the warning \code{"integer64 coerced to double"}
(\pkg{bit64} instead converts the double to integer64; \pkg{rsimd} keeps the
double's fraction but loses integer64 precision beyond \eqn{2^{53}}, so
convert explicitly with \code{\link{simd_as_integer64}()} when that matters).
Raw and complex operands cannot be combined with integer64.

\strong{Overflow.} Checked arithmetic gives \code{NA} with \pkg{bit64}'s
warning \code{"NAs produced by integer64 overflow"} (once per call) for
results outside \eqn{\pm(2^{63} - 1)}, including \eqn{-2^{63}}, the bit
pattern of \code{NA}. The \code{_wrap} variants wrap around modulo
\eqn{2^{64}} instead, never warn, and give \code{NA} for \code{NA} operands;
a result that wraps to \eqn{-2^{63}} reads as \code{NA}.
\code{simd_sum} is exact: it gives \code{NA} with that warning only when the
exact total is out of range, whatever the order of the elements, so every
implementation agrees. \pkg{bit64}'s \code{sum} checks every running total
instead, so it returns \code{NA} for \code{c(MAX, 1, -1)}, where
\code{simd_sum} returns \code{MAX}. \code{simd_cumsum} is \code{NA} from the
first running total that is out of range on, with that warning. A missing
value makes the result \code{NA} without a warning.

\strong{Division.} \code{simd_idiv} and \code{simd_mod} round the quotient
down (\code{-7 \%/\% 2} is \eqn{-4} and \code{-7 \%\% 2} is 1, the remainder
taking the sign of the divisor), as \pkg{bit64} does. A zero divisor gives
\code{NA} with the warning \code{"NAs produced due to division by zero"}
(once per call).

\strong{Empty input.} \code{simd_min} and \code{simd_max} of no values
(after \code{na.rm}) return \eqn{+(2^{63} - 1)} and \eqn{-(2^{63} - 1)} with
\pkg{bit64}'s warnings, such as \code{"no non-NA value, returning the highest
possible integer64 value +9223372036854775807"}; \code{simd_range} returns
both with one warning. The sum of no values is 0.

\strong{Conversions.} \code{simd_as_double} is exact up to \eqn{2^{53}} in
magnitude and rounds to the nearest double (ties to even) beyond, with
\pkg{bit64}'s warning \code{"integer precision lost while converting to
double"} once when any value has a magnitude of \eqn{2^{53}} or more.
\code{simd_as_integer} gives \code{NA} with the warning \code{"NAs produced
by integer overflow"} for values outside the integer range in
\code{"checked"} mode, clamps in \code{"saturating"} mode and keeps the low
32 bits in \code{"truncating"} mode. \code{simd_as_raw} in
\code{"checked"} mode gives 0 with the single warning
\code{"out-of-range values treated as 0 in coercion to raw"}, as
\code{as.raw()} of an integer does; \pkg{bit64}'s \code{as.raw()} converts
to integer first, so it also warns \code{"NAs produced by integer
overflow"} for values outside the integer range.
\code{simd_as_integer64} truncates doubles toward zero; see \code{\link{simd_as_integer64}} for its modes.

\strong{Bit operations} act on the 64-bit patterns. Shift counts must be in
0 to 63 (any other count gives \code{NA} everywhere); rotate counts are taken
modulo 64. A result with the bit pattern of \code{NA} (only the sign bit set,
as in \code{simd_shl(x, 63)} of an odd \code{x}) reads as \code{NA}.

\strong{Missing values.} \code{na_check = FALSE} promises there is no
\code{NA} and skips looking for it. A broken promise counts an
\code{integer64} \code{NA} as \eqn{-2^{63}}. An integer or logical
\code{NA} operand counts as \eqn{-2^{31}} in arithmetic
(\code{simd_add(as.integer64(1), NA_integer_, na_check = FALSE)} is
\eqn{-2147483647}), but the bit operations always convert it to the
\code{integer64} \code{NA}.

\strong{Attributes.} \rsimdAttributes Note that a class vector with more
than one element (such as that of a \code{simd_vec}) is recognised as
integer64 when it includes \code{"integer64"}, but only a class of exactly
\code{"integer64"} is copied to the result.

\strong{Implementations.} The following table shows how each operation runs
on each implementation: \emph{native} uses the tier's own instructions;
\emph{SIMDe} uses SIMDe's emulation of an instruction the tier lacks (64-bit
compares on \code{"sse2"}, which needs SSE4.2 for them); \emph{partial}
builds the operation from other vector instructions (32-bit partial
products, shifts and the sign, a bit-counting sequence); \emph{scalar} runs
a plain loop, also used by \code{"none"} everywhere. Only \code{"neon"} has
been timed when choosing; the choices for \code{"sse2"}, \code{"avx2"},
\code{"avx512"} and \code{"sve"} follow from the instruction sequences and
are not measured.
\tabular{llllll}{
\strong{operation} \tab \strong{sse2} \tab \strong{avx2} \tab \strong{avx512} \tab \strong{neon} \tab \strong{sve, sve2} \cr
sum \tab native, SIMDe NA test \tab native \tab native \tab native \tab native \cr
add, sub \tab native, SIMDe overflow test \tab native \tab native \tab native \tab native \cr
min, max, their magnitudes, pmin, pmax, clamp, comparisons \tab SIMDe \tab native \tab native \tab native \tab native \cr
missing values, predicates, any, all \tab native \tab native \tab native \tab native \tab native \cr
mul (checked) \tab partial \tab native \tab native \tab native \tab native \cr
mul_wrap \tab partial \tab partial \tab native \tab partial \tab native \cr
neg, abs, sign \tab native \tab native \tab native \tab native \tab native \cr
idiv, mod, mul_add, add_mul \tab scalar \tab scalar \tab scalar \tab scalar \tab scalar \cr
cumsum, cummin, cummax \tab scalar \tab scalar \tab scalar \tab scalar \tab scalar \cr
and, or, xor, not, shl, shr, rotl, rotr \tab native \tab native \tab native \tab native \tab native \cr
sar \tab partial \tab partial \tab partial \tab partial \tab native \cr
popcount, hamming_bits \tab partial \tab partial \tab partial \tab native \tab native \cr
is_normal, is_even, is_odd, is_pow2, hamming \tab SIMDe \tab native \tab native \tab native \tab native \cr
lzcnt, tzcnt \tab partial \tab partial \tab partial \tab partial \tab native \cr
integer64 to double \tab partial \tab partial \tab native \tab native \tab native \cr
double to integer64 \tab scalar \tab scalar \tab scalar \tab native* \tab native* \cr
integer64 to integer and logical, integer to integer64 \tab native \tab native \tab native \tab native \tab native \cr
integer64 to and from raw \tab scalar \tab scalar \tab scalar \tab scalar \tab scalar
}
* \code{"checked"} and \code{"saturating"} modes; \code{"truncating"} is
scalar. Checked multiplication takes the vector path for each vector whose
operands are all within \eqn{\pm 2^{31}} (and, on SVE, always); the
elements of a vector with a larger operand are multiplied in scalar
code, and after eight such vectors in a row the rest of the block being
processed (up to \eqn{2^{20}} elements) is. \emph{SIMDe NA test}:
the sum's own additions are native SSE2, and only the test for \code{NA}
(skipped with \code{na_check = FALSE}) uses SIMDe's 64-bit compare.
\emph{SIMDe overflow test}: the additions and subtractions are native,
and only the overflow test of the checked forms uses SIMDe's compares. On 32-bit ARM the \code{"neon"} kernels use NEON's
64-bit lanes but convert to and from double in scalar code.
}
\seealso{
\code{\link{simd_as_integer64}}, \code{\link{simd_sum}},
\code{\link{simd_add}}, \link{rsimd_options}; the \pkg{bit64} package.
}
\examples{
if (requireNamespace("bit64", quietly = TRUE)) {
  x <- simd_as_integer64(c(2^40, 7, NA, -3))
  print(simd_sum(x, na.rm = TRUE))
  print(simd_which_max(x))
  print(simd_gt(simd_add(x, 1L), 7L))
  big <- simd_as_integer64(2^62)
  print(simd_is_na(suppressWarnings(simd_mul(big, 2L)))) # overflow
  print(simd_lt(simd_mul_wrap(big, 3L), 0L)) # wraps around
  print(simd_sum(bit64::as.integer64(c("9223372036854775807", "1", "-1"))))
  print(simd_idiv(bit64::as.integer64(-7), 2L))
}
}
