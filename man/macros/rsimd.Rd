% Shared wording for rsimd's Rd files. R loads macros from man/macros/
% automatically at install time; use these instead of copying the text.

\newcommand{\rsimdAttributes}{Inputs with a class other than \code{"integer64"} and \code{\link{simd_vec}}, such as dates, times, \code{difftime} or factors, are errors. Results are bare vectors: names, dimensions and other attributes of the inputs are dropped. The only exceptions are the \code{"integer64"} class from \pkg{bit64}, which is preserved, and \code{\link{simd_vec}} operands: the function runs on their pinned implementation, and a result that is a vector of values (not a reduction, predicate or comparison) is a \code{simd_vec} again.}

\newcommand{\rsimdSemantics}{\strong{Missing values.} \code{NA} and \code{NaN} are kept
apart. With \code{na.rm = FALSE} a reduction returns \code{NA} if the input
contains an \code{NA}, and otherwise \code{NaN} if it contains a \code{NaN}
(\code{var} and \code{sd} return \code{NA} for either, as base R does).
When \code{na.rm = FALSE} and the input contains both \code{NA} and
\code{NaN}, \pkg{rsimd} always returns \code{NA}; base R's result depends on
the order of elements. With \code{na.rm = TRUE} both are removed (as
\code{is.na(NaN)} is \code{TRUE}). If nothing remains, sums are 0, products
1, means \code{NaN}, \code{var} and \code{sd} \code{NA} (as for fewer than
two values), and \code{min} and \code{max} \code{Inf} and \code{-Inf} with
base R's warning. \code{any} is \code{TRUE} if any element is \code{TRUE},
otherwise \code{NA} if any is \code{NA}, otherwise \code{FALSE}; \code{all}
is \code{FALSE} if any element is \code{FALSE}, otherwise \code{NA} if any is
\code{NA}, otherwise \code{TRUE}. These may stop reading at the first element
that decides the answer.

\code{na_check = FALSE} (default: option \code{rsimd.na_check}) promises
that the input has no missing values and skips looking for them. If the
promise is broken, the result is unspecified (but never a crash): an integer
\code{NA} is counted as \code{-2147483648} and a double \code{NA} may come
back as \code{NaN}. \code{na.rm = TRUE} always looks for missing values.

Elementwise operations on doubles return \code{NA} wherever an operand is
\code{NA}, and \code{NaN} where IEEE arithmetic gives one (base R may give
\code{NaN} for \code{NaN + NA}, depending on the CPU). \code{pmin},
\code{pmax} and \code{clamp} instead return one of their operands, as base
R does. On integers an \code{NA} operand gives \code{NA}.

\strong{Integer overflow.} Integer results outside
\code{[-.Machine$integer.max, .Machine$integer.max]} are \code{NA}, with the
warning \code{"NAs produced by integer overflow"} once per call, as in base
R. The \code{_wrap} variants wrap around in two's complement instead and
never warn, but still return \code{NA} for \code{NA} operands. Beware that
\code{-2147483648} is the bit pattern of \code{NA_integer_}, so a result that
wraps to it reads as \code{NA}: \code{simd_add_wrap(.Machine$integer.max, 1L)}
is \code{NA}. Integer sums are exact, and are returned as a double when the
total does not fit in an integer, as in base R.}

\newcommand{\rsimdPrecision}{Floating-point sums (in \code{sum}, \code{mean},
\code{sum_sq}, \code{sum_abs}, \code{dot}, \code{norm}, \code{dist},
\code{cosine}, \code{var}, \code{sd}, \code{cumsum} and the denominators of
\code{softmax} and \code{log_softmax}) are accumulated in the mode set with
\code{\link{simd_precision}()}, which also sets the order of the
multiplications in \code{prod} of a complex vector (see
\code{\link{simd_prod}}); other operations are not affected.
\describe{
\item{\code{"fast"}}{(the default) keeps several partial sums in the lanes of
the vector registers and adds them at the end. It is the fastest mode. The
result can differ from base R's in the last bits, and it depends on the
implementation in use (\code{\link{simd_current}()}) and on how the input is
split into blocks internally.}
\item{\code{"pairwise"}}{adds blocks of 128 elements as in \code{"fast"} and
combines the block sums pairwise, so rounding errors grow with
\eqn{\log n} rather than \eqn{n}. The result does not depend on how the
input is split internally; implementations can differ in about the last
bit.}
\item{\code{"compensated"}}{uses Neumaier's compensated summation, which
carries the rounding error of every addition along. Its error is about one
rounding of the result unless the sum cancels heavily, implementations agree
to within one unit in the last place, and it needs no \code{long double}, so
it gives the same results on every platform. It is the closest to base R:
\code{sum(c(1e16, 1, -1e16))} is 1 in this mode (and in base R where
\code{long double} is available), while the other modes can give 0.}
}
Base R's \code{sum()} accumulates doubles in \code{long double} where the
platform has it, so no mode promises bit-identical results to base R.

In \code{"pairwise"} and \code{"compensated"} modes, \code{mean} refines the
first estimate with a second pass over the data, as base R does; \code{var}
and \code{sd} always take two passes (mean, then squared deviations).
Sums of products (\code{sum_sq}, \code{dot}, \code{dist}, \code{var},
\code{sd}) use fused multiply-adds in \code{"fast"} mode on implementations
that have them (all but \code{"sse2"} and \code{"none"}). \code{cumsum}
works like \code{"fast"} in \code{"pairwise"} mode: it adds the elements of
each vector register in a tree before adding the running total, so it can
differ from base R's sequential running sum in the last bits. In
\code{"compensated"} mode it is a sequential compensated running sum,
identical on every implementation. \code{norm} and \code{dist} are not
protected against overflow of the sum of squares.}

\newcommand{\rsimdLengths}{\strong{Lengths.} The operands must have the same
length, or length one: a length-one operand is used for every element
(broadcast). Any other combination is an error; there is no general
recycling. A zero-length operand with length-one operands gives a
zero-length result.}

\newcommand{\rsimdBroadcast}{\rsimdLengths

\strong{Types.} Logical operands count as integers. Integer and logical
operands give an integer result where base R's operator does; any double
operand makes the result double. Integers mixed with doubles are converted
as they are read. Raw vectors are rejected with base R's message. Which
functions also take complex or \pkg{bit64} \code{integer64} operands, and
what those give, is stated with each function.}

\newcommand{\rsimdMath}{\strong{Implementations and accuracy.} The SIMD
implementations use the vectorised functions of the bundled SLEEF library:
its 1-ULP (\code{_u10}) variants, and its 0.5-ULP (\code{_u05}) variants for
\code{sinpi}, \code{cospi} and \code{hypot}. Results are within 1 ULP (unit
in the last place) of the exact value (0.5 ULP for the \code{_u05}
functions), unless the function's own description gives another bound,
and so may differ from base R in the last bit. Option
\code{rsimd.math_accuracy = "fast"} trades accuracy for speed with SLEEF's
3.5-ULP variants: see \code{\link{simd_math_accuracy}()} for the bounds of
each function in both modes. The \code{"none"}
implementation calls the C math library, as base R does, and is the
reference the others are tested against (within 2 ULP). An implementation
built without SLEEF (see \code{\link{simd_available}}) uses the next lower
one that has it, ultimately \code{"none"}.

\strong{Missing values and warnings.} A \code{NaN} input gives that
\code{NaN} back (so \code{NA} stays \code{NA}), and functions of two
arguments give \code{NA} if either is \code{NA} and otherwise \code{NaN} if
either is \code{NaN}, as base R's math functions do. When a result is
\code{NaN} although its inputs were not (\code{simd_log(-1)},
\code{simd_sin(Inf)}), the function warns \code{"NaNs produced"} once per
call. Overflow to \code{Inf} and underflow to zero (subnormal results are
produced, never flushed to zero) are not warned about. Integer and logical
input is converted to double, an integer \code{NA} becoming
\code{NA_real_}. Which functions also take complex input is stated with
each function; none takes \pkg{bit64} \code{integer64} input.

\strong{Signed zeros.} Odd functions return \code{-0} for \code{-0}, as C
does and base R did before version 4.6.0; since 4.6.0 base R returns
\code{0} for \code{sin}, \code{tan}, \code{asin}, \code{atan}, \code{sinh},
\code{tanh}, \code{asinh}, \code{atanh}, \code{expm1} and \code{log1p} of
\code{-0}.}
