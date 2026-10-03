% Shared wording for rsimd's Rd files. R loads macros from man/macros/
% automatically at install time; use these instead of copying the text.

\newcommand{\rsimdAttributes}{Results are bare vectors: names, dimensions and other attributes of the inputs are dropped. The only exceptions are the \code{"integer64"} class from \pkg{bit64} and the \code{"simd_vec"} class, which are preserved.}

\newcommand{\rsimdSemantics}{\strong{Missing values.} \code{NA} and \code{NaN} are kept
apart. With \code{na.rm = FALSE} a reduction returns \code{NA} if the input
contains an \code{NA}, and otherwise \code{NaN} if it contains a \code{NaN}.
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
\code{NA}, and \code{NaN} where IEEE arithmetic gives one. On integers an
\code{NA} operand gives \code{NA}.

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
\code{\link{simd_precision}()}; other operations are not affected.
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
\code{cumsum} is sequential, so \code{"pairwise"} works like
\code{"compensated"} there. \code{norm} and \code{dist} are not protected
against overflow of the sum of squares.}
