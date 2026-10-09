\name{simd_vec}
\alias{simd_vec}
\alias{as_simd_vec}
\alias{is_simd_vec}
\alias{simd_impl}
\alias{simd_impl<-}
\alias{simd_na_free}
\alias{simd_unwrap}
\title{SIMD Vectors with a Pinned Implementation}
\usage{
simd_vec(x, impl = NULL, scan_na = FALSE)

as_simd_vec(x, ...)

is_simd_vec(x)

simd_impl(x)
simd_impl(x) <- value

simd_na_free(x)

simd_unwrap(x)
}
\arguments{
\item{x}{For \code{simd_vec()} and \code{as_simd_vec()}, a double,
  integer, logical, raw or complex vector, or a \pkg{bit64}
  \code{integer64} vector. For the accessors, a \code{simd_vec}. For
  \code{simd_unwrap()}, either.}
\item{impl, value}{\code{NULL} (follow the global setting of
  \code{\link{simd_use}}) or the name of an implementation in
  \code{\link{simd_available}()} to pin the object to. \code{"auto"} is
  not allowed: a pin names a concrete implementation.}
\item{scan_na}{If \code{TRUE}, scan \code{x} once for \code{NA} and
  \code{NaN} and record the answer in the object's NA-free flag. This
  is unrelated to the \code{na_check} argument of the \code{simd_}
  functions, which says whether one call looks for missing values
  (\code{\link{rsimd_options}}).}
\item{...}{Passed to \code{simd_vec()}.}
}
\description{
A \code{simd_vec} is an ordinary atomic vector with class
\code{"simd_vec"} on which R's arithmetic, comparison and logical
operators, the \code{Math} and \code{Summary} group generics and
\code{mean()} run the \pkg{rsimd} kernels. It can be pinned to one SIMD
implementation and remembers when it is known to contain no missing
values. The wrapper is optional: every \code{simd_} function takes and
returns plain vectors.
}
\details{
\strong{Structure.} A \code{simd_vec} is the data vector with attributes
\describe{
\item{\code{class}}{\code{"simd_vec"}, or \code{c("simd_vec",
  "integer64")} for \code{integer64} data, so that \pkg{bit64}'s methods
  are reached for what \pkg{rsimd} does not implement.}
\item{\code{rsimd_impl}}{The pinned implementation, absent when the
  object is not pinned.}
\item{\code{rsimd_na_token}}{An external pointer that holds the
  NA-free flag and ties it to this object (see \dQuote{Where the flag
  holds} below), absent when the flag is unknown. The flag is
  \code{TRUE} if the data is known to contain no \code{NA} or
  \code{NaN} and \code{FALSE} if it is known to contain one; read it
  with \code{simd_na_free()}.}
}
\code{simd_vec()} drops every other attribute of \code{x} (names,
dimensions, dimnames). Given a \code{simd_vec}, it replaces the pin only
when \code{impl} is supplied and the flag only when \code{scan_na} is
\code{TRUE}. Character vectors, lists and data with any other class (dates,
times, \code{difftime}, factors, data frames) are errors, and so is
\code{NULL}. Making an
\code{integer64} \code{simd_vec} loads \pkg{bit64} (see
\link{rsimd-integer64}).
\code{simd_na_free()} returns the flag (\code{TRUE}, \code{FALSE} or
\code{NULL}); there is no setter, so the flag can only come from a check or
from the rules below. Raw data has no missing values, so a raw
\code{simd_vec} is always flagged \code{TRUE}.

\strong{Implementation pins.} Any \pkg{rsimd} function, called directly
or through an operator, runs on the implementation its \code{simd_vec}
operands are pinned to, whatever \code{\link{simd_use}()} or
\code{\link{simd_with_impl}()} selected:
\itemize{
\item one pinned operand, or several pinned to the same implementation:
  that implementation, and the result is pinned to it too;
\item operands pinned to different implementations: an error,
  \dQuote{operands pinned to different implementations ('avx2' vs 'sse2');
  unpin one with simd_impl(x) <- NULL};
\item no pinned operand: the global selection, and the result is not
  pinned.
}
A pin to an implementation that is not available on this machine (an
object saved on another one) is an error when the object is used; it
never falls back silently. Unpin it with \code{simd_impl(x) <- NULL}.

\strong{An invalid \code{rsimd.impl} option.} While the option
\code{rsimd.impl} holds a value that cannot be selected (see
\code{\link{simd_use}()}), nearly every operation on a \code{simd_vec}
fails with that error: printing, \code{format()}, \code{str()},
\code{x[i]}, \code{c()}, conversions, \code{simd_na_free()} and
\code{simd_unwrap()} as well as the operators. \code{length()},
\code{is_simd_vec()}, \code{simd_impl()} and \code{unclass()} still
work. Repair the option with \code{simd_use()} or
\code{options(rsimd.impl = "auto")}.

\strong{Results.} When any operand of an \pkg{rsimd} function is a
\code{simd_vec}, a result that is a vector of values (arithmetic,
elementary functions, rounding, bitwise operations, conversions,
cumulative functions, \code{simd_softmax()}, \code{simd_conj()} and so on)
is a \code{simd_vec}. Reductions (\code{sum}, \code{range},
\code{simd_which_max} \dots), predicates, comparisons and the logical
operators \code{&}, \code{|} and \code{!} return plain values.

\strong{The NA-free flag.} When every operand is flagged NA-free, the
functions skip their missing-value checks (as with \code{na_check =
FALSE}), for reductions and elementwise operations alike. The flag is only
set on objects a function returns, since R cannot change an argument in
place. It is \code{TRUE} only when a missing value is impossible:
\tabular{ll}{
\code{simd_vec(x, scan_na = TRUE)} \tab \code{!simd_any_na(x)} \cr
\code{abs}, \code{sign}, \code{floor}, \code{ceiling}, \code{trunc}, \code{round}, unary \code{-} \tab kept when \code{TRUE} \cr
\code{simd_pmin}, \code{simd_pmax} (and \code{_num}), \code{simd_clamp}, \code{cummax}, \code{cummin} \tab kept when every operand is \code{TRUE} \cr
\code{simd_conj}, \code{simd_re}, \code{simd_im}, \code{simd_popcount}, \code{simd_lzcnt}, \code{simd_tzcnt} \tab kept \cr
conversions that cannot overflow (to double, logical or raw; widening; saturating) \tab kept \cr
other arithmetic (\code{+}, \code{*}, \code{/}, \code{^}, the \code{_wrap} forms), elementary functions, other bitwise operations, \code{cumsum}, \code{cumprod}, base fallbacks \tab unknown \cr
\code{x[i]} \tab for a \code{TRUE} \code{x}, \code{TRUE}, or \code{FALSE} when \code{i} selects a missing element (\code{NA}, past the end, a name); otherwise unknown \cr
\code{x[i] <- value} \tab unknown \cr
\code{length(x) <- n} \tab \code{TRUE} kept when shortening; otherwise unknown \cr
\code{c(...)} \tab \code{FALSE} if any part is \code{FALSE}; \code{TRUE} if every \code{simd_vec} part is \code{TRUE} and the plain parts have no missing value; otherwise unknown \cr
}
Double arithmetic can make \code{NaN} from non-missing values
(\code{Inf - Inf}, \code{0 / 0}), and integer arithmetic \code{NA} on
overflow; a wrapped integer result of \code{-2^31} also reads as \code{NA},
so the \code{_wrap} forms do not keep the flag either. An operand of
length 1, plain or not, is checked directly: \code{x + 1} can still skip
the checks, and \code{abs(simd_vec(2))} is flagged \code{TRUE}.

\code{anyNA(x)} answers from the flag when it is known, without reading
the data.

\strong{Where the flag holds.} Base R copies attributes onto new data in
many places (\code{pmin()}, \code{storage.mode<-}, \code{ifelse()},
\code{unclass()}, \code{attributes<-} \dots), so the flag counts only
for the object it was set on. Its token points at that object, and the
flag reads as unknown (\code{simd_na_free()} gives \code{NULL}) on any
object the attributes were copied to, and also on the original once a
copy has been made: R cannot tell which of the two still holds the data
the flag describes. Any change of class voids the flag
too, even \code{class(x) <- NULL} followed by \code{class(x) <-
"simd_vec"}: on an unshared \code{x} the first assignment makes no copy,
so \code{x[i] <- NA} in between would change the data in place. A flag
reloaded by \code{readRDS()} or \code{unserialize()} is unknown, and R
code cannot set one (\code{attr<-} and \code{structure()} cannot make a
valid token). Two in-place routes stay open: putting back the class
vector that was taken off (\code{cl <- class(x); class(x) <- NULL;
x[1] <- NA; class(x) <- cl}), and C code that writes into \code{x} by
reference; take the data out with \code{simd_unwrap()} instead. Assigning into a copy
(\code{y <- x; y[1] <- 0}), passing \code{x} to functions, putting it in
a list, subsetting, sorting and printing leave the flag of \code{x}
alone; \code{attributes(x)}, \code{saveRDS(x)} and \code{serialize(x)}
void it. A lost flag only costs speed:
\code{simd_vec(x, scan_na = TRUE)} sets it again. Raw data is always
NA-free.

\strong{Operators.} \code{+}, \code{-}, \code{*}, \code{/}, \code{^},
\code{\%\%} and \code{\%/\%} call \code{\link{simd_add}},
\code{\link{simd_sub}}, \code{\link{simd_mul}}, \code{\link{simd_div}},
\code{\link{simd_pow}}, \code{\link{simd_mod}} and
\code{\link{simd_idiv}}; unary \code{-} calls \code{\link{simd_neg}}.
Type promotion follows base R (logical and integer give integer,
\code{/} and \code{^} give double), and operands must have equal lengths
or length 1. \code{==}, \code{!=}, \code{<}, \code{>}, \code{<=} and
\code{>=} call \code{\link{simd_eq}} and friends and return plain logical
vectors. Comparing with a character vector (\code{simd_vec(1:3) == "2"})
is an error, because \pkg{rsimd} rejects character input everywhere, where
base R would compare the numbers as strings. \code{&}, \code{|} and \code{!} call \code{\link{simd_and}},
\code{\link{simd_or}} and \code{\link{simd_not}}, or for raw data
\code{\link{simd_bit_and}}, \code{\link{simd_bit_or}} and
\code{\link{simd_bit_not}}. A plain operand may be on either side. An
operand with another class, such as a date, is an error, as in the
\code{simd_} functions: the \code{simd_vec} method is chosen over that
class's operator methods (\code{\link{chooseOpsMethod}}) and rejects it.
\pkg{bit64} always claims operations whose left operand is a plain
\code{integer64}, which then give a plain \code{integer64}.

\strong{Group generics.} The \code{Math} functions with a kernel
(\code{abs}, \code{sign}, \code{sqrt}, \code{floor}, \code{ceiling},
\code{trunc}, \code{round} with a single \code{digits}, \code{exp},
\code{expm1}, \code{log} with a single \code{base} (or, for complex data,
a \code{base} following the broadcast rule), \code{log2}, \code{log10},
\code{log1p}, the trigonometric, \code{pi}-scaled and hyperbolic
functions and their inverses, \code{cumsum}, \code{cumprod},
\code{cummax}, \code{cummin}) call the \code{simd_} functions; a longer
\code{digits}, or a longer \code{base} with non-complex data, is an error,
unlike in base R. The
\code{Summary} functions (\code{sum}, \code{prod}, \code{min},
\code{max}, \code{range}, \code{any}, \code{all}) and \code{mean} call
theirs with \code{na.rm} passed on; several arguments are combined with
\code{c()} first. \code{stats::var()} and \code{stats::sd()} are not
generic: use \code{\link{simd_var}} and \code{\link{simd_sd}}, or
\code{var(x)}, which works on the data as on any numeric vector.

\strong{Base fallbacks.} What has no kernel runs base R on the data and,
for a vector result, wraps it again with the operands' pin and an unknown
NA-free flag: \code{signif}, \code{gamma}, \code{lgamma},
\code{digamma}, \code{trigamma}; for complex data \code{\%\%},
\code{\%/\%}, the comparisons other than \code{==} and \code{!=} (plain
logical results), the \code{Math} functions without a complex kernel, and
\code{min}, \code{max}, \code{range}, \code{any} and \code{all}; for
\code{integer64} data \code{^}, \code{prod}, \code{mean} and the
\code{Math} functions other than \code{abs}, \code{sign}, \code{cumsum},
\code{cummax} and \code{cummin} that \pkg{bit64} has a method for
(\code{sqrt}, \code{log}, \code{floor}, \code{round}, \code{cumprod}
\dots), which need \pkg{bit64};
\code{mean(x, trim)} with \code{trim} other than 0; \code{range(x,
finite = TRUE)}; and \code{Re}, \code{Im}, \code{Mod}, \code{Arg} and
\code{Conj} of non-complex data. The other \code{Math} functions
(\code{exp}, \code{log1p}, the trigonometric functions, \code{gamma}
\dots) convert \code{integer64} data to double first, as \pkg{bit64}'s
\code{sqrt} does, since base R alone would read the 64-bit integers as
doubles. For complex data \code{Re}, \code{Im},
\code{Mod}, \code{Arg} and \code{Conj} call \code{\link{simd_re}},
\code{\link{simd_im}}, \code{\link{simd_abs}}, \code{\link{simd_arg}}
and \code{\link{simd_conj}}.

\strong{Other methods.} \code{x[i]} keeps the class and the pin
(\code{drop} is ignored); \code{x[i] <- value} and \code{x[[i]] <- value}
too (an \code{integer64} object converts \code{value} to
\code{integer64}). \code{value} must be data \pkg{rsimd} takes, as for an
operand: a character, list, factor, \code{Date} or other classed
\code{value} is an error rather than being stored as its codes or
numbers. A \code{simd_vec} \code{value} pinned to another implementation
than \code{x} is an error, as in \code{c()}; one pinned while \code{x} is
not leaves \code{x} unpinned. \code{x[[i]]} returns a plain element.
\code{is.na(x)} calls \code{\link{simd_is_na}} and returns a plain
logical vector.
\code{match()}, \code{\%in\%} and so \code{setdiff()}, \code{union()}
and \code{intersect()} compare the data, as for plain data.
A \code{simd_vec} has no names, dimensions or dimnames. Setting them with a
value other than \code{NULL} (\code{names<-}, \code{dim<-},
\code{dimnames<-}) unwraps: the result is the plain data with the value
set, not a \code{simd_vec}. So \code{quantile(x)}, \code{summary(x)},
\code{setNames(x, nm)} and \code{lm(x ~ z)} give the same results as for
plain data, and \code{table()}, \code{matrix()} and \code{outer()} work
as before. \code{median(x)} runs on the data and also returns what it does
for plain data. Setting \code{NULL} (\code{unname(x)})
leaves a \code{simd_vec}.
\code{c()} combines with base R's type promotion, except that
\code{integer64} absorbs the other parts, as in \pkg{bit64}, and the
parts' pins must agree. \code{c()} uses the \code{simd_vec} method only
when its first argument is a \code{simd_vec}; \code{c(5, x)} returns a
plain vector, and \code{c(1L, x)} for \code{integer64} data gives
\pkg{bit64}'s garbage, as for any \code{integer64} (put the
\code{simd_vec} first, or use \code{c(simd_vec(5), x)}).
\code{length(x) <- n} pads with \code{NA}.
\code{rev()}, \code{sort()}, \code{head()} and \code{tail()} work through
\code{[} and so keep the class; \code{rep()} and \code{unique()} return
plain vectors.

\strong{Unwrapping.} \code{simd_unwrap(x)} returns the data without any
\pkg{rsimd} attributes; \code{integer64} data keeps its class. Given plain
data, it returns it unchanged. \code{bit64::as.integer64(x)} does the same
for \code{integer64} data and converts other data.
\code{as.vector()}, \code{as.double()}, \code{as.numeric()},
\code{as.integer()} and \code{as.logical()} unwrap to plain vectors of
their type: \code{integer64} data becomes double for \code{as.vector()},
as in \pkg{bit64}. \code{as.data.frame()} (so \code{data.frame(v = x)})
unwraps too. Do not unwrap with \code{unclass(x)}: it removes the whole
class, \code{"integer64"} included, and leaves the \pkg{rsimd} attributes
behind (voiding the flag, as above). \code{print()} shows a header
line before the data, \code{format()} formats the data and \code{str()}
gives a one-line summary.

\code{all.equal(x, y)} compares the data, the pins and the flags. The
token makes two separately built flagged objects differ for
\code{identical()}, and \code{dput()} of a flagged object does not parse
back.
}
\value{
\code{simd_vec()}, \code{as_simd_vec()} and \code{simd_impl<-}: a
\code{simd_vec}. \code{is_simd_vec()}: \code{TRUE} or \code{FALSE}.
\code{simd_impl()}: the pinned implementation or \code{NULL}.
\code{simd_na_free()}: \code{TRUE}, \code{FALSE} or \code{NULL}.
\code{simd_unwrap()}: a plain vector, of class \code{"integer64"} for
\code{integer64} data.
\code{all.equal()}: \code{TRUE} or a character vector describing the
differences.
}
\seealso{
\code{\link{simd_use}}, \code{\link{simd_with_impl}},
\code{\link{simd_any_na}}, \link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
x <- simd_vec(c(1.5, -2, 4), scan_na = TRUE)
x
x * 2 + 1
sum(abs(x))
exp(x)
x > 0
simd_na_free(abs(x))
simd_na_free(x + 1)

# Pin to an implementation; results carry the pin.
y <- simd_vec(1:4, impl = simd_available()[1])
simd_impl(y + 1L)
simd_impl(y) <- NULL

str(x)
simd_unwrap(x)

# Setting names unwraps, so base functions that name their result work.
quantile(x)
}
