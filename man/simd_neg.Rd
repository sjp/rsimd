\name{simd_neg}
\alias{simd_neg}
\alias{simd_abs}
\alias{simd_neg_wrap}
\alias{simd_abs_wrap}
\alias{simd_sign}
\alias{simd_copysign}
\alias{simd_recip}
\alias{simd_sqrt}
\title{Negation, Absolute Value, Sign, Reciprocal and Square Root}
\usage{
simd_neg(x)

simd_abs(x)

simd_neg_wrap(x)

simd_abs_wrap(x)

simd_sign(x)

simd_copysign(x, sign, na_check = getOption("rsimd.na_check", TRUE))

simd_recip(x)

simd_sqrt(x)
}
\arguments{
\item{x}{A double, integer or logical vector; integer or logical only for
the \code{_wrap} functions. \code{simd_neg()} also takes a complex
vector, and \code{simd_neg()}, \code{simd_abs()}, \code{simd_sign()} and the
\code{_wrap} functions a \pkg{bit64} \code{integer64} vector.}
\item{sign}{A double, integer or logical vector whose signs are copied to
\code{x}.}
\item{na_check}{\code{TRUE} or \code{FALSE}. \code{FALSE} promises that
\code{x} and \code{sign} have no missing values, so they are not looked for;
see Details.}
}
\description{
\code{-x}, \code{abs(x)}, \code{sign(x)}, \code{1 / x} and \code{sqrt(x)}
computed with SIMD instructions, and \code{simd_copysign(x, sign)}, which is
\code{abs(x)} with the sign of \code{sign}. They correspond to
\code{Negate}, \code{Abs}, \code{Sign}, \code{CopySign}, \code{Reciprocal}
and \code{Sqrt} of .NET's \code{TensorPrimitives}.
}
\details{
\code{simd_neg()} and \code{simd_abs()} keep integer and logical input
integer, like \code{-} and \code{abs()}. They cannot overflow: the only
integer whose negation is out of range is \code{-2147483648}, which is
\code{NA}. The \code{_wrap} spellings are the same functions, provided for
symmetry with \code{\link{simd_add_wrap}()}, and reject doubles.
\code{simd_neg()} of a complex vector negates both parts, as \code{-z}
does (see \code{\link{simd_conj}()} for the conjugate).

\code{simd_sign()}, \code{simd_recip()} and \code{simd_sqrt()} return
doubles for every input, as base R does for \code{sign()} and
\code{sqrt()}. \code{simd_sign()} is \code{-1}, \code{0} or \code{1}
(\code{0}, not \code{-0}, for \code{-0}); \code{simd_recip(0)} is
\code{Inf} and \code{simd_recip(-0)} is \code{-Inf}. \code{simd_sqrt()} is
correctly rounded and warns \code{"NaNs produced"} once when a negative
number gives \code{NaN}.

\code{simd_copysign()} sets the sign bit of \code{x}, including for
\code{0}, infinities and \code{NaN}. A missing \code{sign} gives a missing
result: \code{NA} for \code{NA}, \code{NaN} for \code{NaN}; an \code{NA} in
either operand gives \code{NA}. Its operands follow the broadcast rule:

\rsimdBroadcast

\rsimdSemantics

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM the double operations use the scalar \code{"none"} code.
Integer64 input gives integer64 results (\code{NA} stays \code{NA}; see
\link{rsimd-integer64}).
}
\value{
A vector of the length of \code{x} (for \code{simd_copysign()}, of the
longer operand): integer for \code{simd_neg()} and \code{simd_abs()} of
integer or logical input, complex for \code{simd_neg()} of complex input,
otherwise double.
}
\seealso{
\code{\link{abs}}, \code{\link{sign}}, \code{\link{sqrt}},
\code{\link{simd_add}}, \link{rsimd_options}, \link{rsimd-integer64}.
}
\examples{
simd_neg(c(1L, -2L, NA))
simd_abs(c(-1.5, 2, -Inf))
simd_sign(c(-2, 0, 3))
simd_copysign(c(1, -2, 3), c(-1, -1, 1))
simd_recip(c(2, 0, -0))
simd_sqrt(c(4, 9, 2))
}
