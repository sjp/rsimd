\name{simd_conj}
\alias{simd_conj}
\alias{simd_re}
\alias{simd_im}
\alias{simd_arg}
\title{Complex Conjugate, Real and Imaginary Parts, Argument}
\usage{
simd_conj(z)

simd_re(z)

simd_im(z)

simd_arg(z)
}
\arguments{
\item{z}{A complex, double, integer or logical vector.}
}
\description{
\code{Conj(z)}, \code{Re(z)}, \code{Im(z)} and \code{Arg(z)} computed
with SIMD instructions.
}
\details{
R stores a complex vector as pairs of doubles, the real part followed by
the imaginary part. \code{simd_conj()} flips the sign of each imaginary
part, so an imaginary \code{0} becomes \code{-0}, and the real part is
unchanged, \code{NA} included. \code{simd_re()} and \code{simd_im()} take
every other double; missing parts are returned as they are, so
\code{simd_im(complex(real = NA, imaginary = 1))} is \code{1}.

For double, integer and logical \code{z}, \code{simd_conj()} and
\code{simd_re()} return \code{z} as a double vector and \code{simd_im()}
returns zeros, as base R does; \code{simd_arg()} treats them as complex
numbers with a zero imaginary part, so it is \code{0} or \code{pi}.

\code{simd_arg(z)} is \code{atan2(Im(z), Re(z))}, as base R's
\code{Arg()} computes it with the C library's \code{carg()}. The
\code{"none"} implementation calls the C library; the SIMD
implementations use SLEEF's \code{atan2} (within 1 ULP, or 3.5 ULP in
\link[=simd_math_accuracy]{fast mode}) when both parts are finite, and
the C library otherwise, so infinite and missing parts give exactly base
R's result. \code{\link{simd_abs}()} of a complex vector gives
\code{Mod(z)} in the same way, with SLEEF's \code{hypot} (within 0.5
ULP, or 3.5 ULP in fast mode).

\strong{Complex numbers in rsimd.} Besides these functions,
\code{\link{simd_add}()}, \code{\link{simd_sub}()},
\code{\link{simd_mul}()}, \code{\link{simd_div}()},
\code{\link{simd_neg}()}, \code{\link{simd_abs}()},
\code{\link{simd_sum}()}, \code{\link{simd_prod}()},
\code{\link{simd_mean}()}, \code{\link{simd_cumsum}()},
\code{\link{simd_cumprod}()}, \code{\link{simd_eq}()},
\code{\link{simd_ne}()}, \code{\link{simd_hamming}()}, the
predicates \code{\link{simd_is_na}()}, \code{simd_is_nan()},
\code{simd_is_finite()} and \code{simd_is_infinite()} (with their
\code{_any} and \code{_all} forms) and \code{\link{simd_any_na}()},
\code{simd_count_na()} and \code{simd_which_na()} take complex input.
Multiplication, division, \code{cumsum()} and \code{cumprod()} give
exactly base R's results. The elementary functions and powers reject
complex input; they are planned for a later version. The other
comparisons give base R's error, and \code{min()}, \code{max()} and the
like are not defined for complex numbers.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM they use the scalar \code{"none"} code, as do SIMD
implementations built without SLEEF for \code{simd_arg()}.
}
\value{
\code{simd_conj()}: a complex vector for complex \code{z}, otherwise a
double vector. \code{simd_re()}, \code{simd_im()} and \code{simd_arg()}:
a double vector. All have the length of \code{z}.
}
\seealso{
\code{\link{Conj}}, \code{\link{Re}}, \code{\link{Im}},
\code{\link{Arg}}, \code{\link{simd_add}}, \code{\link{simd_abs}},
\code{\link{simd_sum}}.
}
\examples{
z <- complex(real = c(1, -2, NA), imaginary = c(0.5, 3, 1))
simd_conj(z)
simd_re(z)
simd_im(z)
simd_arg(z)
simd_re(1:3)
simd_arg(c(-1, 2))
}
