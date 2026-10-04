\name{simd_conj}
\alias{simd_conj}
\alias{simd_re}
\alias{simd_im}
\title{Complex Conjugate, Real and Imaginary Parts}
\usage{
simd_conj(z)

simd_re(z)

simd_im(z)
}
\arguments{
\item{z}{A complex, double, integer or logical vector.}
}
\description{
\code{Conj(z)}, \code{Re(z)} and \code{Im(z)} computed with SIMD
instructions.
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
returns zeros, as base R does.

\strong{Complex numbers in rsimd.} Besides these three functions,
\code{\link{simd_add}()}, \code{\link{simd_sub}()}, \code{\link{simd_neg}()},
\code{\link{simd_sum}()}, the predicates \code{\link{simd_is_na}()},
\code{simd_is_nan()}, \code{simd_is_finite()} and
\code{simd_is_infinite()} (with their \code{_any} and \code{_all} forms)
and \code{\link{simd_any_na}()}, \code{simd_count_na()} and
\code{simd_which_na()} take complex input. These are the operations that
treat the two parts separately. The other functions, including
multiplication, division, \code{Mod()}, \code{Arg()}, comparisons and the
elementary functions, reject complex input; they are planned for a later
version.

\rsimdAttributes

All implementations have their own kernels for these functions, except that
on 32-bit ARM they use the scalar \code{"none"} code.
}
\value{
\code{simd_conj()}: a complex vector for complex \code{z}, otherwise a
double vector. \code{simd_re()} and \code{simd_im()}: a double vector. All
have the length of \code{z}.
}
\seealso{
\code{\link{Conj}}, \code{\link{Re}}, \code{\link{Im}},
\code{\link{simd_add}}, \code{\link{simd_sum}}.
}
\examples{
z <- complex(real = c(1, -2, NA), imaginary = c(0.5, 3, 1))
simd_conj(z)
simd_re(z)
simd_im(z)
simd_re(1:3)
}
