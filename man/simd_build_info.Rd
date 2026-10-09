\name{simd_build_info}
\alias{simd_build_info}
\title{How rsimd Was Built}
\usage{
simd_build_info()
}
\description{
Reports how the installed \pkg{rsimd} was built: the versions of the bundled
SIMDe and SLEEF libraries, the compiler, the implementations compiled in,
and which of them this machine can run. Please include its output, together
with that of \code{\link{simd_cpu_features}()}, when reporting a bug.
}
\value{
A list with elements:
\item{simde}{The version of the bundled SIMDe headers, such as
  \code{"0.8.4"}.}
\item{sleef}{The version of the bundled SLEEF headers, such as
  \code{"3.9.0"}, or \code{NA} when no implementation uses them (the
  compiler could not build them, or \env{RSIMD_DISABLE_SLEEF} was set at
  install time; see \link{rsimd_options}).}
\item{compiler}{The C compiler that built the package, such as
  \code{"GCC 14.2.0"} or \code{"Clang 18.1.8"}, or \code{NA} if unknown.}
\item{compiled}{The implementations compiled into the package, in
  \code{\link{simd_tiers}()} order. A tier is compiled when the compiler
  could build it, whatever the machine that built it.}
\item{sleef_tiers}{The compiled implementations that use SLEEF for
  elementary functions. The others take them from a lower tier, and
  \code{none} from the C library.}
\item{available}{\code{\link{simd_available}()}: the compiled
  implementations this CPU and operating system can run.}
\item{current}{\code{\link{simd_current}()}: the implementation in use.}
}
\seealso{
\code{\link{simd_cpu_features}} for what the CPU supports,
\code{\link{simd_use}} to choose an implementation.
}
\examples{
str(simd_build_info())
}
