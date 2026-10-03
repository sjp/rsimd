\name{simd_cpu_features}
\alias{simd_cpu_features}
\title{CPU Features Detected at Run Time}
\usage{
simd_cpu_features()
}
\description{
Reports the SIMD instruction set extensions that the CPU and operating
system support, as detected when \pkg{rsimd} was loaded. This describes the
hardware, not the package: an extension can be reported here even if the
package was built without code for it. \code{simd_available()} lists the
implementations that are both compiled in and supported.
}
\details{
Detection uses \code{cpuid} on x86 and x86-64, \code{getauxval()} on Linux
arm, \code{elf_aux_info()} on FreeBSD arm, \code{sysctlbyname()} on macOS
arm64 and \code{IsProcessorFeaturePresent()} on Windows arm64.

\code{os_avx} and \code{os_avx512} are \code{TRUE} when the operating system
has enabled the AVX and AVX-512 register state. Without that, the CPU
instructions cannot be used.

Features are reported consistently: a feature is \code{TRUE} only if the
features it builds on are too. For example, \code{avx2} needs \code{avx},
\code{avx} needs \code{os_avx}, all \code{avx512*} features need
\code{avx512f}, and \code{sve2} needs \code{sve}.

The environment variable \env{RSIMD_CPU_FEATURES_MASK} is intended for testing
only. Set it before the package is loaded to a comma-separated list of
feature names, for example \code{"avx512f,avx2"}. Those features are reported
as unsupported, together with every feature that builds on them, so the
fallback implementations are used. It never turns a feature on. Unknown
names give a warning when the package is loaded.
}
\value{
A list with elements:
\item{arch}{\code{"x86_64"}, \code{"i686"}, \code{"aarch64"}, \code{"armv7"}
  or \code{"other"}.}
\item{os}{\code{"linux"}, \code{"darwin"}, \code{"windows"},
  \code{"freebsd"} or \code{"other"}.}
\item{method}{How the features were detected: \code{"cpuid"},
  \code{"getauxval"}, \code{"sysctl"}, \code{"win32"},
  \code{"elf_aux_info"} or \code{"none"}.}
\item{features}{A named logical vector with, in this order, \code{sse2},
  \code{sse3}, \code{ssse3}, \code{sse4_1}, \code{sse4_2}, \code{avx},
  \code{avx2}, \code{fma}, \code{avx512f}, \code{avx512bw}, \code{avx512dq},
  \code{avx512vl}, \code{os_avx}, \code{os_avx512} (x86), and \code{neon},
  \code{fp16}, \code{dotprod}, \code{i8mm}, \code{bf16}, \code{sve},
  \code{sve2} (arm).}
\item{sve_vector_length_bits}{The SVE vector length in bits, or 0 when SVE
  is absent or the length is unknown.}
\item{masked}{The names of features that are supported but were switched off
  by \env{RSIMD_CPU_FEATURES_MASK}. Normally \code{character(0)}.}
}
\examples{
f <- simd_cpu_features()
f$arch
names(f$features)[f$features]
}
