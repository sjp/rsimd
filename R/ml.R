# Machine-learning helpers. simd_sigmoid() is an elementwise function
# with the elementary functions' rules; simd_softmax() and
# simd_log_softmax() normalise over the whole vector, summing in the
# precision mode of option rsimd.precision.

simd_sigmoid <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sigmoid", list(x = x))
  .Call(C_simd_math1, x, "sigmoid", NULL)
}

simd_softmax <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_softmax", list(x = x))
  .Call(C_simd_softmax, x, NULL)
}

simd_log_softmax <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_log_softmax", list(x = x))
  .Call(C_simd_log_softmax, x, NULL)
}
