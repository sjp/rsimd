# Many exported functions differ only in the names they pass on: the op
# the C side runs and the function name their checks report. Each such
# family is made by a maker in its own file (.pred(), .ew2(), .math1(), ...)
# from one bquote() template, and each function is defined with one line,
# as in `simd_is_zero <- .pred("zero", 0L)`. The names are put into the
# body, so the function made is an ordinary closure holding a literal
# .Call(): it makes its own .Call(), so warnings and errors from the C side
# name the user's call, it costs nothing extra per call, and printing it
# shows the code that runs. This file sorts before the files that use it.

# A function with formals `args` (an alist()) and body `body`, with the
# package namespace as its environment.
.wrapper <- function(args, body) as.function(c(args, body), topenv(environment()))
