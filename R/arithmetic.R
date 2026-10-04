# Elementwise arithmetic. Each function checks the types the C side cannot
# name the function for (complex, integer64, double for the _wrap ops) and
# calls one of four entry points with the op's name; the C side applies the
# broadcast rule, picks the kernel and warns.

# Functions that will take complex or integer64 operands: their messages
# say "yet".
.ew_later_complex <- c("simd_add", "simd_sub", "simd_neg")
.ew_later_i64 <- c(
  "simd_add", "simd_sub", "simd_mul", "simd_add_wrap", "simd_sub_wrap", "simd_mul_wrap",
  "simd_neg", "simd_abs", "simd_neg_wrap", "simd_abs_wrap", "simd_idiv", "simd_mod",
  "simd_pmin", "simd_pmax", "simd_pmin_num", "simd_pmax_num", "simd_clamp", "simd_sign"
)

# Checks every operand in `args` (a named list) for function `fun`; the
# _wrap ops (`wrap = TRUE`) also reject doubles.
.ew_check <- function(fun, args, wrap = FALSE) {
  unsupported <- c("integer64", "complex", if (wrap) "double")
  later <- c(
    if (fun %in% .ew_later_complex) "complex",
    if (fun %in% .ew_later_i64) "integer64"
  )
  for (arg in names(args)) .check_supported(args[[arg]], fun, unsupported, later, arg)
  invisible()
}

.ew1 <- function(x, op, fun, wrap = FALSE) {
  .sync_impl()
  .ew_check(fun, list(x = x), wrap)
  .Call(C_simd_ew1, x, op)
}

.ew2 <- function(x, y, op, fun, na_check, wrap = FALSE, names = c("x", "y")) {
  .sync_impl()
  args <- list(x, y)
  names(args) <- names
  .ew_check(fun, args, wrap)
  .Call(C_simd_ew2, x, y, op, na_check)
}

.ew3 <- function(x, y, z, op, fun, na_check, names = c("x", "y", "z")) {
  .sync_impl()
  args <- list(x, y, z)
  names(args) <- names
  .ew_check(fun, args)
  .Call(C_simd_ew3, x, y, z, op, na_check)
}

# na.rm of pmin/pmax: TRUE or FALSE.
.arg_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
  }
  x
}

simd_add <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "add", "simd_add", na_check)
}

simd_sub <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "sub", "simd_sub", na_check)
}

simd_mul <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mul", "simd_mul", na_check)
}

simd_div <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "div", "simd_div", na_check)
}

simd_idiv <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "idiv", "simd_idiv", na_check)
}

simd_mod <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mod", "simd_mod", na_check)
}

simd_add_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "add_wrap", "simd_add_wrap", na_check, wrap = TRUE)
}

simd_sub_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "sub_wrap", "simd_sub_wrap", na_check, wrap = TRUE)
}

simd_mul_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mul_wrap", "simd_mul_wrap", na_check, wrap = TRUE)
}

simd_neg <- function(x) .ew1(x, "neg", "simd_neg")

simd_abs <- function(x) .ew1(x, "abs", "simd_abs")

simd_neg_wrap <- function(x) .ew1(x, "neg", "simd_neg_wrap", wrap = TRUE)

simd_abs_wrap <- function(x) .ew1(x, "abs", "simd_abs_wrap", wrap = TRUE)

simd_sign <- function(x) .ew1(x, "sign", "simd_sign")

simd_copysign <- function(x, sign, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, sign, "copysign", "simd_copysign", na_check, names = c("x", "sign"))
}

simd_recip <- function(x) .ew1(x, "recip", "simd_recip")

simd_sqrt <- function(x) .ew1(x, "sqrt", "simd_sqrt")

simd_fma <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "fma", "simd_fma", na_check)
}

simd_mul_add <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "mul_add", "simd_mul_add", na_check)
}

simd_add_mul <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "add_mul", "simd_add_mul", na_check)
}

simd_lerp <- function(x, y, t, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, t, "lerp", "simd_lerp", na_check, names = c("x", "y", "t"))
}

simd_pmin <- function(x, y, na.rm = FALSE) {
  op <- if (.arg_flag(na.rm, "na.rm")) "pmin_num" else "pmin"
  .ew2(x, y, op, "simd_pmin", NULL)
}

simd_pmax <- function(x, y, na.rm = FALSE) {
  op <- if (.arg_flag(na.rm, "na.rm")) "pmax_num" else "pmax"
  .ew2(x, y, op, "simd_pmax", NULL)
}

simd_pmin_num <- function(x, y) .ew2(x, y, "pmin_num", "simd_pmin_num", NULL)

simd_pmax_num <- function(x, y) .ew2(x, y, "pmax_num", "simd_pmax_num", NULL)

simd_clamp <- function(x, lo, hi) {
  .ew3(x, lo, hi, "clamp", "simd_clamp", NULL, names = c("x", "lo", "hi"))
}

simd_floor <- function(x) .ew1(x, "floor", "simd_floor")

simd_ceiling <- function(x) .ew1(x, "ceiling", "simd_ceiling")

simd_trunc <- function(x) .ew1(x, "trunc", "simd_trunc")

simd_round <- function(x, digits = 0) {
  if (!(is.numeric(digits) || is.logical(digits)) || length(digits) != 1L) {
    stop("'digits' must be a single number", call. = FALSE)
  }
  if (!is.na(digits) && digits == 0) return(.ew1(x, "round", "simd_round"))
  .sync_impl()
  .ew_check("simd_round", list(x = x))
  .Call(C_simd_round_digits, x, digits)
}
