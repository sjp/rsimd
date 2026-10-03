# Restyle the package sources with styler, if it is installed.
# Run from the package root: Rscript tools/style.R
if (requireNamespace("styler", quietly = TRUE)) {
  styler::style_pkg()
} else {
  message("styler is not installed; nothing to do.")
}
