## Submission notes

### Pragmas suppressing diagnostics

`R CMD check` notes that some files under `src/vendor/simde` contain pragmas
suppressing diagnostics. These are unmodified headers of the vendored, MIT-licensed
SIMDe library (`hedley.h`, `simde-common.h`, `x86/mmx.h`, `x86/sse.h`, `x86/sse2.h`).
The pragmas are scoped by push/pop pragmas to SIMDe's own code and silence clang-only
or pedantic diagnostics (`-Wvariadic-macros`, `-Wvector-conversion`,
`-Wc11-extensions` and similar); none concerns uninitialised values, formats or array
bounds. rsimd's own sources contain no such pragmas. We keep the vendored headers
identical to upstream so that they can be updated mechanically.
