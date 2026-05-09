library(testthat)
library(fmrilatent)

# Many package tests intentionally exercise internal helpers. `devtools::test()`
# exposes those helpers via load_all(); mirror that lookup behavior for
# installed-package checks without adding a pkgload dependency.
fmrilatent_ns <- asNamespace("fmrilatent")
for (nm in ls(fmrilatent_ns, all.names = TRUE)) {
  if (!exists(nm, envir = .GlobalEnv, inherits = TRUE)) {
    assign(nm, get(nm, envir = fmrilatent_ns), envir = .GlobalEnv)
  }
}

test_check("fmrilatent")
