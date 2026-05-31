skip_if_no_neurosurf_surface_examples <- function() {
  check_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  skip_if(
    any(commandArgs() == "-f") || grepl("[.]Rcheck($|/)", check_path),
    "neurosurf surface examples abort under this local R/BATCH setup"
  )
  skip_if_not_installed("neurosurf")
}
