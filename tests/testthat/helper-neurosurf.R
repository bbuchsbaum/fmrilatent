.neurosurf_namespace_loadable <- local({
  loadable <- NULL
  function() {
    if (is.null(loadable)) {
      expression <- paste0(
        "quit(status = if (requireNamespace('neurosurf', quietly = TRUE)) ",
        "0L else 1L)"
      )
      command <- paste0(
        "(", shQuote(file.path(R.home("bin"), "Rscript")),
        " --vanilla -e ", shQuote(expression),
        " >/dev/null 2>&1) 2>/dev/null"
      )
      status <- suppressWarnings(system(command))
      loadable <<- identical(status, 0L)
    }
    loadable
  }
})

skip_if_no_neurosurf_surface_examples <- function() {
  skip_if_not(
    nzchar(system.file(package = "neurosurf")),
    "neurosurf is not installed"
  )
  skip_if_not(
    .neurosurf_namespace_loadable(),
    "neurosurf is installed but its namespace cannot be loaded"
  )
}
