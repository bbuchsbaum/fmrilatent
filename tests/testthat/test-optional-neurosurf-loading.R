test_that("package load does not probe the optional neurosurf namespace", {
  on_load_source <- paste(deparse(body(fmrilatent:::.onLoad)), collapse = "\n")

  expect_false(grepl("neurosurf", on_load_source, fixed = TRUE))
  expect_false(grepl("requireNamespace", on_load_source, fixed = TRUE))
})

test_that("optional bilateral surface conversion uses registered S3 dispatch", {
  expect_true(is.function(getS3method(
    "as.matrix", "BilatNeuroSurfaceVector", optional = TRUE
  )))
  expect_false(exists(
    ".register_bilat_neurosurfacevector_methods",
    envir = asNamespace("fmrilatent"),
    inherits = FALSE
  ))
})
