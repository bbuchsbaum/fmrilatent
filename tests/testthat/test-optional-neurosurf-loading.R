test_that("package load does not probe the optional neurosurf namespace", {
  on_load_source <- paste(deparse(body(fmrilatent:::.onLoad)), collapse = "\n")

  expect_false(grepl("neurosurf", on_load_source, fixed = TRUE))
  expect_false(grepl("requireNamespace", on_load_source, fixed = TRUE))
})

test_that("bilateral wrapping owns lazy neurosurf method registration", {
  method_source <- paste(
    deparse(body(methods::selectMethod(
      "wrap_decoded", "BilatLatentNeuroSurfaceVector"
    ))),
    collapse = "\n"
  )

  expect_true(grepl(
    ".register_bilat_neurosurfacevector_methods", method_source, fixed = TRUE
  ))
})
