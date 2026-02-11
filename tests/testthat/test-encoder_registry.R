# Tests for encoder registry, utilities, and decoder validation

# --- Encoder registry ---------------------------------------------------------

test_that("register / list / get cycle works", {
  # Clean slate: remove any test entries added by previous tests
  if (exists("__test_enc__", envir = fmrilatent:::.encoder_registry_env, inherits = FALSE)) {
    rm("__test_enc__", envir = fmrilatent:::.encoder_registry_env)
  }

  register_encoder("__test_enc__", identity, "A test encoder", "testpkg")
  enc <- get_encoder("__test_enc__")
  expect_identical(enc$spec_fn, identity)
  expect_equal(enc$description, "A test encoder")
  expect_equal(enc$package, "testpkg")

  tbl <- list_encoders()
  expect_true("__test_enc__" %in% tbl$family)

  # Clean up

  rm("__test_enc__", envir = fmrilatent:::.encoder_registry_env)
})

test_that("built-in encoders are registered after loading", {
  tbl <- list_encoders()
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("family", "description", "package") %in% names(tbl)))

  expected <- c("time_slepian", "time_dct", "time_bspline",
                "space_slepian", "space_heat", "space_hrbf",
                "space_wavelet_active", "st", "hierarchical")
  for (fam in expected) {
    expect_true(fam %in% tbl$family, info = paste("Missing:", fam))
  }
  expect_true(all(tbl$package[tbl$family %in% expected] == "fmrilatent"))
  expect_true(all(nzchar(tbl$description[tbl$family %in% expected])))
})

test_that("duplicate registration warns and overwrites", {
  register_encoder("__dup_test__", identity, "first", "pkg1")
  expect_warning(
    register_encoder("__dup_test__", sum, "second", "pkg2"),
    "already registered"
  )
  enc <- get_encoder("__dup_test__")
  expect_identical(enc$spec_fn, sum)
  expect_equal(enc$description, "second")

  rm("__dup_test__", envir = fmrilatent:::.encoder_registry_env)
})

test_that("get_encoder errors on nonexistent family", {
  expect_error(get_encoder("__nonexistent__"), "not registered")
})

test_that("list_encoders returns empty data.frame when registry is empty", {
  # Save current state
  saved <- as.list(fmrilatent:::.encoder_registry_env)
  rm(list = ls(fmrilatent:::.encoder_registry_env),
     envir = fmrilatent:::.encoder_registry_env)

  tbl <- list_encoders()
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 0L)
  expect_true(all(c("family", "description", "package") %in% names(tbl)))

  # Restore
  for (nm in names(saved)) {
    assign(nm, saved[[nm]], envir = fmrilatent:::.encoder_registry_env)
  }
})

# --- Test data helper ---------------------------------------------------------

test_that("fmrilatent_test_data returns correct structure", {
  td <- fmrilatent_test_data()
  expect_type(td, "list")
  expect_true(all(c("X", "mask", "dims", "n_time") %in% names(td)))
  expect_true(is.matrix(td$X))
  expect_equal(dim(td$X), c(8L, 18L))
  expect_true(is.array(td$mask))
  expect_equal(dim(td$mask), c(3L, 3L, 2L))
  expect_true(all(td$mask))
})

test_that("fmrilatent_test_data respects custom dims", {
  td <- fmrilatent_test_data(dims = c(4, 4, 3), n_time = 12)
  expect_equal(dim(td$X), c(12L, 48L))
  expect_equal(dim(td$mask), c(4L, 4L, 3L))
  expect_equal(td$n_time, 12L)
})

# --- mask_to_array export -----------------------------------------------------

test_that("mask_to_array is exported and works", {
  arr <- array(TRUE, dim = c(3, 3, 2))
  result <- mask_to_array(arr)
  expect_true(is.array(result))
  expect_equal(dim(result), c(3L, 3L, 2L))
})

test_that("mask_to_array errors on bad input", {
  expect_error(mask_to_array(new.env(), location = "test"), "mask must be")
})

# --- roi_subset_columns -------------------------------------------------------

test_that("roi_subset_columns returns unchanged when roi_mask is NULL", {
  mat <- matrix(1:6, nrow = 2, ncol = 3)
  result <- roi_subset_columns(mat, array(TRUE, c(3, 1, 1)), roi_mask = NULL)
  expect_identical(result, mat)
})

test_that("roi_subset_columns subsets correctly", {
  # mask has 3 TRUE voxels at positions 1, 2, 4
  mask <- array(c(TRUE, TRUE, FALSE, TRUE), dim = c(2, 2, 1))
  rec  <- matrix(1:9, nrow = 3, ncol = 3)
  # roi selects positions 1 and 4 → columns 1 and 3

  roi  <- array(c(TRUE, FALSE, FALSE, TRUE), dim = c(2, 2, 1))

  result <- roi_subset_columns(rec, mask, roi)
  expect_equal(ncol(result), 2L)
  expect_equal(result[, 1], rec[, 1])
  expect_equal(result[, 2], rec[, 3])
})

# --- Decoder contract validation ----------------------------------------------

test_that("implicit_latent validates decoder is a function", {
  expect_error(
    implicit_latent(list(), "not_a_function", list(family = "f"),
                    array(TRUE, c(2, 2, 1))),
    "decoder must be a function"
  )
})

test_that("implicit_latent rejects decoder missing required formals", {
  bad_decoder <- function(x) x
  expect_error(
    implicit_latent(list(), bad_decoder, list(family = "f"),
                    array(TRUE, c(2, 2, 1))),
    "Missing.*time_idx.*roi_mask"
  )
})

test_that("implicit_latent accepts decoder with ... formals", {
  ok_decoder <- function(...) NULL
  il <- implicit_latent(list(), ok_decoder, list(family = "f"),
                        array(TRUE, c(2, 2, 1)))
  expect_s3_class(il, "ImplicitLatent")
})

test_that("implicit_latent accepts decoder with named formals", {
  ok_decoder <- function(time_idx = NULL, roi_mask = NULL) NULL
  il <- implicit_latent(list(), ok_decoder, list(family = "f"),
                        array(TRUE, c(2, 2, 1)))
  expect_s3_class(il, "ImplicitLatent")
})

test_that("implicit_latent accepts decoder with extra formals", {
  ok_decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) NULL
  il <- implicit_latent(list(), ok_decoder, list(family = "f"),
                        array(TRUE, c(2, 2, 1)))
  expect_s3_class(il, "ImplicitLatent")
})
