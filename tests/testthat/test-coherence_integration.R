library(testthat)
library(neuroim2)

make_test_data_ci_int <- function(n_time = 8, dims = c(4, 4, 2)) {
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims))
  n_vox <- sum(mask_arr)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  list(X = X, mask_arr = mask_arr, mask_vol = mask_vol, n_vox = n_vox, n_time = n_time)
}

test_that("common latent interface works on real explicit encoders", {
  td <- make_test_data_ci_int(n_time = 8, dims = c(3, 3, 2))
  roi <- array(FALSE, dim = dim(td$mask_arr))
  roi[1:2, , ] <- TRUE

  lv_dct <- encode(td$X, spec_time_dct(k = 3), mask = td$mask_vol, materialize = "matrix")
  lv_pca <- encode(td$X, spec_space_pca(k = 2, backend = "svd"), mask = td$mask_vol, materialize = "matrix")

  expect_true(is_explicit_latent(lv_dct))
  expect_true(is_explicit_latent(lv_pca))
  expect_equal(reconstruct_matrix(lv_dct), as.matrix(lv_dct))
  expect_equal(reconstruct_matrix(lv_pca), as.matrix(lv_pca))
  expect_equal(
    reconstruct_matrix(lv_pca, time_idx = 2:4, roi_mask = roi),
    as.matrix(lv_pca)[2:4, which(as.logical(roi)), drop = FALSE]
  )
  expect_equal(dim(reconstruct_array(lv_dct)), c(dim(td$mask_arr), td$n_time))
  expect_equal(dim(reconstruct_array(lv_pca, time_idx = 2:3)), c(dim(td$mask_arr), 2))
  expect_true(nzchar(latent_meta(lv_dct)$family))
  expect_true(nzchar(latent_meta(lv_pca)$family))
})

test_that("common latent interface works on real implicit encoders", {
  td <- make_test_data_ci_int(n_time = 8, dims = c(4, 4, 2))
  roi <- array(FALSE, dim = dim(td$mask_arr))
  roi[1:2, , ] <- TRUE

  il_wavelet <- encode(
    td$X,
    spec_space_wavelet_active(levels_space = 1, levels_time = 0),
    mask = td$mask_vol,
    materialize = "matrix"
  )
  il_st <- encode(
    td$X,
    spec_st(
      time = spec_time_bspline(k = 4, degree = 3),
      space = spec_space_hrbf(params = list(sigma0 = 4, levels = 2, seed = 42))
    ),
    mask = td$mask_vol,
    materialize = "matrix"
  )

  expect_false(is_explicit_latent(il_wavelet))
  expect_false(is_explicit_latent(il_st))
  expect_equal(reconstruct_matrix(il_wavelet), predict(il_wavelet))
  expect_equal(reconstruct_matrix(il_st), predict(il_st))
  expect_equal(
    reconstruct_matrix(il_st, time_idx = 2:5, roi_mask = roi),
    predict(il_st, time_idx = 2:5, roi_mask = roi)
  )
  expect_equal(dim(reconstruct_array(il_wavelet)), c(dim(td$mask_arr), td$n_time))
  expect_equal(dim(reconstruct_array(il_st, time_idx = 1:3)), c(dim(td$mask_arr), 3))
  expect_true(nzchar(latent_meta(il_wavelet)$family))
  expect_true(nzchar(latent_meta(il_st)$family))
})

test_that("handle-backed explicit spatial encoders support common interface", {
  skip_if_not_installed("rgsp")

  td <- make_test_data_ci_int(n_time = 6, dims = c(3, 3, 2))
  roi <- array(FALSE, dim = dim(td$mask_arr))
  roi[1:2, , ] <- TRUE

  lv_heat <- encode(
    td$X,
    spec_space_heat(scales = c(1, 2), order = 10, k_neighbors = 3),
    mask = td$mask_vol,
    materialize = "handle"
  )

  expect_true(is_explicit_latent(lv_heat))
  expect_true(inherits(lv_heat@loadings, "LoadingsHandle"))
  expect_equal(reconstruct_matrix(lv_heat), as.matrix(lv_heat))
  expect_equal(
    reconstruct_matrix(lv_heat, roi_mask = roi),
    as.matrix(lv_heat)[, which(as.logical(roi)), drop = FALSE]
  )
  expect_equal(dim(reconstruct_array(lv_heat)), c(dim(td$mask_arr), td$n_time))
})

test_that("handle-backed spatial encoders preserve k_neighbors on rematerialization", {
  skip_if_not_installed("rgsp")

  td <- make_test_data_ci_int(n_time = 6, dims = c(3, 3, 2))
  old_opt <- getOption("fmrilatent.registry.enabled")
  options(fmrilatent.registry.enabled = FALSE)
  on.exit(options(fmrilatent.registry.enabled = old_opt), add = TRUE)

  lv_heat_handle <- encode(
    td$X,
    spec_space_heat(scales = c(1, 2), order = 10, k_neighbors = 3),
    mask = td$mask_vol,
    materialize = "handle"
  )
  lv_heat_matrix <- encode(
    td$X,
    spec_space_heat(scales = c(1, 2), order = 10, k_neighbors = 3),
    mask = td$mask_vol,
    materialize = "matrix"
  )

  expect_equal(as.matrix(loadings(lv_heat_handle)), as.matrix(loadings(lv_heat_matrix)))
})

test_that("saved parcel templates replay identical encodings after reload", {
  td <- make_test_data_ci_int(n_time = 10, dims = c(4, 4, 2))
  map <- as.integer(rep_len(1:4, td$n_vox))
  reduction <- make_cluster_reduction(td$mask_vol, map)
  tmpl <- parcel_basis_template(reduction, basis_slepian(k = 2))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  save_template(tmpl, path)
  tmpl_loaded <- load_template(path)

  lv1 <- encode(td$X, spec_space_parcel(tmpl), mask = td$mask_vol)
  lv2 <- encode(td$X, spec_space_parcel(tmpl_loaded), mask = td$mask_vol)

  expect_equal(reconstruct_matrix(lv1), reconstruct_matrix(lv2))
  expect_equal(as.matrix(loadings(lv1)), as.matrix(loadings(lv2)))
  expect_equal(offset(lv1), offset(lv2))
  expect_equal(latent_meta(lv1)$label_map, latent_meta(lv2)$label_map)
})
