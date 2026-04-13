library(testthat)
library(Matrix)
library(neuroim2)

make_test_mask_vol_ci <- function() {
  LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), NeuroSpace(c(2, 2, 1)))
}

make_test_lvec_ci <- function() {
  mask_vol <- make_test_mask_vol_ci()
  basis <- Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2), sparse = FALSE)
  loadings <- Matrix(matrix(c(1, 0,
                              0, 1,
                              1, 1,
                              0, 0), nrow = 4, byrow = TRUE), sparse = FALSE)
  spc <- NeuroSpace(c(2, 2, 1, 3))
  LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = spc,
    mask = mask_vol,
    offset = c(0.5, 1, -0.5, 0),
    meta = list(family = "explicit_test")
  )
}

make_test_implicit_ci <- function() {
  mask_vol <- make_test_mask_vol_ci()
  mask_arr <- as.array(mask_vol)
  full <- matrix(seq_len(12), nrow = 3, ncol = 4)
  implicit_latent(
    coeff = list(full = full),
    decoder = function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL, ...) {
      rec <- full
      if (!is.null(time_idx)) {
        rec <- rec[time_idx, , drop = FALSE]
      }
      roi_subset_columns(rec, mask_arr, roi_mask)
    },
    meta = list(family = "implicit_test"),
    mask = mask_arr
  )
}

make_test_parcel_template_ci <- function(center = FALSE) {
  mask_vol <- make_test_mask_vol_ci()
  loadings <- Matrix(matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE)
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_slepian(k = 1),
      center = center,
      meta = list(
        family = "spec_slepian",
        k = 1L,
        ridge = 1e-8,
        label_map = list(A = 1L),
        cluster_map = list(`1` = 1:2, `2` = 3:4)
      )
    ),
    class = "ParcelBasisTemplate"
  )
}

make_test_hierarchical_template_ci <- function() {
  mask_vol <- make_test_mask_vol_ci()
  loadings <- Matrix(matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE, LDL = TRUE)
  new("HierarchicalBasisTemplate",
      mask = mask_vol,
      space = NeuroSpace(c(2, 2, 1, 1)),
      levels = list(rep(1L, 4)),
      parents = list(integer(0)),
      loadings = loadings,
      gram_factor = gram_factor,
      atoms = data.frame(
        col_id = 1L,
        level = 1L,
        parcel_id = 1L,
        parent_id = NA_integer_,
        mode = 1L,
        label = "L1:P1"
      ),
      meta = list(family = "hierarchical_laplacian", label = "test_hier"))
}

make_test_parcel_template_nonorthogonal_ci <- function() {
  mask_vol <- make_test_mask_vol_ci()
  loadings <- Matrix(matrix(
    c(1, 0,
      0, 1,
      1, 1,
      0, 1),
    nrow = 4, byrow = TRUE
  ), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE)
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_slepian(k = 2),
      center = FALSE,
      meta = list(
        family = "spec_slepian",
        k = 2L,
        ridge = 1e-8,
        label_map = list(A = 1L, B = 2L),
        cluster_map = list(`1` = 1:2, `2` = 3:4)
      )
    ),
    class = "ParcelBasisTemplate"
  )
}

test_that("common latent interface works for LatentNeuroVec", {
  lv <- make_test_lvec_ci()
  roi <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
  mat <- as.matrix(lv)
  arr <- reconstruct_array(lv)
  arr_roi <- reconstruct_array(lv, roi_mask = roi)
  arr_time <- reconstruct_array(lv, time_idx = 2:3)

  expect_true(is_explicit_latent(lv))
  expect_equal(latent_meta(lv)$family, "explicit_test")
  expect_equal(reconstruct_matrix(lv), mat)
  expect_equal(reconstruct_matrix(lv, time_idx = 2:3), mat[2:3, , drop = FALSE])
  expect_equal(reconstruct_matrix(lv, roi_mask = roi), mat[, c(1, 3), drop = FALSE])
  expect_equal(
    reconstruct_matrix(lv, time_idx = 2:3, roi_mask = roi),
    mat[2:3, c(1, 3), drop = FALSE]
  )
  expect_equal(dim(arr), c(2, 2, 1, 3))
  expect_equal(arr[, , 1, 1], matrix(c(mat[1, 1], mat[1, 2], mat[1, 3], mat[1, 4]), nrow = 2))
  expect_equal(arr_time[, , 1, 1], arr[, , 1, 2])
  expect_equal(arr_roi[, , 1, 1], matrix(c(mat[1, 1], 0, mat[1, 3], 0), nrow = 2))
})

test_that("common latent interface works for ImplicitLatent", {
  il <- make_test_implicit_ci()
  roi <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
  mat <- as.matrix(il)
  arr <- reconstruct_array(il)
  arr_roi <- reconstruct_array(il, roi_mask = roi)
  arr_time <- reconstruct_array(il, time_idx = 2:3)

  expect_false(is_explicit_latent(il))
  expect_equal(latent_meta(il)$family, "implicit_test")
  expect_s4_class(mask(il), "LogicalNeuroVol")
  expect_equal(dim(as.array(mask(il))), c(2, 2, 1))
  expect_true(all(as.logical(as.array(mask(il)))))
  expect_equal(reconstruct_matrix(il), mat)
  expect_equal(reconstruct_matrix(il, time_idx = 2:3), mat[2:3, , drop = FALSE])
  expect_equal(reconstruct_matrix(il, roi_mask = roi), mat[, c(1, 3), drop = FALSE])
  expect_equal(
    reconstruct_matrix(il, time_idx = 2:3, roi_mask = roi),
    mat[2:3, c(1, 3), drop = FALSE]
  )
  expect_equal(dim(arr), c(2, 2, 1, 3))
  expect_equal(arr[, , 1, 1], matrix(c(mat[1, 1], mat[1, 2], mat[1, 3], mat[1, 4]), nrow = 2))
  expect_equal(arr_time[, , 1, 1], arr[, , 1, 2])
  expect_equal(arr_roi[, , 1, 1], matrix(c(mat[1, 1], 0, mat[1, 3], 0), nrow = 2))
})

test_that("common latent interface rejects malformed roi_mask inputs", {
  lv <- make_test_lvec_ci()
  il <- make_test_implicit_ci()
  bad_roi <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), dim = c(3, 2, 1))

  expect_error(reconstruct_matrix(lv, roi_mask = bad_roi), "roi_mask dimensions")
  expect_error(reconstruct_array(lv, roi_mask = bad_roi), "roi_mask dimensions")
  expect_error(reconstruct_matrix(il, roi_mask = bad_roi), "roi_mask dimensions")
  expect_error(reconstruct_array(il, roi_mask = bad_roi), "roi_mask dimensions")
})

test_that("template protocol works for parcel templates", {
  tmpl <- make_test_parcel_template_ci(center = TRUE)
  X <- matrix(c(1, 2, 3, 4,
                5, 6, 7, 8), nrow = 2, byrow = TRUE)
  proj <- template_project(tmpl, X)
  centered <- sweep(X, 2, colMeans(X), "-")
  expected_coeff <- centered[, 1, drop = FALSE]

  expect_true(is_template(tmpl))
  expect_equal(nrow(template_loadings(tmpl)), 4L)
  expect_s4_class(template_mask(tmpl), "LogicalNeuroVol")
  expect_equal(template_meta(tmpl)$family, "spec_slepian")
  expect_equal(proj$offset, colMeans(X))
  expect_equal(dim(proj$coefficients), c(2, 1))
  expect_equal(as.matrix(proj$coefficients), expected_coeff)
})

test_that("template protocol works for hierarchical templates", {
  tmpl <- make_test_hierarchical_template_ci()
  X <- matrix(c(1, 2, 3, 4,
                5, 6, 7, 8), nrow = 2, byrow = TRUE)
  proj <- template_project(tmpl, X)
  expected_coeff <- X[, 1, drop = FALSE]

  expect_true(is_template(tmpl))
  expect_equal(nrow(template_loadings(tmpl)), 4L)
  expect_s4_class(template_mask(tmpl), "LogicalNeuroVol")
  expect_equal(template_meta(tmpl)$family, "hierarchical_laplacian")
  expect_equal(proj$offset, numeric(0))
  expect_equal(dim(proj$coefficients), c(2, 1))
  expect_equal(as.matrix(proj$coefficients), expected_coeff)
})

test_that("template-backed explicit latent objects retain shared-basis coordinate semantics", {
  tmpl <- make_test_parcel_template_nonorthogonal_ci()
  payload <- .template_coordinate_payload(
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    default_measure = "unit"
  )
  Z_analysis <- matrix(
    c(1, 2,
      0, 1,
      -1, 3),
    nrow = 3, byrow = TRUE
  )
  X <- Z_analysis %*% t(payload$analysis_loadings)
  lv <- encode(X, spec_space_parcel(tmpl), mask = template_mask(tmpl))

  expect_equal(coef_time(lv, "analysis"), Z_analysis, tolerance = 1e-8)
  expect_equal(
    coef_time(lv, "raw"),
    t(payload$analysis_transform$to_raw(t(Z_analysis))),
    tolerance = 1e-8
  )
  expect_equal(coef_metric(lv, "raw"), payload$raw_metric, tolerance = 1e-8)
  expect_s3_class(basis_asset(lv), "ParcelBasisTemplate")
  expect_equal(as.matrix(loadings(lv)), as.matrix(payload$analysis_loadings), tolerance = 1e-8)
})

test_that("save_template and load_template round-trip parcel templates", {
  tmpl <- make_test_parcel_template_ci(center = FALSE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  save_template(tmpl, path)
  loaded <- load_template(path)

  expect_s3_class(loaded, "ParcelBasisTemplate")
  expect_equal(template_meta(loaded)$family, "spec_slepian")
  expect_equal(template_meta(loaded)$label_map, template_meta(tmpl)$label_map)
  expect_equal(as.matrix(template_loadings(loaded)), as.matrix(template_loadings(tmpl)))
})

test_that("save_template and load_template round-trip hierarchical templates", {
  tmpl <- make_test_hierarchical_template_ci()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  save_template(tmpl, path)
  loaded <- load_template(path)

  expect_s4_class(loaded, "HierarchicalBasisTemplate")
  expect_equal(template_meta(loaded)$family, "hierarchical_laplacian")
  expect_equal(as.matrix(template_loadings(loaded)), as.matrix(template_loadings(tmpl)))
})

test_that("load_template rejects invalid RDS payloads", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(a = 1), path)

  expect_error(load_template(path), "supported template object")
})

test_that("template generics reject unsupported objects", {
  expect_error(template_loadings(list()), "unable to find an inherited method")
  expect_error(template_mask(list()), "unable to find an inherited method")
  expect_error(template_meta(list()), "unable to find an inherited method")
  expect_error(template_project(list(), matrix(1, 1, 1)), "unable to find an inherited method")
  expect_error(save_template(list(), tempfile(fileext = '.rds')), "unable to find an inherited method")
})

test_that("is_template only recognizes supported template objects", {
  expect_false(is_template(list()))
  expect_false(is_template(make_test_lvec_ci()))
  expect_false(is_template(make_test_implicit_ci()))
  expect_true(is_template(make_test_parcel_template_ci()))
  expect_true(is_template(make_test_hierarchical_template_ci()))
})
