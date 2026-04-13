library(testthat)
library(Matrix)
library(neuroim2)

.make_mask_dsp <- function() {
  LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), NeuroSpace(c(2, 2, 1)))
}

.make_lvec_dsp <- function() {
  mask_vol <- .make_mask_dsp()
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

.make_implicit_dsp <- function() {
  mask_vol <- .make_mask_dsp()
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

.make_parcel_template_dsp <- function() {
  mask_vol <- .make_mask_dsp()
  loadings <- Matrix(matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE)
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_slepian(k = 1),
      center = FALSE,
      meta = list(family = "spec_slepian", k = 1L)
    ),
    class = "ParcelBasisTemplate"
  )
}

.make_hierarchical_template_dsp <- function() {
  mask_vol <- .make_mask_dsp()
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
      meta = list(family = "hierarchical_laplacian"))
}

.make_awpt_template_dsp <- function() {
  mask_vol <- .make_mask_dsp()
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  loadings <- diag(4)
  structure(
    list(
      loadings = Matrix(loadings, sparse = FALSE),
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      roughness = Matrix::Diagonal(x = rep(1, 4)),
      reduction = reduction,
      basis_spec = basis_awpt_wavelet(scales = rep(1, 4)),
      atoms = data.frame(
        col_id = seq_len(4),
        cluster_id = c(1L, 1L, 2L, 2L),
        scale = rep(1, 4),
        scale_index = seq_len(4),
        roughness_weight = rep(1, 4),
        stringsAsFactors = FALSE
      ),
      center = FALSE,
      meta = list(family = "awpt_wavelet", method = "AWPT", k = 4L)
    ),
    class = "AWPTBasisTemplate"
  )
}

.make_transport_template_dsp <- function() {
  mask_vol <- .make_mask_dsp()
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
      meta = list(family = "spec_slepian", k = 2L)
    ),
    class = "ParcelBasisTemplate"
  )
}

.make_transport_operator_dsp <- function() {
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  .linear_map_from_matrix(
    A,
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(operator = "toy_transport")
  )
}

test_that("LatentNeuroVec satisfies the domain/support protocol", {
  lv <- .make_lvec_dsp()
  mat <- reconstruct_matrix(lv)
  arr <- reconstruct_array(lv)

  expect_equal(latent_domain(lv), neuroim2::space(mask(lv)))
  expect_equal(as.array(latent_support(lv)), as.array(mask(lv)))
  expect_equal(wrap_decoded(lv, mat), arr)
})

test_that("ImplicitLatent satisfies the domain/support protocol", {
  il <- .make_implicit_dsp()
  mat <- reconstruct_matrix(il)
  arr <- reconstruct_array(il)

  expect_equal(latent_domain(il), neuroim2::space(mask(il)))
  expect_equal(as.array(latent_support(il)), as.array(mask(il)))
  expect_equal(wrap_decoded(il, mat), arr)
})

test_that("transport-backed latent objects expose domain/support without changing decoder numerics", {
  tmpl <- .make_transport_template_dsp()
  op <- .make_transport_operator_dsp()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Z <- matrix(
    c(1, 2,
      0, 1,
      -1, 3),
    nrow = 3, byrow = TRUE
  )
  tl <- encode_transport(
    x = Z %*% t(D),
    basis_asset = tmpl,
    observation_operator = op,
    mask = .make_mask_dsp(),
    center = FALSE,
    lambda = 0
  )

  expect_equal(latent_domain(tl), neuroim2::space(mask(tl)))
  expect_equal(as.array(latent_support(tl)), as.array(mask(tl)))
  expect_equal(wrap_decoded(tl, reconstruct_matrix(tl)), reconstruct_array(tl))
  expect_equal(
    decode_coefficients(tl, c(2, -1), space = "native"),
    as.vector(D %*% matrix(c(2, -1), ncol = 1)),
    tolerance = 1e-8
  )
})

test_that("template_support is consistent with template_mask for current volumetric templates", {
  parcel <- .make_parcel_template_dsp()
  hier <- .make_hierarchical_template_dsp()
  awpt <- .make_awpt_template_dsp()

  for (tmpl in list(parcel, hier, awpt)) {
    expect_equal(as.array(template_support(tmpl)), as.array(template_mask(tmpl)))
    expect_equal(template_measure(tmpl), rep(1, sum(as.array(template_mask(tmpl)))))
    expect_equal(template_domain(tmpl), neuroim2::space(template_mask(tmpl)))
  }
})

test_that("wrap_decoded rejects incompatible volumetric shapes", {
  lv <- .make_lvec_dsp()
  il <- .make_implicit_dsp()

  expect_error(wrap_decoded(lv, matrix(1, nrow = 1, ncol = 3)), "support cardinality")
  expect_error(wrap_decoded(il, matrix(1, nrow = 1, ncol = 3)), "support cardinality")
})
