library(testthat)

test_that("compat profile helpers resolve expected values", {
  old <- options(fmrilatent.compat = NULL)
  on.exit(options(old), add = TRUE)

  expect_equal(fmrilatent_compat_profile(), "native")
  expect_equal(fmrilatent_compat_profile("neuroarchive_0.1.1"), "neuroarchive_0.1.1")
})

test_that("haar option aliases are honored in neuroarchive compat profile", {
  if (!exists("forward_lift_rcpp", mode = "function")) {
    skip("forward_lift_rcpp is not available")
  }

  old <- options(
    fmrilatent.compat = "native",
    fmrilatent.haar.use_rcpp = NULL,
    fmrilatent.hwt.use_rcpp = NULL,
    lna.hwt.use_rcpp = FALSE
  )
  on.exit(options(old), add = TRUE)

  native_rcpp <- fmrilatent:::use_haar_rcpp()
  options(fmrilatent.compat = "neuroarchive_0.1.1")
  compat_rcpp <- fmrilatent:::use_haar_rcpp()

  expect_true(native_rcpp)
  expect_false(compat_rcpp)
})

test_that("hrbf option alias is honored only in neuroarchive compat profile", {
  old <- options(
    fmrilatent.compat = "native",
    fmrilatent.hrbf.use_rcpp = NULL,
    lna.hrbf.use_rcpp_helpers = TRUE
  )
  on.exit(options(old), add = TRUE)

  expect_false(fmrilatent:::use_hrbf_rcpp())
  options(fmrilatent.compat = "neuroarchive_0.1.1")
  expect_true(fmrilatent:::use_hrbf_rcpp())
})

test_that("lna_hrbf_basis_from_params supports full-grid and active-grid outputs", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask_arr[c(1, 2, 3), c(1, 2, 3), 2] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  centres <- matrix(c(1, 1, 1,
                      2, 2, 2), ncol = 3, byrow = TRUE)
  sigmas <- c(1.25, 0.75)
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian")

  B_full <- lna_hrbf_basis_from_params(
    params = params,
    mask = mask,
    centres = centres,
    sigmas = sigmas,
    full_grid = TRUE
  )
  B_active <- lna_hrbf_basis_from_params(
    params = params,
    mask = mask,
    centres = centres,
    sigmas = sigmas,
    full_grid = FALSE
  )

  expect_equal(ncol(B_full), length(mask_arr))
  expect_equal(ncol(B_active), sum(mask_arr))
  expect_equal(nrow(B_full), nrow(B_active))
})

test_that("lna_hrbf_centres_from_params supports explicit and seeded centres", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask_arr[c(1, 2, 3), c(1, 2, 3), 2] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, seed = 5L)
  centres <- matrix(c(1, 1, 1, 2, 2, 2), ncol = 3, byrow = TRUE)
  sigmas <- c(1.25, 0.75)

  cs_explicit <- fmrilatent:::lna_hrbf_centres_from_params(
    params = params,
    mask = mask,
    centres = centres,
    sigmas = sigmas
  )
  expect_equal(cs_explicit$centres, centres)
  expect_equal(cs_explicit$sigmas, sigmas)
  expect_true(isTRUE(cs_explicit$centres_stored))

  cs_seeded <- fmrilatent:::lna_hrbf_centres_from_params(
    params = params,
    mask = mask,
    compat_profile = "neuroarchive_0.1.1"
  )
  expect_true(nrow(cs_seeded$centres) > 0)
  B_seeded <- lna_hrbf_basis_from_params(
    params = params,
    mask = mask,
    compat_profile = "neuroarchive_0.1.1"
  )
  expect_equal(nrow(B_seeded), nrow(cs_seeded$centres))
})

test_that("compat profile disables tiny-mask identity fallback", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 6, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 3L)

  B_native <- lna_hrbf_basis_from_params(
    params = params,
    mask = mask,
    full_grid = TRUE,
    compat_profile = "native"
  )
  B_compat <- lna_hrbf_basis_from_params(
    params = params,
    mask = mask,
    full_grid = TRUE,
    compat_profile = "neuroarchive_0.1.1"
  )

  expect_equal(dim(B_native), c(27, 27))
  expect_lt(nrow(B_compat), 27)
  expect_equal(ncol(B_compat), 27)
})

test_that("low-level haar matrix wrappers roundtrip in Morton space", {
  set.seed(10)
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(3 * sum(mask_arr)), nrow = 3)

  morton_idx <- fmrilatent:::get_morton_ordered_indices(mask_arr, 42L)
  mask_linear <- which(mask_arr)
  perm <- match(morton_idx, mask_linear)
  data_morton <- X[, perm, drop = FALSE]
  full_order <- fmrilatent:::get_morton_ordered_indices(array(TRUE, dim(mask_arr)), 42L)
  mask_flat_morton <- as.logical(mask_arr)[full_order]
  scalings <- fmrilatent:::precompute_haar_scalings(mask_arr, 1L)

  fw <- lna_forward_lift_matrix(
    data_morton = data_morton,
    mask_flat_morton = mask_flat_morton,
    mask_dims = dim(mask_arr),
    levels = 1L,
    scalings = scalings,
    compat_profile = "neuroarchive_0.1.1"
  )
  reco <- lna_inverse_lift_matrix(
    root_coeff = fw$root_coeff,
    detail_coeffs_by_level = fw$detail_coeffs_by_level,
    mask_flat_morton = mask_flat_morton,
    mask_dims = dim(mask_arr),
    levels = 1L,
    scalings = scalings,
    compat_profile = "neuroarchive_0.1.1"
  )

  expect_equal(reco, data_morton, tolerance = 1e-7)
})
