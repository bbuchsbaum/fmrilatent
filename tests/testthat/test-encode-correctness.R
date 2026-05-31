test_that("DCT norm none uses least-squares coefficients", {
  set.seed(101)
  n_time <- 8L
  n_vox <- 5L
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  mask <- array(TRUE, dim = c(1, n_vox, 1))

  lv <- encode(
    X,
    spec_time_dct(k = n_time, norm = "none"),
    mask = mask,
    materialize = "matrix"
  )

  expect_equal(as.matrix(lv), X, tolerance = 1e-10)
})

test_that("spec_st roundtrips with non-orthonormal DCT time basis", {
  set.seed(106)
  n_time <- 5L
  mask <- array(TRUE, dim = c(2, 2, 1))
  X <- matrix(rnorm(n_time * sum(mask)), nrow = n_time)

  spec <- spec_st(
    time = spec_time_dct(k = n_time, norm = "none"),
    space = spec_space_hrbf(params = list(
      sigma0 = 2,
      levels = 1L,
      radius_factor = 2.5,
      kernel_type = "gaussian",
      seed = 42L
    ))
  )

  lv <- encode(X, spec, mask = mask, materialize = "matrix")
  B_t <- lv$coeff$B_t

  expect_s3_class(lv, "ImplicitLatent")
  expect_gt(max(abs(crossprod(B_t) - diag(n_time))), 1)
  expect_equal(predict(lv), X, tolerance = 1e-10)
})

test_that("non-orthonormal B-spline encode matches QR projection", {
  set.seed(102)
  n_time <- 12L
  n_vox <- 4L
  k <- 5L
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  mask <- array(TRUE, dim = c(1, n_vox, 1))
  B <- as.matrix(build_bspline_basis(
    n_time, k,
    include_intercept = TRUE,
    orthonormalize = FALSE
  ))

  lv <- encode(
    X,
    spec_time_bspline(k = k, include_intercept = TRUE, orthonormalize = FALSE),
    mask = mask,
    materialize = "matrix"
  )
  expected <- B %*% qr.coef(qr(B), X)

  expect_equal(as.matrix(lv), expected, tolerance = 1e-8)
})

test_that("slepian_temporal_latent centers before storing offset", {
  set.seed(103)
  n_time <- 8L
  n_vox <- 3L
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  mask <- array(TRUE, dim = c(1, n_vox, 1))

  lv <- slepian_temporal_latent(
    X,
    mask,
    tr = 1,
    bandwidth = 0.49,
    k = n_time,
    denoise = FALSE,
    backend = "tridiag"
  )

  expect_equal(offset(lv), colMeans(X), tolerance = 1e-12)
  expect_equal(as.matrix(lv), X, tolerance = 1e-10)
})

test_that("handle registry rejects same id with different explicit content", {
  fmrilatent_registry_clear()
  on.exit(fmrilatent_registry_clear(), add = TRUE)

  h1 <- methods::new(
    "BasisHandle",
    id = "collision-test",
    dim = as.integer(c(2, 2)),
    kind = "explicit",
    spec = list(matrix = diag(2)),
    label = "first"
  )
  h2 <- methods::new(
    "BasisHandle",
    id = "collision-test",
    dim = as.integer(c(2, 2)),
    kind = "explicit",
    spec = list(matrix = 2 * diag(2)),
    label = "second"
  )

  expect_equal(as.matrix(basis_mat(h1)), diag(2))
  expect_error(
    basis_mat(h2),
    "does not match the requested handle fingerprint"
  )
})

test_that("robust Gram fallback emits structured warnings", {
  expect_warning(
    expect_warning(
      fmrilatent:::.robust_gram_solve(matrix(0, 2, 2), matrix(1, 2, 1), ridge = 0),
      class = "fmrilatent_warning_gram_svd_fallback"
    ),
    class = "fmrilatent_warning_gram_rank_collapse"
  )
})

test_that("explicit spatial encoders honor materialize handle", {
  set.seed(104)
  mask <- array(TRUE, dim = c(2, 2, 1))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  spec <- spec_space_hrbf(params = list(
    sigma0 = 2,
    levels = 0,
    radius_factor = 2,
    kernel_type = "gaussian",
    seed = 1
  ))

  as_matrix <- encode(X, spec, mask = mask, materialize = "matrix")
  as_handle <- encode(X, spec, mask = mask, materialize = "handle")

  expect_false(inherits(as_matrix@loadings, "LoadingsHandle"))
  expect_true(inherits(as_handle@loadings, "LoadingsHandle"))
  expect_equal(as.matrix(as_handle), as.matrix(as_matrix), tolerance = 1e-10)
})

test_that("latent_factory accepts canonical registry family names", {
  set.seed(105)
  mask <- array(TRUE, dim = c(1, 4, 1))
  X <- matrix(rnorm(8 * sum(mask)), nrow = 8)

  canonical <- latent_factory("time_dct", x = X, mask = mask, k = 4, materialize = "matrix")
  legacy <- latent_factory("dct_time", x = X, mask = mask, k = 4, materialize = "matrix")

  expect_equal(canonical@meta$family, "time_dct")
  expect_equal(as.matrix(canonical), as.matrix(legacy), tolerance = 1e-12)
})
