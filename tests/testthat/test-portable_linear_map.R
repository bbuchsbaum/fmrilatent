# Shared fixture: a transport-backed latent built from a parcel template and
# an identity observation operator. Defined at the top so it is visible when
# tests below run (testthat evaluates test_that() eagerly).
.portable_linear_map_td <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dim3 <- c(3L, 3L, 2L)
    mask_vol <- neuroim2::LogicalNeuroVol(array(TRUE, dim = dim3),
                                          neuroim2::NeuroSpace(dim3))
    V <- prod(dim3)
    parcels <- rep(seq_len(3L), length.out = V)
    red <- make_cluster_reduction(mask_vol, parcels)
    tpl <- parcel_basis_template(red, basis_slepian(k = 2))
    k <- template_rank(tpl)
    I_op <- diag(V)
    obs <- list(
      n_source = V, n_target = V,
      forward = function(x, ...) I_op %*% x,
      adjoint_apply = function(y, ...) t(I_op) %*% y,
      target_support = mask_vol,
      source_domain_id = "template",
      target_domain_id = "native"
    )
    T <- 8L
    Y <- matrix(rnorm(T * V), nrow = T, ncol = V)
    lat <- encode_transport(Y, basis_asset = tpl, field_operator = obs,
                            lambda = 1e-3, center = TRUE)
    cache <<- list(Y = Y, mask_vol = mask_vol, V = V, k = k, lat = lat)
    cache
  }
})

test_that("as_portable_linear_map coerces a base matrix", {
  M <- matrix(seq_len(6), nrow = 2, ncol = 3)
  op <- as_portable_linear_map(M)
  expect_identical(op$n_source, 3L)
  expect_identical(op$n_target, 2L)
  expect_true(is.function(op$forward))
  expect_true(is.function(op$adjoint_apply))
  expect_equal(op$forward(c(1, 0, 0)), M[, 1])
  expect_equal(op$adjoint_apply(c(1, 0)), M[1, ])
  expect_equal(op$materialize(), M)
  expect_identical(op$adjoint_convention, "euclidean_discrete")
})

test_that("as_portable_linear_map threads matrix-side metadata", {
  M <- matrix(runif(12), nrow = 4, ncol = 3)
  mask_vol <- neuroim2::LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)),
                                        neuroim2::NeuroSpace(c(2, 2, 1)))
  op <- as_portable_linear_map(
    M,
    source_domain_id = "template",
    target_domain_id = "native",
    source_support = 1:3,
    target_support = mask_vol,
    provenance = list(tag = "test")
  )
  expect_identical(op$source_domain_id, "template")
  expect_identical(op$target_domain_id, "native")
  expect_identical(op$source_support, 1:3)
  expect_true(inherits(op$target_support, "LogicalNeuroVol"))
  expect_identical(op$provenance$tag, "test")
})

test_that("as_portable_linear_map accepts callback lists with adjoint alias", {
  M <- matrix(rnorm(6), nrow = 3, ncol = 2)
  cbk <- list(
    n_source = 2L,
    n_target = 3L,
    forward = function(x, ...) M %*% x,
    adjoint = function(y, ...) t(M) %*% y
  )
  op <- as_portable_linear_map(cbk)
  expect_true(is.function(op$adjoint_apply))
  expect_equal(as.vector(op$adjoint_apply(c(1, 0, 0))), M[1, ])
})

test_that("source_support / target_support are top-level first-class fields", {
  M <- matrix(seq_len(6), nrow = 3, ncol = 2)
  mask_vol <- neuroim2::LogicalNeuroVol(array(TRUE, dim = c(3, 1, 1)),
                                        neuroim2::NeuroSpace(c(3, 1, 1)))
  cbk <- list(
    n_source = 2L,
    n_target = 3L,
    forward = function(x, ...) M %*% x,
    adjoint_apply = function(y, ...) t(M) %*% y,
    source_support = c("a", "b"),
    target_support = mask_vol
  )
  op <- as_portable_linear_map(cbk)
  expect_identical(op$source_support, c("a", "b"))
  expect_true(inherits(op$target_support, "LogicalNeuroVol"))
})

test_that("source_support / target_support fall back to provenance", {
  M <- matrix(seq_len(6), nrow = 3, ncol = 2)
  cbk <- list(
    n_source = 2L,
    n_target = 3L,
    forward = function(x, ...) M %*% x,
    adjoint_apply = function(y, ...) t(M) %*% y,
    provenance = list(source_support = 1:2, target_support = 1:3)
  )
  op <- as_portable_linear_map(cbk)
  expect_identical(op$source_support, 1:2)
  expect_identical(op$target_support, 1:3)
})

test_that("top-level fields win over provenance on collision", {
  M <- matrix(seq_len(6), nrow = 3, ncol = 2)
  cbk <- list(
    n_source = 2L,
    n_target = 3L,
    forward = function(x, ...) M %*% x,
    adjoint_apply = function(y, ...) t(M) %*% y,
    target_support = 100:102,
    provenance = list(target_support = 1:3)
  )
  op <- as_portable_linear_map(cbk)
  expect_identical(op$target_support, 100:102)
})

test_that("as_portable_linear_map rejects missing required fields", {
  bad <- list(forward = function(x) x)
  expect_error(as_portable_linear_map(bad), "required fields")
})

test_that("validate_portable_linear_map returns normalized form on success", {
  res <- validate_portable_linear_map(matrix(seq_len(4), 2, 2))
  expect_true(is.list(res))
  expect_identical(res$n_source, 2L)
  expect_identical(res$n_target, 2L)
  expect_true(is.function(res$forward))
})

test_that("validate_portable_linear_map(error = FALSE) returns FALSE on failure", {
  res <- validate_portable_linear_map(list(forward = function(x) x),
                                      error = FALSE)
  expect_false(isTRUE(res))
  expect_identical(res, FALSE)
})

test_that("validate_portable_linear_map errors by default on failure", {
  expect_error(validate_portable_linear_map(list(forward = function(x) x)),
               "required fields")
})

test_that("composed portable maps preserve outer supports and satisfy adjoint identity", {
  first <- list(
    n_source = 3L, n_target = 4L,
    forward = function(x, ...) matrix(seq_len(12), 4, 3) %*% x,
    adjoint_apply = function(y, ...) t(matrix(seq_len(12), 4, 3)) %*% y,
    source_support = 1:3,
    source_domain_id = "template",
    target_domain_id = "mid"
  )
  second <- list(
    n_source = 4L, n_target = 2L,
    forward = function(x, ...) matrix(seq_len(8), 2, 4) %*% x,
    adjoint_apply = function(y, ...) t(matrix(seq_len(8), 2, 4)) %*% y,
    target_support = c(10L, 20L),
    source_domain_id = "mid",
    target_domain_id = "native"
  )
  composed <- fmrilatent:::.compose_linear_maps(first, second)
  expect_identical(composed$n_source, 3L)
  expect_identical(composed$n_target, 2L)
  expect_identical(composed$source_support, 1:3)
  expect_identical(composed$target_support, c(10L, 20L))

  set.seed(1)
  xv <- rnorm(3)
  yv <- rnorm(2)
  lhs <- sum(composed$forward(xv) * yv)
  rhs <- sum(xv * composed$adjoint_apply(yv))
  expect_equal(lhs, rhs, tolerance = 1e-10)

  M1 <- matrix(seq_len(12), 4, 3)
  M2 <- matrix(seq_len(8), 2, 4)
  expect_equal(as.matrix(composed$materialize()), M2 %*% M1)
})

test_that("decode_coefficients wrap='none' returns raw numeric vector", {
  td <- .portable_linear_map_td()
  gamma <- rnorm(ncol(td$Y))
  gamma <- gamma[seq_len(td$k)]
  vals <- decode_coefficients(td$lat, gamma, space = "native", wrap = "none")
  expect_true(is.numeric(vals))
  expect_length(vals, td$V)
})

test_that("decode_coefficients wrap='auto' returns a NeuroVol for volumetric targets", {
  td <- .portable_linear_map_td()
  gamma <- seq_len(td$k) / td$k
  wrapped <- decode_coefficients(td$lat, gamma, space = "native", wrap = "auto")
  expect_true(methods::is(wrapped, "NeuroVol"))
  expect_equal(dim(wrapped), dim(td$mask_vol))
})

test_that("decode_coefficients wrap='auto' on template space errors cleanly", {
  td <- .portable_linear_map_td()
  gamma <- seq_len(td$k) / td$k
  expect_error(
    decode_coefficients(td$lat, gamma, space = "template", wrap = "auto"),
    "only implemented for space = \"native\""
  )
})

test_that("decode_coefficients wrap='auto' handles multi-column gamma", {
  td <- .portable_linear_map_td()
  gamma <- matrix(rnorm(td$k * 2), nrow = td$k, ncol = 2)
  wrapped <- decode_coefficients(td$lat, gamma, space = "native", wrap = "auto")
  expect_true(methods::is(wrapped, "NeuroVec") || methods::is(wrapped, "NeuroVol"))
  raw <- decode_coefficients(td$lat, gamma, space = "native", wrap = "none")
  expect_true(is.matrix(raw))
  expect_identical(dim(raw), c(as.integer(td$V), 2L))
})

