# Tests for hierarchical_template.R
# Note: build_hierarchical_template has issues with small matrices (diag error)
# Tests that require building templates are skipped until underlying bug is fixed

# -----------------------------------------------------------------------------
# Helper functions for testing
# -----------------------------------------------------------------------------

create_test_parcellations <- function(n_vox) {
  coarse <- rep(c(1L, 2L), each = ceiling(n_vox / 2))[1:n_vox]
  fine <- rep(c(1L, 2L, 3L, 4L), each = ceiling(n_vox / 4))[1:n_vox]
  list(coarse, fine)
}

# -----------------------------------------------------------------------------
# Tests for is_hierarchical_template
# -----------------------------------------------------------------------------

test_that("is_hierarchical_template returns FALSE for non-template objects", {
  expect_false(is_hierarchical_template(NULL))
  expect_false(is_hierarchical_template(list()))
  expect_false(is_hierarchical_template("string"))
  expect_false(is_hierarchical_template(42))
  expect_false(is_hierarchical_template(matrix(1:4, 2, 2)))
})

# -----------------------------------------------------------------------------
# Tests for build_hierarchical_template - input validation
# -----------------------------------------------------------------------------

test_that("build_hierarchical_template rejects invalid mask", {
  parcellations <- list(c(1L, 1L, 2L, 2L))
  k_per_level <- c(2L)

  expect_error(
    build_hierarchical_template("not_a_mask", parcellations, k_per_level)
  )
})

test_that("build_hierarchical_template rejects non-list parcellations", {
  skip_if_not_installed("neuroim2")

  mask <- array(TRUE, dim = c(2, 2, 2))
  k_per_level <- c(2L)

  expect_error(
    build_hierarchical_template(mask, "not_a_list", k_per_level),
    "non-empty list"
  )

  expect_error(
    build_hierarchical_template(mask, c(1, 2, 3, 4), k_per_level),
    "non-empty list"
  )
})

test_that("build_hierarchical_template rejects empty parcellations list", {
  skip_if_not_installed("neuroim2")

  mask <- array(TRUE, dim = c(2, 2, 2))
  k_per_level <- c(2L)

  expect_error(
    build_hierarchical_template(mask, list(), k_per_level),
    "non-empty list"
  )
})

test_that("build_hierarchical_template rejects mismatched k_per_level length", {
  skip_if_not_installed("neuroim2")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(1L, n_vox), rep(c(1L, 2L), each = n_vox / 2))

  expect_error(
    build_hierarchical_template(mask, parcellations, c(2L, 2L, 2L)),
    "must match"
  )
})

test_that("build_hierarchical_template rejects non-positive scalar parameters", {
  skip_if_not_installed("neuroim2")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(1L, n_vox), rep(c(1L, 2L), each = n_vox / 2))

  expect_error(
    build_hierarchical_template(mask, parcellations, c(0L, 2L)),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_hierarchical_template(mask, parcellations, c(2L, 2L), k_neighbors = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_hierarchical_template(mask, parcellations, c(2L, 2L), ridge = NA_real_),
    class = "fmrilatent_error_invalid_scalar"
  )
})

test_that("build_hierarchical_template rejects parcellation with wrong voxel count", {
  skip_if_not_installed("neuroim2")

  mask <- array(TRUE, dim = c(2, 2, 2))  # 8 voxels
  parcellations <- list(c(1L, 1L, 2L, 2L))  # Only 4 elements
  k_per_level <- c(2L)

  expect_error(
    build_hierarchical_template(mask, parcellations, k_per_level),
    "has length.*but mask has"
  )
})

test_that("build_hierarchical_template uses base eigen when rank equals parcel size", {
  skip_if_not_installed("rgsp")

  mask <- array(TRUE, dim = c(2, 2, 1))
  parcellations <- list(rep(1L, sum(mask)))
  template <- build_hierarchical_template(mask, parcellations, k_per_level = sum(mask),
                                          k_neighbors = 2L)
  expect_s4_class(template, "HierarchicalBasisTemplate")
  expect_equal(dim(template_loadings(template)), c(sum(mask), sum(mask)))
})

test_that("build_hierarchical_template warns with classed condition on Cholesky fallback", {
  skip_if_not_installed("rgsp")
  skip_if_not(exists("local_mocked_bindings", envir = asNamespace("testthat")),
              "testthat local_mocked_bindings unavailable")

  local_mocked_bindings(
    Cholesky = function(...) stop("forced Cholesky failure"),
    .package = "Matrix"
  )

  mask <- array(TRUE, dim = c(2, 2, 1))
  parcellations <- list(rep(1L, sum(mask)))

  expect_warning(
    template <- build_hierarchical_template(
      mask,
      parcellations,
      k_per_level = sum(mask),
      k_neighbors = 2L,
      solver = "chol"
    ),
    class = "fmrilatent_warning_cholesky_fallback"
  )
  expect_s4_class(template, "HierarchicalBasisTemplate")
})

# -----------------------------------------------------------------------------
# Tests for build_hierarchical_template - successful construction
# -----------------------------------------------------------------------------

test_that("hierarchical Laplacian eigensolver returns finite low modes", {
  skip_if_not_installed("RSpectra")

  L <- Matrix::Matrix(
    matrix(
      c(1, -1, 0, 0,
        -1, 2, -1, 0,
        0, -1, 2, -1,
        0, 0, -1, 1),
      nrow = 4L,
      byrow = TRUE
    ),
    sparse = TRUE
  )

  eig <- fmrilatent:::.hierarchical_laplacian_eigs(L, k = 2L)

  expect_equal(dim(eig$vectors), c(4L, 2L))
  expect_true(all(is.finite(eig$values)))
  expect_equal(eig$values, sort(eig$values), tolerance = 1e-8)
  expect_lt(abs(eig$values[[1L]]), 1e-6)
})

test_that("build_hierarchical_template creates valid template with single level", {
    skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(c(1L, 2L), each = n_vox / 2))
  k_per_level <- c(2L)

  template <- build_hierarchical_template(
    mask, parcellations, k_per_level,
    k_neighbors = 3L
  )

  expect_true(is_hierarchical_template(template))
  expect_s4_class(template, "HierarchicalBasisTemplate")
})

test_that("HierarchicalBasisTemplate validity enforces cross-slot dimensions", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  loadings <- Matrix::Matrix(
    matrix(
      c(1, 0,
        0, 1,
        1, 1,
        0, 1),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = FALSE
  )
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE, LDL = TRUE)

  expect_error(
    new("HierarchicalBasisTemplate",
      mask = mask_vol,
      space = neuroim2::NeuroSpace(c(2, 2, 1, 1)),
      levels = list(rep(1L, 4)),
      parents = list(integer(0)),
      loadings = loadings,
      gram_factor = gram_factor,
      atoms = data.frame(col_id = 1L, level = 1L, parcel_id = 1L),
      meta = list(family = "hierarchical_laplacian")
    ),
    "@atoms rows must match"
  )

  expect_error(
    new("HierarchicalBasisTemplate",
      mask = mask_vol,
      space = neuroim2::NeuroSpace(c(2, 2, 1, 1)),
      levels = list(rep(1L, 3)),
      parents = list(integer(0)),
      loadings = loadings,
      gram_factor = gram_factor,
      atoms = data.frame(col_id = 1:2, level = 1L, parcel_id = 1L),
      meta = list(family = "hierarchical_laplacian")
    ),
    "@levels\\[\\[1\\]\\] length must match"
  )
})

test_that("build_hierarchical_template creates valid template with two levels", {
    skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  # Use larger mask to avoid RSpectra dimension requirements (min 3 per parcel)
  mask <- array(TRUE, dim = c(4, 4, 4))
  n_vox <- sum(mask)
  parcellations <- create_test_parcellations(n_vox)
  k_per_level <- c(2L, 1L)

  template <- build_hierarchical_template(
    mask, parcellations, k_per_level,
    k_neighbors = 6L
  )

  expect_true(is_hierarchical_template(template))
  expect_length(template@levels, 2)
})

test_that("build_hierarchical_template clamps k_neighbors for two-voxel parcels", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")

  mask <- array(TRUE, dim = c(2, 1, 1))
  template <- build_hierarchical_template(
    mask = mask,
    parcellations = list(c(1L, 1L)),
    k_per_level = 1L,
    k_neighbors = 6L
  )

  expect_true(is_hierarchical_template(template))
  expect_equal(nrow(template_loadings(template)), 2L)
})

# -----------------------------------------------------------------------------
# Tests for encode_hierarchical
# -----------------------------------------------------------------------------

test_that("encode_hierarchical rejects non-template input", {
  X <- matrix(rnorm(10), nrow = 2)

  expect_error(
    encode_hierarchical(X, "not_a_template"),
    "must be a HierarchicalBasisTemplate"
  )

  expect_error(
    encode_hierarchical(X, list()),
    "must be a HierarchicalBasisTemplate"
  )
})

test_that("encode_hierarchical returns LatentNeuroVec with correct dimensions", {
    skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(c(1L, 2L), each = n_vox / 2))
  k_per_level <- c(2L)

  template <- build_hierarchical_template(
    mask, parcellations, k_per_level,
    k_neighbors = 3L
  )

  n_time <- 10
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  result <- encode_hierarchical(X, template)

  expect_s4_class(result, "LatentNeuroVec")
})

test_that("encode_hierarchical preserves custom mask geometry", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  spc3 <- NeuroSpace(c(2, 2, 1), spacing = c(2, 3, 4), origin = c(10, 20, 30))
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), spc3)
  parcellations <- list(c(1L, 1L, 2L, 2L))
  template <- build_hierarchical_template(mask_vol, parcellations, k_per_level = 1L, k_neighbors = 2L)
  X <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)

  result <- encode_hierarchical(X, template, mask = mask_vol)

  expect_equal(neuroim2::spacing(neuroim2::space(mask(result))), c(2, 3, 4))
  expect_equal(neuroim2::origin(neuroim2::space(mask(result))), c(10, 20, 30))
})

test_that("hierarchical template encoding rejects mismatched masks", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask1_arr <- array(c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE), dim = c(2, 2, 2))
  mask2_arr <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 2))
  mask1 <- LogicalNeuroVol(mask1_arr, NeuroSpace(dim(mask1_arr)))
  mask2 <- LogicalNeuroVol(mask2_arr, NeuroSpace(dim(mask2_arr)))
  template <- build_hierarchical_template(mask1, list(rep(1L, sum(mask1_arr))), k_per_level = 1L, k_neighbors = 2L)
  X <- matrix(rnorm(5 * sum(mask1_arr)), nrow = 5)

  expect_error(
    encode(X, spec_hierarchical_template(template), mask = mask2, materialize = "matrix"),
    "template mask"
  )
})

# -----------------------------------------------------------------------------
# Tests for project_hierarchical
# -----------------------------------------------------------------------------

test_that("project_hierarchical rejects non-template input", {
  X <- matrix(rnorm(10), nrow = 2)

  expect_error(
    project_hierarchical("not_a_template", X),
    "must be a HierarchicalBasisTemplate"
  )
})

test_that("project_hierarchical returns coefficient matrix", {
  skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(c(1L, 2L), each = n_vox / 2))
  k_per_level <- c(2L)

  template <- build_hierarchical_template(
    mask, parcellations, k_per_level,
    k_neighbors = 3L
  )

  n_time <- 10
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  coeff <- project_hierarchical(template, X)

  expect_s4_class(coeff, "Matrix")
  expect_equal(nrow(coeff), n_time)
})

test_that("project_hierarchical returns coefficients compatible with stored analysis loadings", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  loadings <- Matrix::Matrix(
    matrix(
      c(1, 0,
        1, 1,
        0, 1,
        1, 0),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = FALSE
  )
  template <- new("HierarchicalBasisTemplate",
    mask = mask_vol,
    space = neuroim2::NeuroSpace(c(2, 2, 1, 1)),
    levels = list(rep(1L, 4)),
    parents = list(integer(0)),
    loadings = loadings,
    gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE, LDL = TRUE),
    atoms = data.frame(
      col_id = 1:2,
      level = 1L,
      parcel_id = 1L,
      parent_id = NA_integer_,
      mode = 1:2
    ),
    meta = list(family = "hierarchical_laplacian")
  )
  X <- matrix(
    c(1, 2, 3, 4,
      2, 3, 4, 5),
    nrow = 2,
    byrow = TRUE
  )

  coeff <- project_hierarchical(template, X)
  encoded <- encode_hierarchical(X, template, materialize = "matrix")

  expect_equal(as.matrix(coeff), as.matrix(basis(encoded)), tolerance = 1e-8)
  expect_equal(
    as.matrix(coeff) %*% t(as.matrix(loadings(encoded))),
    as.matrix(encoded),
    tolerance = 1e-8
  )
})

# -----------------------------------------------------------------------------
# Tests for save_hierarchical_template and load_hierarchical_template
# -----------------------------------------------------------------------------

test_that("save_hierarchical_template rejects non-template input", {
  expect_error(
    save_hierarchical_template("not_a_template", tempfile()),
    "must be a HierarchicalBasisTemplate"
  )
})

test_that("load_hierarchical_template rejects invalid RDS", {
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  saveRDS(list(a = 1, b = 2), temp_file)

  expect_error(
    load_hierarchical_template(temp_file),
    "does not contain a HierarchicalBasisTemplate"
  )
})

test_that("save and load hierarchical template round-trip", {
    skip_if_not_installed("neuroim2")
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask)
  parcellations <- list(rep(c(1L, 2L), each = n_vox / 2))
  k_per_level <- c(2L)

  template <- build_hierarchical_template(
    mask, parcellations, k_per_level,
    k_neighbors = 3L, label = "test_template"
  )

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  result_path <- save_hierarchical_template(template, temp_file)
  expect_true(file.exists(temp_file))

  loaded <- load_hierarchical_template(temp_file)
  expect_true(is_hierarchical_template(loaded))
  expect_equal(loaded@meta$label, "test_template")
})
