# Tests for saveRDS/readRDS roundtrip of LatentNeuroVec, especially with
# BasisHandle and LoadingsHandle slots whose materialization cache lives in
# .fmrilatent_cache_env. A handle that survives serialization must be able to
# rematerialize from its spec, even when the in-memory registry is empty or
# disabled.

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

skip_if_no_neuroim2 <- function() {
  testthat::skip_if_not_installed("neuroim2")
}

make_small_mask <- function(dims = c(3L, 3L, 3L)) {
  arr <- array(TRUE, dim = dims)
  neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims))
}

make_small_space <- function(spatial = c(3L, 3L, 3L), n_time = 8L) {
  neuroim2::NeuroSpace(c(spatial, n_time))
}

# Round-trip an object through tempfile to genuinely exercise R's serialization
# path. Returns the deserialized object.
roundtrip_rds <- function(obj) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(obj, path)
  # Force the cached blob out of memory and back in from disk.
  out <- readRDS(path)
  out
}

# -----------------------------------------------------------------------------
# 1. Dense basis + dense loadings
# -----------------------------------------------------------------------------

test_that("LatentNeuroVec with dense basis/loadings round-trips through saveRDS", {
  skip_if_no_neuroim2()

  set.seed(1)
  n_time <- 8L
  k      <- 3L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  basis    <- Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k),  n_vox,  k), sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = basis,
    loadings = loadings,
    space    = space,
    mask     = mask,
    label    = "dense-roundtrip"
  )

  lvec2 <- roundtrip_rds(lvec)

  expect_s4_class(lvec2, "LatentNeuroVec")
  expect_equal(
    as.matrix(fmrilatent:::basis_mat(lvec2@basis)),
    as.matrix(basis),
    tolerance = 1e-12
  )
  expect_equal(
    as.matrix(fmrilatent:::loadings_mat(lvec2@loadings)),
    as.matrix(loadings),
    tolerance = 1e-12
  )
  expect_identical(lvec2@label, "dense-roundtrip")
})

# -----------------------------------------------------------------------------
# 2. BasisHandle (kind = "dct") survives round-trip with empty cache
# -----------------------------------------------------------------------------

test_that("BasisHandle(kind='dct') rematerializes after saveRDS+empty-cache readRDS", {
  skip_if_no_neuroim2()

  set.seed(2)
  n_time <- 12L
  k      <- 4L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  bh <- dct_basis_handle(n_time = n_time, k = k, id = "dct-rt-test")
  expected_basis <- as.matrix(build_dct_basis(n_time = n_time, k = k))

  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k), n_vox, k), sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = loadings,
    space    = space,
    mask     = mask,
    label    = "dct-handle"
  )

  # Pre-warm the cache, then save, then wipe the cache so reload cannot cheat.
  fmrilatent:::basis_mat(lvec@basis)
  expect_true(fmrilatent:::.latent_has_matrix("dct-rt-test", "basis"))

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(lvec, path)

  fmrilatent_registry_clear()
  expect_false(fmrilatent:::.latent_has_matrix("dct-rt-test", "basis"))

  lvec2 <- readRDS(path)
  expect_s4_class(lvec2@basis, "BasisHandle")
  expect_identical(lvec2@basis@id,   "dct-rt-test")
  expect_identical(lvec2@basis@kind, "dct")

  # Materialization must succeed from spec alone with the cache empty.
  reconstructed <- as.matrix(fmrilatent:::basis_mat(lvec2@basis))
  expect_equal(reconstructed, expected_basis, tolerance = 1e-10)

  # And the cache must now hold the rehydrated entry.
  expect_true(fmrilatent:::.latent_has_matrix("dct-rt-test", "basis"))
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# 3. BasisHandle (kind = "bspline") round-trip
# -----------------------------------------------------------------------------

test_that("BasisHandle(kind='bspline') rematerializes after round-trip", {
  skip_if_no_neuroim2()

  set.seed(3)
  n_time <- 16L
  k      <- 5L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  bh <- bspline_basis_handle(n_time = n_time, k = k, id = "bspline-rt-test")
  expected_basis <- as.matrix(build_bspline_basis(n_time = n_time, k = k))

  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k), n_vox, k), sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = loadings,
    space    = space,
    mask     = mask
  )

  fmrilatent_registry_clear()
  lvec2 <- roundtrip_rds(lvec)

  expect_s4_class(lvec2@basis, "BasisHandle")
  expect_identical(lvec2@basis@kind, "bspline")

  reconstructed <- as.matrix(fmrilatent:::basis_mat(lvec2@basis))
  expect_equal(reconstructed, expected_basis, tolerance = 1e-10)
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# 4. BasisHandle (kind = "explicit") round-trip
# -----------------------------------------------------------------------------

test_that("BasisHandle(kind='explicit') round-trips with embedded matrix", {
  skip_if_no_neuroim2()

  set.seed(4)
  n_time <- 6L
  k      <- 2L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  payload <- matrix(rnorm(n_time * k), n_time, k)
  bh <- new("BasisHandle",
            id    = "explicit-basis-rt",
            dim   = as.integer(c(n_time, k)),
            kind  = "explicit",
            spec  = list(matrix = payload),
            label = "explicit basis")

  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k), n_vox, k), sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = loadings,
    space    = space,
    mask     = mask
  )

  fmrilatent_registry_clear()
  lvec2 <- roundtrip_rds(lvec)

  expect_s4_class(lvec2@basis, "BasisHandle")
  reconstructed <- as.matrix(fmrilatent:::basis_mat(lvec2@basis))
  expect_equal(reconstructed, payload, tolerance = 1e-12)
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# 5. LoadingsHandle (kind = "explicit") round-trip
# -----------------------------------------------------------------------------

test_that("LoadingsHandle(kind='explicit') round-trips with embedded matrix", {
  skip_if_no_neuroim2()

  set.seed(5)
  n_time <- 8L
  k      <- 3L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  payload <- matrix(rnorm(n_vox * k), n_vox, k)
  lh <- new("LoadingsHandle",
            id    = "explicit-loadings-rt",
            dim   = as.integer(c(n_vox, k)),
            kind  = "explicit",
            spec  = list(matrix = payload),
            label = "explicit loadings")

  basis <- Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  space <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = basis,
    loadings = lh,
    space    = space,
    mask     = mask
  )

  fmrilatent_registry_clear()
  lvec2 <- roundtrip_rds(lvec)

  expect_s4_class(lvec2@loadings, "LoadingsHandle")
  reconstructed <- as.matrix(fmrilatent:::loadings_mat(lvec2@loadings))
  expect_equal(reconstructed, payload, tolerance = 1e-12)
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# 6. Both basis and loadings as handles round-trip together
# -----------------------------------------------------------------------------

test_that("LatentNeuroVec with both BasisHandle and LoadingsHandle round-trips", {
  skip_if_no_neuroim2()

  set.seed(6)
  n_time <- 10L
  k      <- 3L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  bh <- dct_basis_handle(n_time = n_time, k = k, id = "dual-rt-basis")
  expected_basis <- as.matrix(build_dct_basis(n_time = n_time, k = k))

  loadings_payload <- matrix(rnorm(n_vox * k), n_vox, k)
  lh <- new("LoadingsHandle",
            id    = "dual-rt-loadings",
            dim   = as.integer(c(n_vox, k)),
            kind  = "explicit",
            spec  = list(matrix = loadings_payload),
            label = "dual loadings")

  space <- make_small_space(n_time = n_time)
  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = lh,
    space    = space,
    mask     = mask
  )

  fmrilatent_registry_clear()
  lvec2 <- roundtrip_rds(lvec)

  expect_s4_class(lvec2@basis,    "BasisHandle")
  expect_s4_class(lvec2@loadings, "LoadingsHandle")
  expect_equal(
    as.matrix(fmrilatent:::basis_mat(lvec2@basis)),
    expected_basis,
    tolerance = 1e-10
  )
  expect_equal(
    as.matrix(fmrilatent:::loadings_mat(lvec2@loadings)),
    loadings_payload,
    tolerance = 1e-12
  )
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# 7. Round-trip survives a disabled registry on the receiving side
# -----------------------------------------------------------------------------

test_that("Handle materialization works after readRDS even when registry is disabled", {
  skip_if_no_neuroim2()

  set.seed(7)
  n_time <- 8L
  k      <- 3L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  bh <- dct_basis_handle(n_time = n_time, k = k, id = "disabled-registry-rt")
  expected_basis <- as.matrix(build_dct_basis(n_time = n_time, k = k))

  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k), n_vox, k), sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = loadings,
    space    = space,
    mask     = mask
  )

  path <- tempfile(fileext = ".rds")
  on.exit({
    unlink(path)
    fmrilatent_registry_enable()
    fmrilatent_registry_clear()
  }, add = TRUE)

  fmrilatent_registry_clear()
  saveRDS(lvec, path)

  fmrilatent_registry_disable()
  on.exit(fmrilatent_registry_enable(), add = TRUE)

  lvec2 <- readRDS(path)
  expect_s4_class(lvec2@basis, "BasisHandle")

  reconstructed <- as.matrix(fmrilatent:::basis_mat(lvec2@basis))
  expect_equal(reconstructed, expected_basis, tolerance = 1e-10)

  # Registry stays empty when disabled.
  expect_false(fmrilatent_registry_enabled())
  expect_false(fmrilatent:::.latent_has_matrix("disabled-registry-rt", "basis"))
})

# -----------------------------------------------------------------------------
# 8. Cross-session simulation: drop all R refs, force GC, then reload
# -----------------------------------------------------------------------------

test_that("LatentNeuroVec round-trip survives full reference drop + gc", {
  skip_if_no_neuroim2()

  set.seed(8)
  n_time <- 8L
  k      <- 2L
  mask   <- make_small_mask()
  n_vox  <- sum(as.array(mask))

  bh <- dct_basis_handle(n_time = n_time, k = k, id = "gc-rt-test")
  expected_basis <- as.matrix(build_dct_basis(n_time = n_time, k = k))

  loadings_payload <- matrix(rnorm(n_vox * k), n_vox, k)
  loadings <- Matrix::Matrix(loadings_payload, sparse = FALSE)
  space    <- make_small_space(n_time = n_time)

  lvec <- LatentNeuroVec(
    basis    = bh,
    loadings = loadings,
    space    = space,
    mask     = mask
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  fmrilatent_registry_clear()
  fmrilatent:::basis_mat(lvec@basis)  # populate cache
  saveRDS(lvec, path)

  rm(lvec, bh, loadings)
  fmrilatent_registry_clear()
  gc(verbose = FALSE)

  lvec2 <- readRDS(path)
  expect_s4_class(lvec2, "LatentNeuroVec")
  expect_equal(
    as.matrix(fmrilatent:::basis_mat(lvec2@basis)),
    expected_basis,
    tolerance = 1e-10
  )
  expect_equal(
    as.matrix(fmrilatent:::loadings_mat(lvec2@loadings)),
    loadings_payload,
    tolerance = 1e-12
  )
  fmrilatent_registry_clear()
})
