test_that("extract_mask_array accepts list masks and rejects dimensionless masks", {
  mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))

  expect_equal(
    fmrilatent:::.extract_mask_array(list(arr = mask_arr), "list-mask-test"),
    mask_arr
  )
  expect_error(
    fmrilatent:::.extract_mask_array(list(arr = c(TRUE, FALSE)), "bad-list-mask"),
    "mask must have dimensions",
    class = "fmrilatent_error_mask_extract"
  )
})

test_that("extract_mask_array preserves mask conversion parent error", {
  err <- rlang::catch_cnd(
    fmrilatent:::.extract_mask_array(new.env(), "parent-mask-test"),
    classes = "error"
  )

  expect_s3_class(err, "fmrilatent_error_mask_extract")
  expect_s3_class(err$parent, "fmrilatent_error_invalid_mask")
  expect_match(conditionMessage(err$parent), "Underlying error")
})

test_that("sparse dictionary builder skips empty clusters and preserves columns", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  reduction <- new(
    "ClusterReduction",
    mask = mask_vol,
    info = list(),
    map = as.integer(c(1L, 1L, 3L, 3L)),
    cluster_ids = as.integer(c(1L, 2L, 3L))
  )
  seen <- integer()

  out <- fmrilatent:::.build_sparse_dictionary_from_clusters(
    reduction,
    function(vox_idx, cid) {
      seen <<- c(seen, cid)
      Matrix::sparseMatrix(
        i = seq_along(vox_idx),
        j = seq_along(vox_idx),
        x = rep(as.numeric(cid), length(vox_idx)),
        dims = c(length(vox_idx), length(vox_idx))
      )
    }
  )

  expect_identical(seen, c(1L, 3L))
  expect_s4_class(out, "sparseMatrix")
  expect_equal(dim(out), c(4L, 4L))
  expect_equal(as.matrix(out)[1:2, 1:2], diag(1, 2))
  expect_equal(as.matrix(out)[3:4, 3:4], diag(3, 2))
})

test_that("sparse dictionary builder returns zero-column matrix when all clusters are empty", {
  mask_arr <- array(TRUE, dim = c(3, 1, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 1, 1)))
  reduction <- new(
    "ClusterReduction",
    mask = mask_vol,
    info = list(),
    map = as.integer(c(1L, 1L, 1L)),
    cluster_ids = as.integer(c(2L, 3L))
  )
  called <- FALSE

  out <- fmrilatent:::.build_sparse_dictionary_from_clusters(
    reduction,
    function(vox_idx, cid) {
      called <<- TRUE
      Matrix::Diagonal(n = length(vox_idx))
    }
  )

  expect_false(called)
  expect_s4_class(out, "sparseMatrix")
  expect_equal(dim(out), c(3L, 0L))
})

test_that("make_latent_neurovector synthesizes basis from data and loadings", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  loadings <- Matrix::Diagonal(n = sum(mask_arr))
  X <- matrix(seq_len(12), nrow = 3L)

  lvec <- fmrilatent:::.make_latent_neurovector(
    X,
    mask = mask_vol,
    loadings = loadings,
    basis = NULL,
    label = "synth-basis",
    meta = list(family = "test"),
    location = "make-latent-test",
    expect_dense = TRUE
  )

  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(as.matrix(basis(lvec)), as.matrix(X %*% loadings), tolerance = 1e-8)
  expect_equal(as.matrix(lvec), X, tolerance = 1e-8)
  expect_equal(lvec@label, "synth-basis")
  expect_equal(lvec@meta$family, "test")
})
