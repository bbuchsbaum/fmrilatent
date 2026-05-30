# Regression tests for classed conditions in non-encoder infrastructure.
# These assert that bare stop()/warning() calls converted to classed
# conditions (mote bd-01KSWB3Y55C6XAM44AHHXQB51M) carry the expected
# fmrilatent_error_* / fmrilatent_warning_* class AND preserve the
# original human-readable message text.

test_that("LatentNeuroVec rejects a non-NeuroSpace space (type)", {
  expect_error(
    LatentNeuroVec(basis = matrix(0, 2, 1), loadings = matrix(0, 1, 1),
                   space = "not a space", mask = TRUE),
    class = "fmrilatent_error_type"
  )
  expect_error(
    LatentNeuroVec(basis = matrix(0, 2, 1), loadings = matrix(0, 1, 1),
                   space = "not a space", mask = TRUE),
    "'space' must be a NeuroSpace object"
  )
})

test_that("BasisHandle with unknown kind aborts (value)", {
  h <- BasisHandle(kind = "definitely-not-a-kind", spec = list(),
                   dim = c(2L, 1L))
  expect_error(
    fmrilatent:::basis_mat(h),
    class = "fmrilatent_error_value"
  )
  expect_error(
    fmrilatent:::basis_mat(h),
    "Unknown BasisHandle kind"
  )
})

test_that("LoadingsHandle with unknown kind aborts (value)", {
  h <- LoadingsHandle(kind = "definitely-not-a-kind", spec = list(),
                      dim = c(1L, 1L))
  expect_error(
    fmrilatent:::loadings_mat(h),
    class = "fmrilatent_error_value"
  )
  expect_error(
    fmrilatent:::loadings_mat(h),
    "Unknown LoadingsHandle kind"
  )
})

test_that("BilatLatentNeuroSurfaceVector validates its inputs (type)", {
  expect_error(
    BilatLatentNeuroSurfaceVector(1, 2),
    class = "fmrilatent_error_type"
  )
  expect_error(
    BilatLatentNeuroSurfaceVector(1, 2),
    "'left' must be a LatentNeuroSurfaceVector"
  )
})

test_that("LatentNeuroSurfaceVector rejects a bad geometry", {
  skip_if_not_installed("neurosurf")
  expect_error(
    LatentNeuroSurfaceVector(geometry = "nope",
                             basis = matrix(0, 2, 1),
                             loadings = matrix(0, 1, 1)),
    class = "fmrilatent_error_type"
  )
})

test_that("SurfaceBasisTemplate requires neurosurf when absent (missing dependency)", {
  skip_if(requireNamespace("neurosurf", quietly = TRUE),
          "neurosurf is installed; cannot exercise the missing-dependency path")
  expect_error(
    SurfaceBasisTemplate(geometry = "nope", loadings = matrix(0, 1, 1)),
    class = "fmrilatent_error_missing_dependency"
  )
  expect_error(
    SurfaceBasisTemplate(geometry = "nope", loadings = matrix(0, 1, 1)),
    "requires the 'neurosurf' package"
  )
})

test_that("BlockLatentNeuroVector rejects an empty block list (value)", {
  expect_error(
    BlockLatentNeuroVector(list()),
    class = "fmrilatent_error_value"
  )
  expect_error(
    BlockLatentNeuroVector(list()),
    "requires a non-empty list of latent blocks"
  )
})

test_that("ClusterReduction validates cluster_map length (dim)", {
  expect_error(
    ClusterReduction(cluster_map = c(1L, 2L), n_voxels = 5L),
    class = "fmrilatent_error_dim"
  )
  expect_error(
    ClusterReduction(cluster_map = c(1L, 2L), n_voxels = 5L),
    "ClusterReduction requires a cluster_map of length n_voxels"
  )
})

test_that("graph_from_reduction rejects a non-GraphReduction (type)", {
  expect_error(
    graph_from_reduction(list(a = 1)),
    class = "fmrilatent_error_type"
  )
  expect_error(
    graph_from_reduction(list(a = 1)),
    "'reduction' must be a GraphReduction object"
  )
})

test_that("register_encoder validates its family argument (type)", {
  expect_error(
    register_encoder("", identity),
    class = "fmrilatent_error_type"
  )
  expect_error(
    register_encoder("", identity),
    "requires a non-empty character `family`"
  )
})

test_that("get_encoder reports an unknown family (value)", {
  expect_error(
    get_encoder("no-such-encoder-family-xyz"),
    class = "fmrilatent_error_value"
  )
  expect_error(
    get_encoder("no-such-encoder-family-xyz"),
    "No encoder registered for family"
  )
})

test_that("encode_template raw_metric squareness check is classed (dim)", {
  # .rotate_quadratic_to_analysis lives in encode_template.R; the square-matrix
  # guard on raw_metric fires before quad / analysis_transform are touched.
  expect_error(
    fmrilatent:::.rotate_quadratic_to_analysis(NULL, matrix(0, 2, 3), NULL),
    class = "fmrilatent_error_dim"
  )
  expect_error(
    fmrilatent:::.rotate_quadratic_to_analysis(NULL, matrix(0, 2, 3), NULL),
    "raw_metric must be square"
  )
})
