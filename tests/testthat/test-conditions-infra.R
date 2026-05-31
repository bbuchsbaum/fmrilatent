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
  h <- dct_basis_handle(n_time = 2L, k = 1L)
  h@kind <- "definitely-not-a-kind"
  expect_error(
    fmrilatent:::basis_mat(h),
    class = "fmrilatent_error_value"
  )
  expect_error(
    fmrilatent:::basis_mat(h),
    "Unknown basis handle kind"
  )
})

test_that("LoadingsHandle with unknown kind aborts (value)", {
  h <- methods::new(
    "LoadingsHandle",
    id = "bad-loadings-kind",
    dim = as.integer(c(1L, 1L)),
    kind = "explicit",
    spec = list(matrix = matrix(1, nrow = 1L, ncol = 1L)),
    label = ""
  )
  h@kind <- "definitely-not-a-kind"
  expect_error(
    fmrilatent:::loadings_mat(h),
    class = "fmrilatent_error_value"
  )
  expect_error(
    fmrilatent:::loadings_mat(h),
    "Unknown loadings handle kind"
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
  skip_if_no_neurosurf_surface_examples()
  expect_error(
    LatentNeuroSurfaceVector(geometry = "nope",
                             basis = matrix(0, 2, 1),
                             loadings = matrix(0, 1, 1)),
    class = "fmrilatent_error_type"
  )
})

test_that("SurfaceBasisTemplate requires neurosurf when absent (missing dependency)", {
  skip_if_no_neurosurf_surface_examples()
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
    make_cluster_reduction(array(TRUE, dim = c(5L, 1L, 1L)), c(1L, 2L)),
    class = "fmrilatent_error_dim"
  )
  expect_error(
    make_cluster_reduction(array(TRUE, dim = c(5L, 1L, 1L)), c(1L, 2L)),
    "map must have length"
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
  expect_error(
    register_encoder("bad_spec_fn", list()),
    class = "fmrilatent_error_type"
  )
})

test_that("get_encoder reports an unknown family (value)", {
  expect_error(
    get_encoder(character()),
    class = "fmrilatent_error_type"
  )
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
  expect_error(
    fmrilatent:::.analysis_transform_from_metric(matrix(0, 2, 3)),
    class = "fmrilatent_error_dim"
  )
  expect_error(
    fmrilatent:::.analysis_transform_from_metric(matrix(0, 2, 3)),
    "raw_metric must be square"
  )
})
