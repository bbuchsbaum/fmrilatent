library(testthat)

with_clean_loadings_registry <- function(code) {
  old_enabled <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = old_enabled), add = TRUE)
  fmrilatent_registry_enable()
  fmrilatent_registry_clear("loadings")
  on.exit(fmrilatent_registry_clear("loadings"), add = TRUE)
  force(code)
}

expect_constructor_uses_cached_loadings <- function(constructor,
                                                    family,
                                                    kind,
                                                    id,
                                                    label,
                                                    include_data = TRUE) {
  reduction <- list(marker = paste0(family, "-reduction"))
  basis_spec <- list(marker = paste0(family, "-basis"))
  cached <- Matrix::Matrix(matrix(seq_len(6L), nrow = 2L), sparse = FALSE)
  spec_payload <- list(
    family = family,
    reduction = reduction,
    basis_spec = basis_spec,
    k_neighbors = 6L
  )
  if (include_data) {
    spec_payload <- append(spec_payload, list(data = NULL), after = 3L)
  }
  expected <- new(
    "LoadingsHandle",
    id = id,
    dim = as.integer(dim(cached)),
    kind = kind,
    spec = spec_payload,
    label = label
  )
  fmrilatent:::.latent_register_matrix(
    id,
    cached,
    type = "loadings",
    overwrite = TRUE,
    fingerprint = fmrilatent:::.latent_handle_fingerprint(expected)
  )

  handle <- constructor(
    reduction = reduction,
    basis_spec = basis_spec,
    id = id,
    label = label
  )

  expect_s4_class(handle, "LoadingsHandle")
  expect_identical(handle@dim, as.integer(dim(cached)))
  expect_identical(handle@kind, kind)
  expect_equal(as.matrix(fmrilatent:::loadings_mat(handle)), as.matrix(cached))
}

test_that("heat_wavelet_loadings_handle returns cached loadings without lift", {
  with_clean_loadings_registry({
    expect_constructor_uses_cached_loadings(
      heat_wavelet_loadings_handle,
      family = "heat_wavelet",
      kind = "lifted",
      id = "cached-heat-loadings",
      label = "heat-wavelet"
    )
  })
})

test_that("diffusion_wavelet_loadings_handle returns cached loadings without lift", {
  with_clean_loadings_registry({
    expect_constructor_uses_cached_loadings(
      diffusion_wavelet_loadings_handle,
      family = "diffusion_wavelet",
      kind = "lifted",
      id = "cached-diffusion-loadings",
      label = "diffusion-wavelet",
      include_data = FALSE
    )
  })
})

test_that("diffusion_wavelet_loadings_handle id ignores unused data", {
  skip_if_not_installed("rgsp")
  with_clean_loadings_registry({
    mask <- array(TRUE, dim = c(2, 2, 1))
    red <- make_cluster_reduction(mask, seq_len(sum(mask)))
    spec <- basis_diffusion_wavelet(
      target_rank = 1L,
      oversample = 0L,
      threshold = 0,
      max_scales = 1L,
      seed = 7L
    )

    h1 <- diffusion_wavelet_loadings_handle(red, spec, data = matrix(1), k_neighbors = 2L)
    h2 <- diffusion_wavelet_loadings_handle(red, spec, data = matrix(2), k_neighbors = 2L)

    expect_identical(h1@id, h2@id)
    expect_false("data" %in% names(h1@spec))
    expect_false("data" %in% names(h2@spec))
  })
})

test_that("slepian_spatial_loadings_handle returns cached loadings without lift", {
  with_clean_loadings_registry({
    expect_constructor_uses_cached_loadings(
      slepian_spatial_loadings_handle,
      family = "slepian_spatial",
      kind = "slepian_spatial",
      id = "cached-slepian-loadings",
      label = "slepian-spatial"
    )
  })
})
