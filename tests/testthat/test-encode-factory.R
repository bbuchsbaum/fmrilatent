test_that("encode with spec_time_slepian returns LatentNeuroVec", {
  mask <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(2, 2, 1)))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1, k = 4)
  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(nrow(X), 4))
})

test_that("latent_factory slepian_time mirrors encode", {
  mask <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(2, 2, 1)))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  lv1 <- latent_factory("slepian_time", x = X, mask = mask_vol, tr = 2, bandwidth = 0.1, k = 3, materialize = "matrix")
  lv2 <- encode(X, spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3), mask = mask_vol, materialize = "matrix")
  expect_equal(as.matrix(lv1), as.matrix(lv2), tolerance = 1e-8)
})

test_that("latent_factory slepian_st builds ImplicitLatent", {
  skip_if_not_installed("RSpectra")
  mask <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(2, 2, 1)))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  lv <- latent_factory("slepian_st", x = X, mask = mask_vol,
                       tr = 2, bandwidth = 0.1, k_time = 3, k_space = 2, k_neighbors = 2)
  expect_true(is_implicit_latent(lv))
  # small sanity: decode dims
  reco <- predict(lv)
  expect_equal(dim(reco), dim(X))
})
