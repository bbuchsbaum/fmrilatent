library(testthat)

# =============================================================================
# C1: opt-in core_mode = "explicit" for spec_st (DCT x HRBF).
# Default ("auto") must remain ImplicitLatent; "explicit" must return a
# LatentNeuroVec whose reconstruction matches the ImplicitLatent decoder.
# =============================================================================

make_st_inputs <- function() {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  set.seed(11)
  n_time <- 6L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  list(X = X, mask = mask_vol, n_time = n_time)
}

st_spec <- function(core_mode = c("auto", "explicit")) {
  spec_st(
    time = spec_time_dct(k = 4L),
    space = spec_space_hrbf(params = list(
      sigma0 = 2, levels = 1L, radius_factor = 2.5,
      kernel_type = "gaussian", seed = 42L
    )),
    core_mode = match.arg(core_mode)
  )
}

test_that("spec_st default (auto) still returns an ImplicitLatent", {
  io <- make_st_inputs()
  lv <- encode(io$X, st_spec("auto"), mask = io$mask, materialize = "matrix")
  expect_true(inherits(lv, "ImplicitLatent"))
  expect_false(is(lv, "LatentNeuroVec"))
})

test_that("spec_st core_mode='explicit' returns a LatentNeuroVec", {
  io <- make_st_inputs()
  lv <- encode(io$X, st_spec("explicit"), mask = io$mask, materialize = "matrix")
  expect_true(is(lv, "LatentNeuroVec"))
  expect_false(inherits(lv, "ImplicitLatent"))
})

test_that("explicit spec_st reconstruction matches the ImplicitLatent decoder", {
  io <- make_st_inputs()
  imp <- encode(io$X, st_spec("auto"), mask = io$mask, materialize = "matrix")
  exp_lv <- encode(io$X, st_spec("explicit"), mask = io$mask, materialize = "matrix")

  rec_imp <- predict(imp)            # time x voxel
  rec_exp <- as.matrix(exp_lv)       # time x voxel

  expect_equal(dim(rec_exp), dim(rec_imp))
  expect_lt(max(abs(rec_exp - rec_imp)), 1e-10)
})
