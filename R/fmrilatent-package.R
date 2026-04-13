#' fmrilatent: Latent Space Representations of fMRI Data
#'
#' `fmrilatent` owns latent representations for neuroimaging data, including
#' explicit factorized objects, shared basis assets, operator-backed latent
#' objects, coefficient recovery, and coefficient-to-map projection.
#'
#' The package deliberately stops at the latent-method layer. It does not own
#' first-level GLM fitting, contrasts, temporal autocorrelation modeling, or
#' statistical inference. Downstream modeling packages should consume
#' `coef_time(x, "analysis")` plus the decoder and projection generics exposed
#' here.
#'
#' @useDynLib fmrilatent, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom neuroim2 NeuroSpace LogicalNeuroVol SparseNeuroVol IndexLookupVol NeuroVecSeq
#' @importFrom neuroim2 drop_dim space spacing origin trans axes indices
#'   lookup matricized_access linear_access series concat mask
#' @importFrom graphics image legend matplot
#' @importFrom stats aggregate predict rnorm
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Declare global variables used in NSE contexts (ggplot2 aes, data.table)
# to avoid R CMD check NOTEs
utils::globalVariables(c(
  "i", "j", "val",                    # plot_basis_gram

  "method", "median_ms",              # plot_benchmark_roundtrip

  "time", "value", "component",       # plot_slepian_temporal
  "x", "y",                           # plot_spatial_atom

  "forward_lift_matrix_rcpp",          # haar_wavelet Rcpp (conditional)
  "inverse_lift_matrix_rcpp"           # haar_wavelet Rcpp (conditional)
))
