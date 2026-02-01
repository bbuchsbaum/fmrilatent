#' @useDynLib fmrilatent, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom neuroim2 NeuroSpace LogicalNeuroVol SparseNeuroVol IndexLookupVol NeuroVecSeq
#' @importFrom neuroim2 drop_dim space spacing origin trans axes indices lookup matricized_access linear_access series concat mask
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
  "x", "y"                            # plot_spatial_atom
))
