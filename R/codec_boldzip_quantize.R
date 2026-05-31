# BOLDZip-SR codec: reliability-shaped quantization and the unified noise-scale policy
#
# Split out of R/codec_boldzip.R (see that file's header for the full map).
# Functions in this file:
#   .boldzip_quantize_values
#   boldzip_quantization

# ---------------------------------------------------------------------------
# Unified noise-scale policy for reliability-shaped quantization
#
# BOLDZip-SR quantizes a coefficient vector with a step size that is shaped by
# two quantities: a global `noise_scale` (the typical magnitude of the values
# being quantized, i.e. the amount of "noise" a quantization step may safely
# absorb) and a per-element `reliability` in [0, 1]. The step formula is
#
#   step = base_step * noise_scale / sqrt(max(reliability, 0) + epsilon)
#
# so that more reliable coefficients get a finer step (smaller denominator -> 
# wait: larger reliability -> larger denominator -> smaller step) and less
# reliable coefficients get a coarser step. `noise_scale` sets the overall
# magnitude of the step before reliability shaping.
#
# `.boldzip_noise_scale()` is the single canonical definition of how that
# global scale is derived from a set of values: the sample standard deviation
# of the supplied values when at least two values are available, and zero for
# degenerate length-0/1 inputs. Every call site (texture fit residuals, event
# amplitudes, carrier coefficients) and the default inside
# `.boldzip_quantize_values()` route through this helper so the policy lives in
# exactly one place. Non-finite or non-positive scales are treated downstream
# in `.boldzip_quantize_values()` (replaced with 1).
# ---------------------------------------------------------------------------
.boldzip_noise_scale <- function(values) {
  values <- as.numeric(values)
  if (length(values) < 2L) {
    return(0)
  }
  scale <- stats::sd(values)
  if (!is.finite(scale)) 0 else scale
}


.boldzip_quantize_values <- function(values, reliability, quantization,
                                     noise_scale = NULL) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    return(values)
  }
  base_step <- quantization$base_step %||% 0
  if (!is.finite(base_step) || base_step <= 0) {
    return(values)
  }
  reliability <- rep(as.numeric(reliability), length.out = length(values))
  noise_scale <- rep(as.numeric(noise_scale %||% .boldzip_noise_scale(values)), length.out = length(values))
  noise_scale[!is.finite(noise_scale) | noise_scale <= 0] <- 1
  step <- base_step * noise_scale / sqrt(pmax(reliability, 0) + quantization$epsilon)
  step <- pmax(step, .Machine$double.eps)
  round(values / step) * step
}

#' Construct reliability-aware quantization settings for BOLDZip-SR
#'
#' @param base_step Base quantization step. Set to 0 to disable quantization.
#' @param epsilon Small positive value used in reliability-shaped step sizes.
#' @return A `BoldZipSRQuantization` settings object.
#' @export
boldzip_quantization <- function(base_step = 0, epsilon = 1e-6) {
  base_step <- .boldzip_check_scalar_number(base_step, "base_step", min = 0)
  epsilon <- .boldzip_check_scalar_number(epsilon, "epsilon", min = .Machine$double.eps)
  structure(
    list(base_step = base_step, epsilon = epsilon),
    class = "BoldZipSRQuantization"
  )
}
