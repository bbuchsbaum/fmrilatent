#' Benchmark encode/decode round-trips
#'
#' @param mask_dims Integer vector length 3 for spatial dims.
#' @param n_time Integer time points.
#' @param methods Character vector of families: "slepian_space", "hrbf", "wavelet_active", "bspline_hrbf_st".
#' @param iterations Integer iterations per method (default 1 to stay fast).
#' @return A data.frame/tibble with timings and reconstruction error.
#' @export
benchmark_roundtrip <- function(mask_dims = c(16, 16, 8),
                                n_time = 10L,
                                methods = c("slepian_space", "hrbf", "wavelet_active", "bspline_hrbf_st"),
                                iterations = 1L) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("bench package required for benchmark_roundtrip (Suggests).", call. = FALSE)
  }
  mask_arr <- array(TRUE, dim = mask_dims)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(mask_dims))
  n_vox <- sum(mask_arr)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)

  res_list <- list()
  for (m in methods) {
    mark_res <- withCallingHandlers(
      bench::mark(
        iterations = iterations,
        encode_decode = {
          if (m == "slepian_space") {
            lv <- encode(X, spec_space_slepian(k = 3L, k_neighbors = 6L), mask = mask_vol, materialize = "matrix")
            rec <- as.matrix(lv)
          } else if (m == "hrbf") {
            lv <- encode(X, spec_space_hrbf(params = list(
              sigma0 = 2, levels = 0L, radius_factor = 2.5,
              kernel_type = "gaussian", seed = 1L)),
              mask = mask_vol, materialize = "matrix")
            rec <- as.matrix(lv)
          } else if (m == "wavelet_active") {
            lv <- encode(X,
              spec_space_wavelet_active(levels_space = 1L,
                levels_time = 0L, threshold = 0),
              mask = mask_vol, materialize = "matrix")
            rec <- predict(lv)
          } else if (m == "bspline_hrbf_st") {
            lv <- latent_factory("bspline_hrbf_st", x = X, mask = mask_vol,
                                 k_time = min(n_time, 4L),
                                 params = list(sigma0 = 2, levels = 0L,
                                   radius_factor = 2.5,
                                   kernel_type = "gaussian", seed = 2L),
                                 materialize = "matrix")
            rec <- predict(lv)
          } else {
            stop("Unknown method: ", m, call. = FALSE)
          }
          list(rec = rec)
        }
      ),
      warning = function(w) {
        if (grepl("GC in every iteration", conditionMessage(w), ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    rec <- mark_res$result[[1]]$rec
    err <- sqrt(mean((rec - X)^2))
    res_list[[length(res_list) + 1L]] <- data.frame(
      method = m,
      median_ms = as.numeric(mark_res$median) / 1e6,
      itr = iterations,
      rmse = err,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res_list)
}

#' Plot benchmark results
#'
#' @param df Data frame from benchmark_roundtrip.
#' @return ggplot object if ggplot2 available; otherwise prints df.
#' @export
plot_benchmark_roundtrip <- function(df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    print(df)
    return(invisible(df))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = median_ms, color = method)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(y = "Median time (ms)", x = NULL, title = "Encode+decode round-trip benchmark") +
    ggplot2::theme_minimal()
}
