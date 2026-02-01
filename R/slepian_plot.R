# Lightweight plotting helpers for Slepian bases

#' Plot temporal Slepians (DPSS)
#'
#' @param basis Matrix (time x k) or BasisHandle (kind = "slepian_temporal").
#' @param max_components Maximum components to display (default 6).
#' @return A ggplot object (if ggplot2 available); otherwise invisibly plots with base graphics.
#' @export
plot_slepian_temporal <- function(basis, max_components = 6L) {
  if (inherits(basis, "BasisHandle")) {
    basis <- basis_mat(basis)
  }
  if (!is.matrix(basis) && !inherits(basis, "Matrix")) {
    stop("basis must be a matrix or BasisHandle (slepian_temporal).", call. = FALSE)
  }
  k <- min(ncol(basis), max_components)
  B <- as.matrix(basis[, seq_len(k), drop = FALSE])

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    matplot(B, type = "l", lty = 1, xlab = "Time", ylab = "Amplitude",
            main = "Temporal Slepians")
    legend("topright", legend = paste0("k", seq_len(k)), col = seq_len(k), lty = 1, cex = 0.8)
    return(invisible(NULL))
  }

  df <- data.frame(
    time = rep(seq_len(nrow(B)), times = k),
    component = factor(rep(seq_len(k), each = nrow(B))),
    value = as.vector(B)
  )
  ggplot2::ggplot(df, ggplot2::aes(time, value, color = component)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Time", y = "Amplitude", color = "Component",
                  title = "Temporal Slepians (DPSS)") +
    ggplot2::theme_minimal()
}

#' Plot Gram matrix of a basis (orthogonality check)
#'
#' @param basis Matrix or BasisHandle.
#' @return ggplot heatmap (if ggplot2 available); otherwise shows base image.
#' @export
plot_basis_gram <- function(basis) {
  if (inherits(basis, "BasisHandle")) {
    basis <- basis_mat(basis)
  }
  if (!is.matrix(basis) && !inherits(basis, "Matrix")) {
    stop("basis must be a matrix or BasisHandle.", call. = FALSE)
  }
  G <- crossprod(as.matrix(basis))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    image(G, main = "Basis Gram matrix", xlab = "Component", ylab = "Component")
    return(invisible(NULL))
  }
  k <- ncol(G)
  df <- data.frame(
    i = rep(seq_len(k), each = k),
    j = rep(seq_len(k), times = k),
    val = as.vector(G)
  )
  ggplot2::ggplot(df, ggplot2::aes(j, i, fill = val)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = "Component", y = "Component", fill = "Gij",
                  title = "Basis Gram matrix") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal()
}
