# B-spline temporal basis and handle

#' Build a B-spline temporal basis matrix
#'
#' @param n_time Integer, number of time points.
#' @param k      Integer, number of spline basis functions (df).
#' @param degree Degree of the spline (default 3 = cubic).
#' @param knots  Optional numeric vector of interior knots in [0, 1].
#' @param boundary_knots Optional numeric vector length 2 giving boundary knots in [0, 1].
#' @param include_intercept Logical; passed to \code{splines::bs()} (default FALSE).
#'
#' @return Dense Matrix (n_time x k).
#' @keywords internal
build_bspline_basis <- function(n_time,
                                k,
                                degree = 3L,
                                knots = NULL,
                                boundary_knots = NULL,
                                include_intercept = FALSE,
                                orthonormalize = TRUE) {
  n_time <- as.integer(n_time)
  k <- as.integer(k)
  degree <- as.integer(degree)
  if (k > n_time) stop("k (", k, ") cannot exceed n_time (", n_time, ").")
  if (degree < 1) stop("degree must be >= 1")
  if (k <= degree && !include_intercept) {
    stop("k (", k, ") must be greater than degree (", degree,
      ") for a valid B-spline basis.")
  }

  # Ensure boundary knots defined
  if (is.null(boundary_knots)) boundary_knots <- c(0, 1)
  boundary_knots <- as.numeric(boundary_knots)

  # If user didn't supply interior knots, create evenly spaced ones when needed
  if (is.null(knots)) {
    n_int <- max(0L, k - degree - (if (include_intercept) 0L else 1L))
    if (n_int > 0L) {
      knots <- seq(boundary_knots[1], boundary_knots[2], length.out = n_int + 2L)[-c(1L, n_int + 2L)]
    }
  }

  t_scaled <- seq_len(n_time) / n_time
  df_use <- if (include_intercept) k else k + 1L
  bs_mat <- splines::bs(
    t_scaled,
    df = df_use,
    degree = degree,
    knots = knots,
    Boundary.knots = boundary_knots,
    intercept = include_intercept
  )
  # Align to requested column count k
  if (ncol(bs_mat) > k) {
    bs_mat <- bs_mat[, seq_len(k), drop = FALSE]
  } else if (ncol(bs_mat) < k) {
    extra <- matrix(0, nrow = nrow(bs_mat), ncol = k - ncol(bs_mat))
    bs_mat <- cbind(bs_mat, extra)
  }

  bs_mat <- as.matrix(bs_mat)
  if (orthonormalize) {
    q <- qr(bs_mat)
    bs_mat <- qr.Q(q)
  }

  Matrix::Matrix(bs_mat, sparse = FALSE)
}

#' Create a BasisHandle for a B-spline temporal basis
#'
#' @param n_time Integer, number of time points.
#' @param k      Integer, number of spline basis functions (df).
#' @param degree Degree of the spline (default 3).
#' @param knots  Optional interior knots (scaled 0-1).
#' @param boundary_knots Optional boundary knots (scaled 0-1).
#' @param include_intercept Logical; include intercept column (default FALSE).
#' @param id     Optional registry key (generated if NULL).
#' @param label  Optional human-readable label.
#'
#' @return A \code{BasisHandle}.
#' @export
bspline_basis_handle <- function(n_time,
                                 k,
                                 degree = 3L,
                                 knots = NULL,
                                 boundary_knots = NULL,
                                 include_intercept = FALSE,
                                 id = NULL,
                                 label = NULL) {
  n_time <- as.integer(n_time)
  k <- as.integer(k)
  degree <- as.integer(degree)

  if (is.null(id)) {
    id <- sprintf("bspline-%d-%d-deg%d-int%s", n_time, k, degree, include_intercept)
  }
  if (is.null(label)) {
    label <- sprintf("B-spline(%d,%d,deg=%d)", n_time, k, degree)
  }

  new("BasisHandle",
      id    = id,
      dim   = as.integer(c(n_time, k)),
      kind  = "bspline",
      spec  = list(
        n_time = n_time,
        k = k,
        degree = degree,
        knots = knots,
        boundary_knots = boundary_knots,
        include_intercept = include_intercept
      ),
      label = label)
}
