# Encoder registry for discovering and registering spec constructors

# Private environment backing the encoder registry
.encoder_registry_env <- new.env(parent = emptyenv())

#' Register an encoder family
#'
#' Registers a spec constructor so it can be discovered via \code{\link{list_encoders}}
#' and retrieved via \code{\link{get_encoder}}.
#'
#' @param family Character string identifying the encoder family (e.g. "time_slepian").
#' @param spec_fn The spec constructor function (e.g. \code{spec_time_slepian}).
#' @param description Character string describing the encoder.
#' @param package Character string naming the package that provides this encoder.
#'
#' @details
#' External packages can register their encoders in their \code{.onLoad} hook:
#' \preformatted{
#' .onLoad <- function(libname, pkgname) {
#'   if (requireNamespace("fmrilatent", quietly = TRUE)) {
#'     fmrilatent::register_encoder("my_encoder", spec_my_encoder,
#'                                   description = "My custom encoder",
#'                                   package = pkgname)
#'   }
#' }
#' }
#'
#' Registering the same family twice issues a warning and overwrites the previous entry.
#'
#' @return Invisibly, \code{TRUE}.
#' @export
#' @examples
#' register_encoder("test_enc", identity, "A test encoder", "mypkg")
#' list_encoders()
#' get_encoder("test_enc")
register_encoder <- function(family, spec_fn, description = "", package = "") {
  stopifnot(is.character(family), length(family) == 1L, nzchar(family))
  stopifnot(is.function(spec_fn))
  if (exists(family, envir = .encoder_registry_env, inherits = FALSE)) {
    warning("Encoder '", family, "' is already registered; overwriting.",
            call. = FALSE)
  }
  assign(family, list(spec_fn = spec_fn, description = description, package = package),
         envir = .encoder_registry_env)
  invisible(TRUE)
}

#' List registered encoders
#'
#' Returns a data.frame describing all registered encoder families,
#' including both built-in encoders and those added by external packages.
#'
#' @return A data.frame with columns: \code{family}, \code{description}, \code{package}.
#' @export
#' @examples
#' list_encoders()
list_encoders <- function() {
  families <- ls(.encoder_registry_env)
  if (length(families) == 0L) {
    return(data.frame(family = character(0), description = character(0),
                      package = character(0), stringsAsFactors = FALSE))
  }
  entries <- lapply(families, function(f) {
    e <- get(f, envir = .encoder_registry_env, inherits = FALSE)
    data.frame(family = f, description = e$description, package = e$package,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, entries)
}

#' Get a registered encoder
#'
#' Retrieves the registration entry for a given encoder family.
#'
#' @param family Character string identifying the encoder family.
#' @return A list with elements \code{spec_fn}, \code{description}, and \code{package}.
#' @export
#' @examples
#' \dontrun{
#' enc <- get_encoder("time_dct")
#' spec <- enc$spec_fn(k = 5)
#' }
get_encoder <- function(family) {
  stopifnot(is.character(family), length(family) == 1L)
  if (!exists(family, envir = .encoder_registry_env, inherits = FALSE)) {
    available <- ls(.encoder_registry_env)
    msg <- paste0("Encoder '", family, "' is not registered.")
    if (length(available) > 0L) {
      msg <- paste0(msg, " Available encoders: ",
                    paste(available, collapse = ", "))
    }
    stop(msg, call. = FALSE)
  }
  get(family, envir = .encoder_registry_env, inherits = FALSE)
}

#' Generate test data for encoder development
#'
#' Creates a small synthetic dataset suitable for testing \code{encode_spec} methods.
#' Useful for extension packages that implement custom encoders.
#'
#' @param dims Integer vector of spatial dimensions (default \code{c(3, 3, 2)}).
#' @param n_time Number of time points (default 8).
#' @return A list with elements:
#'   \item{X}{Numeric matrix (n_time x n_voxels) of random data.}
#'   \item{mask}{Logical array of dimensions \code{dims}, all \code{TRUE}.}
#'   \item{dims}{The spatial dimensions used.}
#'   \item{n_time}{The number of time points used.}
#' @export
#' @importFrom stats rnorm
#' @examples
#' td <- fmrilatent_test_data()
#' dim(td$X)        # 8 x 18
#' dim(td$mask)     # 3 x 3 x 2
fmrilatent_test_data <- function(dims = c(3L, 3L, 2L), n_time = 8L) {
  dims <- as.integer(dims)
  n_time <- as.integer(n_time)
  n_vox <- prod(dims)
  X <- matrix(stats::rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  mask <- array(TRUE, dim = dims)
  list(X = X, mask = mask, dims = dims, n_time = n_time)
}
