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
#' @section Dispatch model:
#' The registry is for \strong{discovery and introspection only}. It does
#' \emph{not} participate in actual encoding. \code{\link{encode}} routes work
#' through the S3 generic \code{encode_spec}, which dispatches on the class of
#' the spec object. That means an external package contributing a new encoder
#' family must do \strong{both}:
#' \enumerate{
#'   \item Define a spec constructor whose return value carries a distinctive
#'         class (e.g. \code{class(spec) <- c("spec_my_encoder", "list")}).
#'   \item Define an S3 method \code{encode_spec.spec_my_encoder(x, spec, ...)}
#'         that builds the actual \code{LatentNeuroVec} (or other latent type).
#'   \item Optionally call \code{register_encoder()} so that
#'         \code{list_encoders()} surfaces the family for users.
#' }
#' Calling \code{register_encoder()} alone is \strong{not enough}: without an
#' \code{encode_spec.spec_*} method, \code{encode()} will fall through to
#' \code{encode_spec.default} and raise an error.
#'
#' Transport-backed AWPT is an intentional parallel API, not a registry-dispatch
#' family. \code{basis_awpt_wavelet()} describes a shared template basis, while
#' subject fitting also needs a template/basis asset plus subject field or
#' observation operators. Use \code{\link{encode_awpt}} or
#' \code{\link{encode_operator}} for those fits; do not expect
#' \code{register_encoder()} or \code{encode()} to route AWPT subject encoding.
#'
#' @details
#' External packages typically register their encoders in their \code{.onLoad}
#' hook so the registry is populated when the package is loaded:
#' \preformatted{
#' .onLoad <- function(libname, pkgname) {
#'   if (requireNamespace("fmrilatent", quietly = TRUE)) {
#'     fmrilatent::register_encoder(
#'       "my_encoder", spec_my_encoder,
#'       description = "My custom encoder",
#'       package = pkgname
#'     )
#'   }
#' }
#' }
#' The S3 method itself is exported via the package's NAMESPACE in the usual
#' way (\code{S3method(encode_spec, spec_my_encoder)}).
#'
#' Registering the same family twice issues a warning and overwrites the
#' previous entry.
#'
#' @return Invisibly, \code{TRUE}.
#' @seealso \code{\link{encode}}, \code{\link{list_encoders}}, \code{\link{get_encoder}}
#' @export
#' @examples
#' register_encoder("test_enc", identity, "A test encoder", "mypkg")
#' list_encoders()
#' get_encoder("test_enc")
register_encoder <- function(family, spec_fn, description = "", package = "") {
  if (!is.character(family) || length(family) != 1L || is.na(family) || !nzchar(family)) {
    .encoder_cli_abort(
      "register_encoder() requires a non-empty character `family`.",
      class = "fmrilatent_error_type",
      call = rlang::caller_env()
    )
  }
  if (!is.function(spec_fn)) {
    .encoder_cli_abort(
      "register_encoder() requires `spec_fn` to be a function.",
      class = "fmrilatent_error_type",
      call = rlang::caller_env()
    )
  }
  if (exists(family, envir = .encoder_registry_env, inherits = FALSE)) {
    .encoder_cli_warn(
      paste0("Encoder '", family, "' is already registered; overwriting."),
      class = "fmrilatent_warning_encoder", call = rlang::caller_env())
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
#' enc <- get_encoder("time_dct")
#' spec <- enc$spec_fn(k = 5)
get_encoder <- function(family) {
  if (!is.character(family) || length(family) != 1L || is.na(family) || !nzchar(family)) {
    .encoder_cli_abort(
      "get_encoder() requires a non-empty character `family`.",
      class = "fmrilatent_error_type",
      call = rlang::caller_env()
    )
  }
  if (!exists(family, envir = .encoder_registry_env, inherits = FALSE)) {
    available <- ls(.encoder_registry_env)
    msg <- paste0("No encoder registered for family '", family, "'.")
    if (length(available) > 0L) {
      msg <- paste0(msg, " Available encoders: ",
                    paste(available, collapse = ", "))
    }
    .encoder_cli_abort(msg,
                       class = "fmrilatent_error_value", call = rlang::caller_env())
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
  if (length(dims) != 3L || anyNA(dims) || any(dims < 1L)) {
    .encoder_cli_abort(
      "fmrilatent_test_data() requires `dims` to be a length-3 positive integer vector.",
      class = "fmrilatent_error_dim",
      call = rlang::caller_env()
    )
  }
  if (length(n_time) != 1L || is.na(n_time) || n_time < 1L) {
    .encoder_cli_abort(
      "fmrilatent_test_data() requires `n_time` to be a positive integer.",
      class = "fmrilatent_error_invalid_count",
      call = rlang::caller_env()
    )
  }
  n_vox <- prod(dims)
  X <- matrix(stats::rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  mask <- array(TRUE, dim = dims)
  list(X = X, mask = mask, dims = dims, n_time = n_time)
}
