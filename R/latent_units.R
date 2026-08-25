.latent_units_contract <- "fmrilatent.latent-units.v1"
.latent_units_semantic_fields <- c(
  "response_scaling",
  "coefficient_units",
  "loading_normalization",
  "loading_metric",
  "analysis_coordinate_metric",
  "sign_convention"
)

.latent_units_scalar_string <- function(x, what, allow_null = FALSE) {
  if (isTRUE(allow_null) && is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .encoder_cli_abort(
      paste0(what, " must be one non-empty character string."),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  x
}

.latent_units_hash <- function(payload) {
  digest::digest(payload, algo = "sha256", serialize = TRUE)
}

.latent_units_payload <- function(x) {
  c(
    list(
      contract = x$contract,
      status = x$status
    ),
    x[.latent_units_semantic_fields],
    list(
      notes = x$notes,
      reason = x$reason,
      child_compatibility_ids = x$child_compatibility_ids
    )
  )
}

.new_latent_units <- function(status,
                              response_scaling = NA_character_,
                              coefficient_units = NA_character_,
                              loading_normalization = NA_character_,
                              loading_metric = NA_character_,
                              analysis_coordinate_metric = NA_character_,
                              sign_convention = NA_character_,
                              notes = NULL,
                              reason = NULL,
                              child_compatibility_ids = NULL) {
  semantic <- list(
    response_scaling = response_scaling,
    coefficient_units = coefficient_units,
    loading_normalization = loading_normalization,
    loading_metric = loading_metric,
    analysis_coordinate_metric = analysis_coordinate_metric,
    sign_convention = sign_convention
  )
  compatibility_id <- switch(
    status,
    declared = paste0(
      "latent-units:",
      .latent_units_hash(c(list(contract = .latent_units_contract), semantic))
    ),
    undeclared = "latent-units:undeclared",
    incompatible = paste0(
      "latent-units:incompatible:",
      .latent_units_hash(child_compatibility_ids)
    ),
    .encoder_cli_abort(
      "Unknown latent-units status.",
      class = "fmrilatent_error_invalid_argument"
    )
  )
  value <- c(
    list(
      contract = .latent_units_contract,
      status = status,
      compatibility_id = compatibility_id
    ),
    semantic,
    list(
      notes = notes,
      reason = reason,
      child_compatibility_ids = child_compatibility_ids
    )
  )
  value$integrity_id <- paste0(
    "latent-units-record:",
    .latent_units_hash(.latent_units_payload(value))
  )
  structure(value, class = c("fmrilatent_units", "list"))
}

#' Declare a Latent Units Contract
#'
#' Constructs the unit and normalization record captured by [encode()] and the
#' transport encoders. Every semantic field is required: the constructor does
#' not guess response scaling or basis conventions from numeric data.
#'
#' @param response_scaling Scaling applied to the response before encoding,
#'   such as `"raw_signal"`, `"percent_signal_change"`, or another precise
#'   caller-defined convention.
#' @param coefficient_units Physical or statistical units of analysis-coordinate
#'   coefficients.
#' @param loading_normalization Normalization convention for decoder/loading
#'   columns.
#' @param loading_metric Inner product or weighting under which loading
#'   normalization is defined.
#' @param analysis_coordinate_metric Metric convention used by analysis
#'   coordinates.
#' @param sign_convention Rule fixing coefficient/loading signs, or an explicit
#'   declaration that signs are arbitrary.
#' @param notes Optional human-readable note. Notes are integrity checked but do
#'   not change the compatibility ID.
#'
#' @return An immutable declared `fmrilatent_units` record.
#' @export
latent_units_record <- function(response_scaling,
                                coefficient_units,
                                loading_normalization,
                                loading_metric,
                                analysis_coordinate_metric,
                                sign_convention,
                                notes = NULL) {
  fields <- list(
    response_scaling = response_scaling,
    coefficient_units = coefficient_units,
    loading_normalization = loading_normalization,
    loading_metric = loading_metric,
    analysis_coordinate_metric = analysis_coordinate_metric,
    sign_convention = sign_convention
  )
  fields <- lapply(names(fields), function(name) {
    .latent_units_scalar_string(fields[[name]], name)
  })
  names(fields) <- .latent_units_semantic_fields
  notes <- .latent_units_scalar_string(notes, "notes", allow_null = TRUE)
  do.call(.new_latent_units, c(list(status = "declared"), fields, list(notes = notes)))
}

.undeclared_latent_units <- function(reason = "no units record was retained") {
  .new_latent_units(status = "undeclared", reason = reason)
}

.validate_latent_units <- function(x, what = "latent units") {
  if (!inherits(x, "fmrilatent_units") || !is.list(x)) {
    .encoder_cli_abort(
      paste0(what, " must be created by latent_units_record()."),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  required <- c(
    "contract", "status", "compatibility_id",
    .latent_units_semantic_fields,
    "notes", "reason", "child_compatibility_ids", "integrity_id"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) || !identical(x$contract, .latent_units_contract)) {
    .encoder_cli_abort(
      paste0(what, " is malformed or has an unsupported contract."),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  expected <- .new_latent_units(
    status = x$status,
    response_scaling = x$response_scaling,
    coefficient_units = x$coefficient_units,
    loading_normalization = x$loading_normalization,
    loading_metric = x$loading_metric,
    analysis_coordinate_metric = x$analysis_coordinate_metric,
    sign_convention = x$sign_convention,
    notes = x$notes,
    reason = x$reason,
    child_compatibility_ids = x$child_compatibility_ids
  )
  if (!identical(x$compatibility_id, expected$compatibility_id) ||
      !identical(x$integrity_id, expected$integrity_id)) {
    .encoder_cli_abort(
      paste0(what, " failed its integrity check."),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  x
}

.stored_latent_units <- function(x) {
  meta <- if (methods::is(x, "ExplicitLatent")) {
    x@meta
  } else {
    x$meta %||% list()
  }
  value <- meta[["latent_units"]] %||% NULL
  if (is.null(value)) {
    return(NULL)
  }
  .validate_latent_units(value, "stored latent units")
}

.combine_child_latent_units <- function(children) {
  records <- lapply(children, latent_units)
  statuses <- vapply(records, `[[`, character(1), "status")
  ids <- stats::setNames(
    vapply(records, `[[`, character(1), "compatibility_id"),
    names(children)
  )
  if (all(statuses == "declared") && length(unique(ids)) == 1L) {
    return(records[[1L]])
  }
  if (any(statuses != "declared")) {
    return(.new_latent_units(
      status = "undeclared",
      reason = "one or more child representations have undeclared units",
      child_compatibility_ids = ids
    ))
  }
  .new_latent_units(
    status = "incompatible",
    reason = "child representations declare incompatible units",
    child_compatibility_ids = ids
  )
}

.with_latent_units <- function(x, units) {
  if (is.null(units)) {
    return(x)
  }
  units <- .validate_latent_units(units, "units")
  if (!identical(units$status, "declared")) {
    .encoder_cli_abort(
      "Encode-time units must have status 'declared'.",
      class = "fmrilatent_error_invalid_argument"
    )
  }
  if (methods::is(x, "ExplicitLatent")) {
    x@meta$latent_units <- units
    methods::validObject(x)
    return(x)
  }
  if (inherits(x, "ImplicitLatent")) {
    x$meta$latent_units <- units
    return(x)
  }
  .encoder_cli_abort(
    paste0("Cannot attach latent units to class ", class(x)[1L], "."),
    class = "fmrilatent_error_unsupported_operation"
  )
}

#' @export
#' @rdname latent_units
setMethod("latent_units", "ExplicitLatent", function(x, ...) {
  .stored_latent_units(x) %||% .undeclared_latent_units()
})

#' @export
#' @rdname latent_units
setMethod("latent_units", "BilatLatentNeuroSurfaceVector", function(x, ...) {
  .stored_latent_units(x) %||%
    .combine_child_latent_units(list(left = x@left, right = x@right))
})

#' @export
#' @rdname latent_units
setMethod("latent_units", "BlockLatentNeuroVector", function(x, ...) {
  .stored_latent_units(x) %||% .combine_child_latent_units(x@blocks)
})

#' @export
#' @rdname latent_units
setMethod("latent_units", "ImplicitLatent", function(x, ...) {
  .stored_latent_units(x) %||% .undeclared_latent_units()
})

#' @export
print.fmrilatent_units <- function(x, ...) {
  x <- .validate_latent_units(x)
  cat("Latent units\n")
  cat("  status: ", x$status, "\n", sep = "")
  cat("  compatibility_id: ", x$compatibility_id, "\n", sep = "")
  if (identical(x$status, "declared")) {
    for (field in .latent_units_semantic_fields) {
      cat("  ", field, ": ", x[[field]], "\n", sep = "")
    }
  } else if (!is.null(x$reason)) {
    cat("  reason: ", x$reason, "\n", sep = "")
  }
  invisible(x)
}

#' @export
`$<-.fmrilatent_units` <- function(x, name, value) {
  .encoder_cli_abort(
    "latent-units records are immutable; construct a new record.",
    class = "fmrilatent_error_invalid_argument"
  )
}

#' @export
`[[<-.fmrilatent_units` <- function(x, i, value) {
  .encoder_cli_abort(
    "latent-units records are immutable; construct a new record.",
    class = "fmrilatent_error_invalid_argument"
  )
}

#' @export
`[<-.fmrilatent_units` <- function(x, i, value) {
  .encoder_cli_abort(
    "latent-units records are immutable; construct a new record.",
    class = "fmrilatent_error_invalid_argument"
  )
}

#' @export
`names<-.fmrilatent_units` <- function(x, value) {
  .encoder_cli_abort(
    "latent-units records are immutable; construct a new record.",
    class = "fmrilatent_error_invalid_argument"
  )
}
