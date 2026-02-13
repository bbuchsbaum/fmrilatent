#' Compatibility Profile for External Integrations
#'
#' Returns the active fmrilatent compatibility profile. This is used to
#' opt into strict behavior for external packages that require historical
#' semantics.
#'
#' @param profile Optional explicit profile name. If `NULL`, uses
#'   `getOption("fmrilatent.compat", "native")`.
#'
#' @return A single character string profile identifier.
#' @export
fmrilatent_compat_profile <- function(profile = NULL) {
  prof <- profile
  if (is.null(prof)) {
    prof <- getOption("fmrilatent.compat", "native")
  }
  if (is.null(prof) || !length(prof) || is.na(prof[[1]]) || !nzchar(prof[[1]])) {
    return("native")
  }
  as.character(prof[[1]])
}

#' @keywords internal
fmrilatent_is_neuroarchive_compat <- function(profile = NULL) {
  identical(fmrilatent_compat_profile(profile), "neuroarchive_0.1.1")
}

#' @keywords internal
fmrilatent_option_alias <- function(primary,
                                    aliases = character(),
                                    default = NULL,
                                    compat_profile = NULL) {
  val <- getOption(primary, NULL)
  if (!is.null(val)) return(val)

  if (fmrilatent_is_neuroarchive_compat(compat_profile) && length(aliases) > 0L) {
    for (nm in aliases) {
      alt <- getOption(nm, NULL)
      if (!is.null(alt)) return(alt)
    }
  }

  default
}
