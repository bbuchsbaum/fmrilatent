# Hierarchical encode_spec methods

#' @include encode.R
NULL

#' Create hierarchical template spec
#' @param template HierarchicalBasisTemplate object
#' @param template_file Path to saved template file
#' @return Spec object of class spec_hierarchical
#' @export
spec_hierarchical_template <- function(template = NULL, template_file = NULL) {
  if (is.null(template) && is.null(template_file)) {
    .encoder_cli_abort(
      "Provide either a template object or template_file",
      class = "fmrilatent_error_invalid_hierarchical_template"
    )
  }
  structure(list(template = template, template_file = template_file), class = "spec_hierarchical")
}

#' @exportS3Method
encode_spec.spec_hierarchical <- function(x, spec, mask, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_hierarchical")
  tmpl <- spec$template
  if (is.null(tmpl)) {
    if (is.null(spec$template_file)) {
      .encoder_cli_abort(
        "template_file is required when template is NULL",
        class = "fmrilatent_error_invalid_hierarchical_template"
      )
    }
    tmpl <- load_hierarchical_template(spec$template_file)
  }
  .assert_template_mask_match(mask, template_mask(tmpl), "encode_spec.spec_hierarchical")
  encode_hierarchical(x, tmpl, label = label, mask = mask, materialize = materialize)
}
