# Validate the reusable template protocol

Checks that a template asset exposes the method-neutral surface expected
by shared-structure consumers: loadings, rank, support/domain, measure,
roughness, decoder, and projection behavior.

## Usage

``` r
validate_template_protocol(template, error = TRUE)
```

## Arguments

- template:

  Template object.

- error:

  If \`TRUE\`, throw on failure; otherwise return \`FALSE\`.

## Value

Invisibly, a manifest describing the verified protocol.
