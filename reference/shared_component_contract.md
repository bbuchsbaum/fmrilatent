# Build a method-neutral shared component contract

A shared component contract describes a reusable spatial component
matrix without committing to a particular decomposition method. The
contract records dimensions, domain/support metadata, a digest, and
optional advisory fields.

## Usage

``` r
shared_component_contract(
  loadings,
  id = NULL,
  family = "shared_component",
  domain_id = "",
  support = NULL,
  measure = NULL,
  meta = list()
)
```

## Arguments

- loadings:

  Matrix-like component loadings (\`features x components\`), a
  \`SharedReference\` to such a matrix, or a \`GroupDeltaLoadings\`
  object.

- id:

  Optional component id.

- family:

  Method-neutral or method-family label.

- domain_id:

  Stable domain identifier.

- support:

  Optional support descriptor whose cardinality must match
  \`nrow(loadings)\` when supplied.

- measure:

  Optional support measure. Vectors must have one value per feature;
  matrices must be square over features.

- meta:

  Optional advisory metadata.

## Value

A \`SharedComponentContract\` object.
