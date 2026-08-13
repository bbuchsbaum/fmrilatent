# Create an in-session shared object reference

\`shared_reference()\` records a reusable object behind a stable
in-session identifier. It deliberately does not create files, archive
locators, HDF5 datasets, checksums, or persistent registries; those
responsibilities belong to archive packages such as neuroarchive.

## Usage

``` r
shared_reference(
  value = NULL,
  id = NULL,
  kind = "object",
  materialize = NULL,
  meta = list(),
  register = TRUE
)
```

## Arguments

- value:

  Optional object to cache immediately.

- id:

  Optional reference id. If omitted, a deterministic digest is used when
  \`value\` is present; otherwise a unique session id is generated.

- kind:

  Human-readable object kind.

- materialize:

  Optional zero-argument function used to create the value on first
  resolution.

- meta:

  Optional advisory metadata.

- register:

  If \`TRUE\`, cache \`value\` in the session registry.

## Value

A \`SharedReference\` descriptor.
