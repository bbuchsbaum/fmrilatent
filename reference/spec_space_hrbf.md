# Spatial HRBF spec

Spatial HRBF spec

## Usage

``` r
spec_space_hrbf(
  params = list(),
  sigma0 = NULL,
  levels = NULL,
  radius_factor = NULL,
  num_extra_fine_levels = NULL,
  kernel_type = NULL,
  kernel_type_fine_levels = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- params:

  Optional named list for advanced or compatibility HRBF parameters.
  Explicit formals override entries with the same names.

- sigma0:

  Base kernel width. \`NULL\` uses the HRBF default (\`6\`).

- levels:

  Number of dyadic resolution levels. \`NULL\` uses the HRBF default
  (\`3L\`).

- radius_factor:

  Atom support radius as a multiple of level width. \`NULL\` uses the
  HRBF default (\`2.5\`).

- num_extra_fine_levels:

  Optional extra fine levels. \`NULL\` uses the HRBF default (\`0L\`).

- kernel_type:

  Radial kernel family. \`NULL\` uses \`"gaussian"\`.

- kernel_type_fine_levels:

  Kernel family for extra fine levels. \`NULL\` uses the HRBF default.

- seed:

  Random seed. \`NULL\` uses the HRBF default (\`1L\`).

- ...:

  Additional named HRBF parameters for advanced use.

## Value

A \`spec_space_hrbf\` object.
