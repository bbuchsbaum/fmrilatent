# Create a LoadingsHandle for spatial Slepians (graph Laplacian)

Create a LoadingsHandle for spatial Slepians (graph Laplacian)

## Usage

``` r
slepian_spatial_loadings_handle(
  reduction,
  basis_spec = basis_slepian(),
  data = NULL,
  id = NULL,
  label = "slepian-spatial"
)
```

## Arguments

- reduction:

  Graph reduction (e.g., ClusterReduction).

- basis_spec:

  Slepian basis spec (from \`basis_slepian()\`).

- data:

  Optional data passed to \`lift()\` (if needed).

- id:

  Optional registry id; generated if NULL.

- label:

  Optional label.

## Value

A `LoadingsHandle`.
