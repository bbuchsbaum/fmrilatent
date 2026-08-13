# Create a LoadingsHandle for spatial Slepians (graph Laplacian)

Create a LoadingsHandle for spatial Slepians (graph Laplacian)

## Usage

``` r
slepian_spatial_loadings_handle(
  reduction,
  basis_spec = basis_slepian(),
  data = NULL,
  k_neighbors = 6L,
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

- k_neighbors:

  Number of neighbors used for local graph construction when
  materializing the lifted basis.

- id:

  Optional registry id; generated if NULL.

- label:

  Optional label.

## Value

A `LoadingsHandle`.

## Details

This constructor lifts the spatial dictionary eagerly so the returned
handle records the realized dimensions and registers a fingerprinted
cache entry. Repeated constructor calls may therefore recompute the lift
even when later \`loadings_mat()\` calls can reuse the registry cache.
