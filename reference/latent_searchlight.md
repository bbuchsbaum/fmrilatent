# Apply a user-defined function in latent space over neighborhoods

Runs a user-supplied function \`fun\` for each neighborhood using only
latent quantities. \`fun\` is called with arguments \`(B, L_V, M_V, idx,
...)\`, where: - \`B\` is the basis matrix (time x k) - \`L_V\` is the
loadings restricted to the neighborhood (\|V\| x k) - \`M_V = t(L_V) -
\`idx\` is the voxel indices of the neighborhood \`fun\` should return
any R object; results are collected in a list.

## Usage

``` r
latent_searchlight(basis, loadings, neighborhoods, fun, ...)
```

## Arguments

- basis:

  Matrix or BasisHandle (time x k) from a LatentNeuroVec.

- loadings:

  Matrix or LoadingsHandle (voxels x k) from a LatentNeuroVec.

- neighborhoods:

  List of integer vectors of voxel indices (mask order).

- fun:

  Function(B, L_V, M_V, idx, ...) returning a result per neighborhood.

- ...:

  Passed through to \`fun\`.

## Value

List of results, one per neighborhood.
