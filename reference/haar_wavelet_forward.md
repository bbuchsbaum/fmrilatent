# Forward Haar wavelet transform (mask-adaptive, Morton order)

Forward Haar wavelet transform (mask-adaptive, Morton order)

## Usage

``` r
haar_wavelet_forward(
  X,
  mask,
  levels = NULL,
  z_seed = 42L,
  threshold = list(type = "none", value = 0)
)
```

## Arguments

- X:

  Numeric matrix with time in rows and voxels (masked) in columns.

- mask:

  3D logical array or LogicalNeuroVol defining the mask.

- levels:

  Number of decomposition levels; defaults to
  ceil(log2(max(dim(mask)))).

- z_seed:

  Integer seed for Morton ordering (tie-breaking).

- threshold:

  List with fields \`type\` ("none" \| "absolute" \|
  "relative_to_root_std") and \`value\`.

## Value

List with elements \`coeff\` (list(root, detail)), \`meta\` (counts,
Morton hash, etc.).
