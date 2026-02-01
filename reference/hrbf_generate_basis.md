# Hierarchical radial basis functions (HRBF) for latent fMRI

Utilities to generate analytic HRBF bases, project data, reconstruct,
and build \`LatentNeuroVec\` objects. Parameters are kept simple and
in-R only (no descriptors or HDF5).

## Usage

``` r
hrbf_generate_basis(params, mask)

hrbf_project_matrix(X, mask, params)

hrbf_reconstruct_matrix(coeff, mask, params)
```

## Arguments

- params:

  List with fields: - \`sigma0\` (numeric, default 6) - \`levels\`
  (integer, default 3) - \`radius_factor\` (numeric, default 2.5) -
  \`num_extra_fine_levels\` (integer, default 0) - \`kernel_type\`
  (\\gaussian\\ or \\wendland_c6\\, alias \\wendland_c4\\) - \`seed\`
  (integer) for deterministic Poisson sampling

- mask:

  \`LogicalNeuroVol\` mask defining voxel locations.

- X:

  Numeric matrix with time in rows and voxels in columns.

- coeff:

  Coefficient matrix with rows = time points.

## Value

For \`hrbf_generate_basis\`, a sparse matrix with one row per HRBF atom
and columns matching mask voxels.
