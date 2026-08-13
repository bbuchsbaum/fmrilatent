# Build a matrix spatial basis descriptor for BOLDZip-SR

Build a matrix spatial basis descriptor for BOLDZip-SR

## Usage

``` r
boldzip_spatial_basis(
  phi_c = NULL,
  phi_d = NULL,
  label = NULL,
  basis_asset = NULL
)
```

## Arguments

- phi_c:

  Optional coarse basis matrix with rows equal to voxels.

- phi_d:

  Optional detail basis matrix with rows equal to voxels. If \`NULL\`,
  the detail basis is the identity basis.

- label:

  Optional label stored in metadata.

- basis_asset:

  Optional source template or shared-basis object used to build this
  descriptor.

## Value

A \`BoldZipSRSpatialBasis\` object.
