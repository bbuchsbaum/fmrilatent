# List registered encoders

Returns a data.frame describing all registered encoder families,
including both built-in encoders and those added by external packages.

## Usage

``` r
list_encoders()
```

## Value

A data.frame with columns: `family`, `description`, `package`.

## Examples

``` r
list_encoders()
#>                  family                                         description
#> 1          hierarchical               Hierarchical Laplacian template basis
#> 2            space_heat          Spatial graph diffusion heat-wavelet basis
#> 3            space_hrbf         Spatial hierarchical radial basis functions
#> 4          space_parcel         Shared parcel basis template (project-only)
#> 5             space_pca    Spatial PCA eigenvectors within clusters/parcels
#> 6         space_slepian Spatial Slepian eigenvectors on clustered reduction
#> 7  space_wavelet_active       Spatial lifting wavelet (active pencil) basis
#> 8                    st    Separable spatiotemporal encoding (time x space)
#> 9          time_bspline                             Temporal B-spline basis
#> 10             time_dct            Temporal Discrete Cosine Transform basis
#> 11         time_slepian                         Temporal DPSS/Slepian basis
#>       package
#> 1  fmrilatent
#> 2  fmrilatent
#> 3  fmrilatent
#> 4  fmrilatent
#> 5  fmrilatent
#> 6  fmrilatent
#> 7  fmrilatent
#> 8  fmrilatent
#> 9  fmrilatent
#> 10 fmrilatent
#> 11 fmrilatent
```
