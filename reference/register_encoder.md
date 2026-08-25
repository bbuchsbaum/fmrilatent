# Register an encoder family

Registers a spec constructor so it can be discovered via
[`list_encoders`](list_encoders.md) and retrieved via
[`get_encoder`](get_encoder.md).

## Usage

``` r
register_encoder(family, spec_fn, description = "", package = "")
```

## Arguments

- family:

  Character string identifying the encoder family (e.g. "time_slepian").

- spec_fn:

  The spec constructor function (e.g. `spec_time_slepian`).

- description:

  Character string describing the encoder.

- package:

  Character string naming the package that provides this encoder.

## Value

Invisibly, `TRUE`.

## Details

External packages typically register their encoders in their `.onLoad`
hook so the registry is populated when the package is loaded:


    .onLoad <- function(libname, pkgname) {
      if (requireNamespace("fmrilatent", quietly = TRUE)) {
        fmrilatent::register_encoder(
          "my_encoder", spec_my_encoder,
          description = "My custom encoder",
          package = pkgname
        )
      }
    }

The S3 method itself is exported via the package's NAMESPACE in the
usual way (`S3method(encode_spec, spec_my_encoder)`).

Registering the same family twice issues a warning and overwrites the
previous entry.

## Dispatch model

The registry is for **discovery and introspection only**. It does *not*
participate in actual encoding. [`encode`](encode.md) routes work
through the S3 generic `encode_spec`, which dispatches on the class of
the spec object. That means an external package contributing a new
encoder family must do **both**:

1.  Define a spec constructor whose return value carries a distinctive
    class (e.g. `class(spec) <- c("spec_my_encoder", "list")`).

2.  Define an S3 method `encode_spec.spec_my_encoder(x, spec, ...)` that
    builds the actual `LatentNeuroVec` (or other latent type).

3.  Optionally call `register_encoder()` so that
    [`list_encoders()`](list_encoders.md) surfaces the family for users.

Calling `register_encoder()` alone is **not enough**: without an
`encode_spec.spec_*` method, [`encode()`](encode.md) will fall through
to `encode_spec.default` and raise an error.

Transport-backed AWPT is an intentional parallel API, not a
registry-dispatch family.
[`basis_awpt_wavelet()`](basis_awpt_wavelet.md) describes a shared
template basis, while subject fitting also needs a template/basis asset
plus subject field or observation operators. Use
[`encode_awpt`](encode_awpt.md) or
[`encode_operator`](encode_operator.md) for those fits; do not expect
`register_encoder()` or [`encode()`](encode.md) to route AWPT subject
encoding.

## See also

[`encode`](encode.md), [`list_encoders`](list_encoders.md),
[`get_encoder`](get_encoder.md)

## Examples

``` r
register_encoder("test_enc", identity, "A test encoder", "mypkg")
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
#> 9              test_enc                                      A test encoder
#> 10         time_bspline                             Temporal B-spline basis
#> 11             time_dct            Temporal Discrete Cosine Transform basis
#> 12         time_slepian                         Temporal DPSS/Slepian basis
#>       package
#> 1  fmrilatent
#> 2  fmrilatent
#> 3  fmrilatent
#> 4  fmrilatent
#> 5  fmrilatent
#> 6  fmrilatent
#> 7  fmrilatent
#> 8  fmrilatent
#> 9       mypkg
#> 10 fmrilatent
#> 11 fmrilatent
#> 12 fmrilatent
get_encoder("test_enc")
#> $spec_fn
#> function (x) 
#> x
#> <bytecode: 0x55aef8df29a8>
#> <environment: namespace:base>
#> 
#> $description
#> [1] "A test encoder"
#> 
#> $package
#> [1] "mypkg"
#> 
```
