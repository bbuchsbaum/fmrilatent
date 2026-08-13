# Spatial parcel-basis spec (shared/template-based)

Creates a spec for projecting data onto a pre-built
`"ParcelBasisTemplate"`. The loadings are fixed (shared across
subjects); only the temporal scores and per-subject offset vary.

## Usage

``` r
spec_space_parcel(template)
```

## Arguments

- template:

  A `"ParcelBasisTemplate"` object built by
  [`parcel_basis_template()`](parcel_basis_template.md).

## Value

A `spec_space_parcel` object for [`encode()`](encode.md).

## See also

[`parcel_basis_template`](parcel_basis_template.md),
[`encode`](encode.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tmpl <- parcel_basis_template(atlas, basis_slepian(k = 8))
lvec_s1 <- encode(bold_s1, spec_space_parcel(tmpl))
lvec_s2 <- encode(bold_s2, spec_space_parcel(tmpl))
# Same loadings; different basis (scores) and offset
} # }
```
