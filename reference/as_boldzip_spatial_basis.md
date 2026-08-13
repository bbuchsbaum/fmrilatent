# Coerce a spatial object to a BOLDZip-SR spatial basis

This helper lets the standalone BOLDZip-SR codec consume matrix-like
shared basis assets without registering BOLDZip as an \`encode()\`
family. Matrix and template inputs are used as the detail basis by
default and are orthonormalized because BOLDZip projection currently
uses the transpose as the analysis operator.

## Usage

``` r
as_boldzip_spatial_basis(x, ...)
```

## Arguments

- x:

  A \`BoldZipSRSpatialBasis\`, matrix-like object, shared reference, or
  template object supporting \[template_loadings()\].

- ...:

  Additional arguments reserved for methods. The default method accepts
  \`label\`, \`role\`, and \`orthonormalize\`.

## Value

A \`BoldZipSRSpatialBasis\` object.
