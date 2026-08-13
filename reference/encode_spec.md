# Dispatch standard encoding based on spec type

Dispatch standard encoding based on spec type

## Usage

``` r
encode_spec(x, spec, ...)
```

## Arguments

- x:

  Data matrix.

- spec:

  Spec object.

- ...:

  Additional arguments passed to methods.

## Value

Encoded representation.

## Details

\`encode_spec()\` is the S3 dispatch seam for standard \`encode()\`
specs: temporal, spatial, hierarchical, parcel, and spatiotemporal
specs. AWPT's \`basis_awpt_wavelet()\` object is a template-construction
spec, not a complete subject-encoding spec, because AWPT fitting
additionally needs a shared template/basis asset and subject field or
observation operators. AWPT therefore uses the parallel transport API
\[encode_awpt()\] / \[encode_operator()\].
