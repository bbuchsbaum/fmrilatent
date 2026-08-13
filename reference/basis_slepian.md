# Create a Slepian basis specification

\`basis_slepian()\` creates a lightweight descriptor for graph/Slepian
basis construction during \`lift()\` or spatial encoding. It records the
requested component count and Slepian flavor; the actual basis is
computed later by a \`lift()\` method for the supplied reduction.

## Usage

``` r
basis_slepian(k = 3, type = "laplacian")
```

## Arguments

- k:

  Positive integer number of Slepian components.

- type:

  Character scalar naming the Slepian basis flavor. The built-in spatial
  methods use \`"laplacian"\`.

## Value

A list with class \`spec_slepian\` containing \`k\` and \`type\`.
