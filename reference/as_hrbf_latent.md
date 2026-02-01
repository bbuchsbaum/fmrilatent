# Attach HRBF metadata to an existing LatentNeuroVec

Attach HRBF metadata to an existing LatentNeuroVec

## Usage

``` r
as_hrbf_latent(lvec, params, centres = NULL, sigmas = NULL)
```

## Arguments

- lvec:

  A LatentNeuroVec object

- params:

  HRBF parameters list

- centres:

  Optional matrix of HRBF centres

- sigmas:

  Optional vector of HRBF sigmas

## Value

The LatentNeuroVec with HRBF metadata attached
