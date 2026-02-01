# Diffusion wavelet basis specification

Diffusion wavelet basis specification

## Usage

``` r
basis_diffusion_wavelet(
  target_rank = 2000L,
  oversample = 20L,
  threshold = 1e-05,
  max_scales = 1L,
  epsilon = NULL
)
```

## Arguments

- target_rank:

  Cap on retained components per scale (keeps runtime bounded).

- oversample:

  Oversampling for randomized range finder.

- threshold:

  Absolute value threshold to enforce sparsity in compressed ops.

- max_scales:

  Maximum diffusion scales to compute.

- epsilon:

  Optional precision (unused in capped-rank path; kept for API parity).
