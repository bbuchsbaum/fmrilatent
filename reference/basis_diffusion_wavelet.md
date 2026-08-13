# Diffusion wavelet basis specification

Diffusion wavelet basis specification

## Usage

``` r
basis_diffusion_wavelet(
  target_rank = 2000L,
  oversample = 20L,
  threshold = NULL,
  max_scales = 1L,
  epsilon = NULL,
  seed = 1L,
  sparsify_eps = NULL
)
```

## Arguments

- target_rank:

  Cap on retained components per scale (keeps runtime bounded).

- oversample:

  Oversampling for randomized range finder.

- threshold:

  Deprecated alias for \`sparsify_eps\`.

- max_scales:

  Maximum diffusion scales to compute.

- epsilon:

  Optional precision (unused in capped-rank path; kept for API parity).

- seed:

  Optional integer seed for deterministic randomized range finding.

- sparsify_eps:

  Absolute value threshold to enforce sparsity in compressed operators.
  Stored as \`threshold\` for compatibility.

## Value

A \`spec_diffusion_wavelet\` basis specification descriptor.
