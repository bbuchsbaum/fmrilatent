# Simulate data with BOLDZip-SR carrier, texture, and event structure

Simulate data with BOLDZip-SR carrier, texture, and event structure

## Usage

``` r
boldzip_sr_simulate(
  n_voxels = 40L,
  n_time = 80L,
  k_carriers = 3L,
  q_texture = 1L,
  n_events = 8L,
  noise_sd = 0.05,
  seed = NULL
)
```

## Arguments

- n_voxels:

  Number of spatial rows.

- n_time:

  Number of time points.

- k_carriers:

  Number of latent carrier time courses.

- q_texture:

  Number of non-zero carrier loadings per voxel.

- n_events:

  Number of paired impulse events to add.

- noise_sd:

  Independent Gaussian noise standard deviation.

- seed:

  Optional random seed.

## Value

List containing \`X\`, \`mu\`, \`carriers\`, \`texture\`, and
\`events\`.
