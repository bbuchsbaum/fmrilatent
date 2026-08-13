# Get a decoder view for a latent object

Get a decoder view for a latent object

## Usage

``` r
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'ImplicitLatent'
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'LatentNeuroVec'
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
decoder(
  x,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)
```

## Arguments

- x:

  A latent object.

- space:

  Output space to decode into.

- coordinates:

  Coordinate system consumed by the decoder.

- ...:

  Additional arguments passed to methods.

## Value

Decoder view object.
