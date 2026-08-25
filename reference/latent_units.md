# Describe Latent Response and Coordinate Units

Returns the truthful unit and normalization contract retained on a
latent object. Legacy objects that predate encode-time capture return a
structured record with \`status = "undeclared"\`; no units are inferred
from numeric values or representation class.

## Usage

``` r
latent_units(x, ...)

# S4 method for class 'ExplicitLatent'
latent_units(x, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
latent_units(x, ...)

# S4 method for class 'BlockLatentNeuroVector'
latent_units(x, ...)

# S4 method for class 'ImplicitLatent'
latent_units(x, ...)
```

## Arguments

- x:

  A latent object.

- ...:

  Additional arguments passed to methods.

## Value

An \`fmrilatent_units\` record with a deterministic \`compatibility_id\`
and human-readable fields.
