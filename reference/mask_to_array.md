# Convert mask to array

Safely converts a `LogicalNeuroVol` or array-like mask to a plain
logical array, with informative error messages on failure.

## Usage

``` r
mask_to_array(mask, location = "unknown function")
```

## Arguments

- mask:

  A `LogicalNeuroVol` or logical array.

- location:

  Character string used in error messages to identify the caller.

## Value

A logical array.
