# Get a registered encoder

Retrieves the registration entry for a given encoder family.

## Usage

``` r
get_encoder(family)
```

## Arguments

- family:

  Character string identifying the encoder family.

## Value

A list with elements `spec_fn`, `description`, and `package`.

## Examples

``` r
enc <- get_encoder("time_dct")
spec <- enc$spec_fn(k = 5)
```
