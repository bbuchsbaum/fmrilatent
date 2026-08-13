# Construct a reusable event shape dictionary

Construct a reusable event shape dictionary

## Usage

``` r
shared_event_dictionary(
  shapes = list(impulse = 1),
  n_time = NULL,
  meta = list()
)
```

## Arguments

- shapes:

  Named list of numeric shape vectors.

- n_time:

  Optional temporal length used for validation.

- meta:

  Optional advisory metadata.

## Value

A \`SharedEventDictionary\`.
