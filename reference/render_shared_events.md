# Render sparse events with a shared event dictionary

Render sparse events with a shared event dictionary

## Usage

``` r
render_shared_events(dictionary, events, n_atoms, n_time = NULL)
```

## Arguments

- dictionary:

  A \`SharedEventDictionary\`.

- events:

  Data frame with columns \`atom\`, \`time\`, \`amplitude\`, and
  optional \`shape_id\`.

- n_atoms:

  Number of event rows/atoms.

- n_time:

  Number of time points. Defaults to the dictionary value.

## Value

Sparse \`n_atoms x n_time\` event matrix.
