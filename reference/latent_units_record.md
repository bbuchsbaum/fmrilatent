# Declare a Latent Units Contract

Constructs the unit and normalization record captured by \[encode()\]
and the transport encoders. Every semantic field is required: the
constructor does not guess response scaling or basis conventions from
numeric data.

## Usage

``` r
latent_units_record(
  response_scaling,
  coefficient_units,
  loading_normalization,
  loading_metric,
  analysis_coordinate_metric,
  sign_convention,
  notes = NULL
)
```

## Arguments

- response_scaling:

  Scaling applied to the response before encoding, such as
  \`"raw_signal"\`, \`"percent_signal_change"\`, or another precise
  caller-defined convention.

- coefficient_units:

  Physical or statistical units of analysis-coordinate coefficients.

- loading_normalization:

  Normalization convention for decoder/loading columns.

- loading_metric:

  Inner product or weighting under which loading normalization is
  defined.

- analysis_coordinate_metric:

  Metric convention used by analysis coordinates.

- sign_convention:

  Rule fixing coefficient/loading signs, or an explicit declaration that
  signs are arbitrary.

- notes:

  Optional human-readable note. Notes are integrity checked but do not
  change the compatibility ID.

## Value

An immutable declared \`fmrilatent_units\` record.
