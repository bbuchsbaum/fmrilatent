# Compatibility Profile for External Integrations

Returns the active fmrilatent compatibility profile. This is used to opt
into strict behavior for external packages that require historical
semantics.

## Usage

``` r
fmrilatent_compat_profile(profile = NULL)
```

## Arguments

- profile:

  Optional explicit profile name. If \`NULL\`, uses
  \`getOption("fmrilatent.compat", "native")\`.

## Value

A single character string profile identifier.
