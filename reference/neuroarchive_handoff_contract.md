# Create the fmrilatent-to-neuroarchive handoff contract

The handoff contract records the boundary between in-memory
representation responsibilities owned by fmrilatent and persistent
archive responsibilities owned by neuroarchive. It is intentionally a
manifest, not a lazy archive reader or file locator.

## Usage

``` r
neuroarchive_handoff_contract(
  representation = NULL,
  components = list(),
  templates = list(),
  references = list(),
  meta = list()
)
```

## Arguments

- representation:

  Optional fmrilatent representation object.

- components:

  Optional list of shared component contracts.

- templates:

  Optional list of reusable template assets.

- references:

  Optional list of in-session shared references.

- meta:

  Optional advisory metadata.

## Value

A \`NeuroarchiveHandoffContract\`.
