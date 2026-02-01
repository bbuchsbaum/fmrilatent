# Run spectral+Ward hierarchical clustering on a parcel graph

Produces an hclust over parcels using a spectral embedding of the
similarity matrix, then Ward linkage. Real data glue: supply W from
\`parcel_similarity_matrix()\`, \`hemi\` from parcel hemisphere labels,
and \`network\` from Yeo17 labels. Penalties bias merges to stay within
hemisphere/network until higher in the tree.

## Usage

``` r
spectral_ward_hclust(W, k_embed = 3, hemi = NULL, network = NULL)
```

## Arguments

- W:

  Similarity matrix (symmetric).

- k_embed:

  Embedding dimension (2 or 3).

- hemi:

  Optional factor/character vector (n_parc) to enforce within-hemisphere
  merges first.

- network:

  Optional factor/character vector (n_parc) (e.g., Yeo17) to enforce
  within-network merges until relaxed.

## Value

hclust object.
