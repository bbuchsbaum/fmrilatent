# Graph reduction scaffolds (abstract)

These classes describe how voxels are grouped or coarsened before a
basis is computed and lifted back to ambient space. Implementations of
\`lift()\` for specific combinations (e.g., supervoxel + Slepian, parcel
PCA) live in external packages or downstream code; fmrilatent ships only
the contracts.

## Slots

- `mask`:

  \`LogicalNeuroVol\` defining the ambient domain.

- `info`:

  Optional list for implementation-specific metadata.
