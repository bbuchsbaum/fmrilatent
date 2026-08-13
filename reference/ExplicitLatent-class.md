# Virtual ExplicitLatent S4 marker

Common parent for all latent objects that store data as an explicit
`basis x loadings + offset` factorisation. Concrete subclasses
(`LatentNeuroVec`, `LatentNeuroSurfaceVector`, `BlockLatentNeuroVector`,
`BilatLatentNeuroSurfaceVector`) carry their domain-specific slots and
inherit the [`is_explicit_latent()`](is_explicit_latent.md) predicate
from this class.

## Details

Decoder-backed latents (`ImplicitLatent` and its subclass
`TransportLatent`) do NOT inherit from `ExplicitLatent` — they are S3
and the predicate returns `FALSE` for them.
