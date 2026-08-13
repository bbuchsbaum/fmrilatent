# Latent object class union

A common S4 dispatch target spanning the package's two latent branches:
explicit S4 latent objects that inherit from `ExplicitLatent`, and
decoder-backed S3 objects registered through `ImplicitLatent` (including
`TransportLatent`). This union is intentionally broader than
`ExplicitLatent`: use [`is_explicit_latent()`](is_explicit_latent.md)
when the caller specifically requires stored basis/loadings, and use the
`Latent` union for cross-cutting S4 methods that can operate on either
explicit or decoder-backed latent objects.
