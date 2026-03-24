.onLoad <- function(libname, pkgname) {
  register_encoder("time_slepian", spec_time_slepian,
                   description = "Temporal DPSS/Slepian basis",
                   package = pkgname)
  register_encoder("time_dct", spec_time_dct,
                   description = "Temporal Discrete Cosine Transform basis",
                   package = pkgname)
  register_encoder("time_bspline", spec_time_bspline,
                   description = "Temporal B-spline basis",
                   package = pkgname)
  register_encoder("space_slepian", spec_space_slepian,
                   description = "Spatial Slepian eigenvectors on clustered reduction",
                   package = pkgname)
  register_encoder("space_pca", spec_space_pca,
                   description = "Spatial PCA eigenvectors within clusters/parcels",
                   package = pkgname)
  register_encoder("space_heat", spec_space_heat,
                   description = "Spatial graph diffusion heat-wavelet basis",
                   package = pkgname)
  register_encoder("space_hrbf", spec_space_hrbf,
                   description = "Spatial hierarchical radial basis functions",
                   package = pkgname)
  register_encoder("space_wavelet_active", spec_space_wavelet_active,
                   description = "Spatial lifting wavelet (active pencil) basis",
                   package = pkgname)
  register_encoder("st", spec_st,
                   description = "Separable spatiotemporal encoding (time x space)",
                   package = pkgname)
  register_encoder("hierarchical", spec_hierarchical_template,
                   description = "Hierarchical Laplacian template basis",
                   package = pkgname)
  register_encoder("space_parcel", spec_space_parcel,
                   description = "Shared parcel basis template (project-only)",
                   package = pkgname)
}
