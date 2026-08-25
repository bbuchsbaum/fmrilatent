# Package index

## All functions

- [`BilatLatentNeuroSurfaceVector-class`](BilatLatentNeuroSurfaceVector-class.md)
  : BilatLatentNeuroSurfaceVector Class
- [`BilatLatentNeuroSurfaceVector()`](BilatLatentNeuroSurfaceVector.md)
  : Construct a bilateral surface latent object
- [`BlockLatentNeuroVector-class`](BlockLatentNeuroVector-class.md) :
  BlockLatentNeuroVector Class
- [`BlockLatentNeuroVector()`](BlockLatentNeuroVector.md) : Construct a
  block-domain latent object
- [`ClusterReduction-class`](ClusterReduction-class.md) : Cluster-based
  reduction (e.g., supervoxels or atlas)
- [`CoarsenedReduction-class`](CoarsenedReduction-class.md) : Coarsened
  graph reduction (e.g., prolongation from coarse to fine)
- [`GraphReduction-class`](GraphReduction-class.md) : Graph reduction
  scaffolds (abstract)
- [`HierarchicalBasisTemplate-class`](HierarchicalBasisTemplate-class.md)
  : HierarchicalBasisTemplate Class
- [`LatentNeuroSurfaceVector-class`](LatentNeuroSurfaceVector-class.md)
  : LatentNeuroSurfaceVector Class
- [`LatentNeuroSurfaceVector()`](LatentNeuroSurfaceVector.md) :
  Construct a LatentNeuroSurfaceVector
- [`LatentNeuroVec-class`](LatentNeuroVec-class.md) : LatentNeuroVec
  Class
- [`LatentNeuroVec()`](LatentNeuroVec.md) : Create a Latent Space
  Representation of Neuroimaging Data
- [`analysis_transform()`](analysis_transform.md) : Describe the
  transform from raw to analysis coordinates
- [`as.array(`*`<LatentNeuroVec>`*`)`](as.array-LatentNeuroVec-method.md)
  : Reconstruct LatentNeuroVec as a 4D array
- [`as.array(`*`<ImplicitLatent>`*`)`](as.array.ImplicitLatent.md) :
  Reconstruct an ImplicitLatent as an array
- [`as.matrix(`*`<LatentNeuroSurfaceVector>`*`)`](as.matrix-LatentNeuroVec-method.md)
  [`as.matrix(`*`<BilatLatentNeuroSurfaceVector>`*`)`](as.matrix-LatentNeuroVec-method.md)
  [`as.matrix(`*`<BlockLatentNeuroVector>`*`)`](as.matrix-LatentNeuroVec-method.md)
  [`as.matrix(`*`<LatentNeuroVec>`*`)`](as.matrix-LatentNeuroVec-method.md)
  : Reconstruct LatentNeuroVec as a matrix (time x voxels)
- [`as.matrix(`*`<GroupDeltaLoadings>`*`)`](as.matrix.GroupDeltaLoadings.md)
  : Coerce group-plus-delta loadings to a dense matrix
- [`as.matrix(`*`<HaarLatent>`*`)`](as.matrix.HaarLatent.md) : Convert
  HaarLatent to matrix
- [`as.matrix(`*`<ImplicitLatent>`*`)`](as.matrix.ImplicitLatent.md) :
  Reconstruct an ImplicitLatent as a matrix
- [`as_boldzip_spatial_basis()`](as_boldzip_spatial_basis.md) : Coerce a
  spatial object to a BOLDZip-SR spatial basis
- [`as_cluster_reduction()`](as_cluster_reduction.md) : Convert a
  ClusteredNeuroVol to a ClusterReduction
- [`as_haar_latent()`](as_haar_latent.md) : Convert to HaarLatent class
- [`as_hrbf_latent()`](as_hrbf_latent.md) : Attach HRBF metadata to an
  existing LatentNeuroVec
- [`as_implicit_latent(`*`<BoldZipSR>`*`)`](as_implicit_latent.BoldZipSR.md)
  : Coerce a BOLDZip-SR payload to an ImplicitLatent
- [`as_implicit_latent()`](as_implicit_latent.md) : Coerce an object to
  an ImplicitLatent decoder
- [`as_portable_linear_map()`](as_portable_linear_map.md) : Coerce an
  operator to the portable linear-map contract
- [`awpt_basis_template()`](awpt_basis_template.md) : Build an AWPT
  basis template
- [`awpt_mean_conductance()`](awpt_mean_conductance.md) : Average
  subject conductance matrices on a shared template graph
- [`awpt_operator_from_conductance()`](awpt_operator_from_conductance.md)
  : Build an AWPT field operator from a conductance matrix
- [`awpt_operator_from_subject_conductances()`](awpt_operator_from_subject_conductances.md)
  : Build an AWPT field operator from subject conductance summaries
- [`awpt_surface_basis_template()`](awpt_surface_basis_template.md) :
  Build a surface AWPT basis template
- [`basis()`](basis-methods.md) : Get the basis matrix (temporal
  components)
- [`basis_asset()`](basis_asset.md) : Extract the basis asset behind a
  latent object
- [`basis_awpt_wavelet()`](basis_awpt_wavelet.md) : AWPT wave-packet
  basis specification
- [`basis_decoder()`](basis_decoder.md) : Extract a basis decoder from a
  template asset
- [`basis_diffusion_wavelet()`](basis_diffusion_wavelet.md) : Diffusion
  wavelet basis specification
- [`basis_flat()`](basis_flat.md) : Create a flat basis specification
- [`basis_heat_wavelet()`](basis_heat_wavelet.md) : Heat wavelet basis
  specification
- [`basis_pca()`](basis_pca.md) : Create a PCA basis specification
- [`basis_slepian()`](basis_slepian.md) : Create a Slepian basis
  specification
- [`benchmark_roundtrip()`](benchmark_roundtrip.md) : Benchmark
  encode/decode round-trips
- [`boldzip_events()`](boldzip_events.md) : Construct sparse innovation
  event settings for BOLDZip-SR
- [`boldzip_graph_spatial_basis()`](boldzip_graph_spatial_basis.md) :
  Build a spectral graph spatial basis for BOLDZip-SR
- [`boldzip_parcel_reconstruct()`](boldzip_parcel_reconstruct.md) :
  Reconstruct a matrix from parcel-average time series
- [`boldzip_quantization()`](boldzip_quantization.md) : Construct
  reliability-aware quantization settings for BOLDZip-SR
- [`boldzip_reliability()`](boldzip_reliability.md) : Construct
  split-half reliability settings for BOLDZip-SR
- [`boldzip_spatial_basis()`](boldzip_spatial_basis.md) : Build a matrix
  spatial basis descriptor for BOLDZip-SR
- [`boldzip_sr_decode()`](boldzip_sr_decode.md) : Decode an experimental
  BOLDZip-SR object
- [`boldzip_sr_encode()`](boldzip_sr_encode.md) : Encode a matrix with
  experimental BOLDZip-SR compression
- [`boldzip_sr_payload_summary()`](boldzip_sr_payload_summary.md) :
  Summarize an experimental BOLDZip-SR payload
- [`boldzip_sr_simulate()`](boldzip_sr_simulate.md) : Simulate data with
  BOLDZip-SR carrier, texture, and event structure
- [`boldzip_svd_reconstruct()`](boldzip_svd_reconstruct.md) :
  Reconstruct a matrix from a low-rank SVD baseline
- [`bspline_basis_handle()`](bspline_basis_handle.md) : Create a
  BasisHandle for a B-spline temporal basis
- [`build_hierarchical_template()`](build_hierarchical_template.md) :
  Build a hierarchical Laplacian template (offline)
- [`build_schaefer_hierarchical_template()`](build_schaefer_hierarchical_template.md)
  : Build hierarchical template from Schaefer surface atlas
- [`build_schaefer_levels()`](build_schaefer_levels.md) : Build
  hierarchical parcellation levels from Schaefer surface atlas
- [`coef_metric()`](coef_metric.md) : Extract coefficient-space metric
  from a latent object
- [`coef_time()`](coef_time.md) : Extract coefficient time series from a
  latent object
- [`compare_boldzip_sr()`](compare_boldzip_sr.md) : Compare BOLDZip-SR
  against simple reconstruction baselines
- [`concat(`*`<LatentNeuroVec>`*`,`*`<LatentNeuroVec>`*`)`](concat-methods.md)
  : Concatenate LatentNeuroVec Objects
- [`cut_hclust_nested()`](cut_hclust_nested.md) : Cut an hclust into
  nested label vectors
- [`dct_basis_handle()`](dct_basis_handle.md) : Create a BasisHandle for
  a DCT temporal basis
- [`decode_coefficients()`](decode_coefficients.md) : Decode
  coefficient-space vectors into an output space
- [`decode_covariance()`](decode_covariance.md) : Push coefficient
  covariance into an output space
- [`decoder()`](decoder.md) : Get a decoder view for a latent object
- [`diffusion_wavelet_latent()`](diffusion_wavelet_latent.md) :
  Diffusion wavelet latent constructor (explicit basis)
- [`diffusion_wavelet_loadings_handle()`](diffusion_wavelet_loadings_handle.md)
  : Construct a shared LoadingsHandle via diffusion-wavelet lifting
- [`dpss_time_basis()`](dpss_time_basis.md) : DPSS temporal basis
  (Slepian sequences)
- [`encode()`](encode.md) : Encode data into a latent representation
  using a spec
- [`encode_awpt()`](encode_awpt.md) : Encode data using the AWPT model
- [`encode_hierarchical()`](encode_hierarchical.md) : Encode data using
  a hierarchical template
- [`encode_operator()`](encode_operator.md) : Encode data using a shared
  basis asset and subject field operator
- [`encode_spec()`](encode_spec.md) : Dispatch standard encoding based
  on spec type
- [`encode_transport()`](encode_transport.md) : Encode data using
  transport-backed latent semantics
- [`evaluate_boldzip_sr()`](evaluate_boldzip_sr.md) : Evaluate
  BOLDZip-SR reconstruction quality
- [`` `[[`( ``*`<LatentNeuroVec>`*`,`*`<numeric>`*`)`](extract-methods.md)
  [`` `[`( ``*`<LatentNeuroVec>`*`,`*`<numeric>`*`,`*`<numeric>`*`,`*`<ANY>`*`)`](extract-methods.md)
  [`` `[`( ``*`<LatentNeuroVec>`*`,`*`<matrix>`*`,`*`<missing>`*`,`*`<ANY>`*`)`](extract-methods.md)
  [`` `[`( ``*`<LatentNeuroVec>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](extract-methods.md)
  : Extract Elements from LatentNeuroVec
- [`fmrilatent_compat_profile()`](fmrilatent_compat_profile.md) :
  Compatibility Profile for External Integrations
- [`fmrilatent_registry_clear()`](fmrilatent_registry_clear.md) : Clear
  the fmrilatent handle registry
- [`fmrilatent_registry_enable()`](fmrilatent_registry_enable.md)
  [`fmrilatent_registry_disable()`](fmrilatent_registry_enable.md)
  [`fmrilatent_registry_enabled()`](fmrilatent_registry_enable.md) :
  Enable or disable the fmrilatent handle registry
- [`fmrilatent_registry_list()`](fmrilatent_registry_list.md) : List
  entries in the fmrilatent handle registry
- [`fmrilatent_registry_stats()`](fmrilatent_registry_stats.md) : Get
  registry statistics
- [`fmrilatent_test_data()`](fmrilatent_test_data.md) : Generate test
  data for encoder development
- [`get_encoder()`](get_encoder.md) : Get a registered encoder
- [`group_delta_loadings()`](group_delta_loadings.md) : Represent
  subject loadings as group shared loadings plus a delta
- [`haar_latent()`](haar_latent.md) : Build Haar latent representation
- [`haar_meta()`](haar_meta.md) : Get metadata from Haar latent object
- [`haar_wavelet_forward()`](haar_wavelet_forward.md) : Forward Haar
  wavelet transform (mask-adaptive, Morton order)
- [`haar_wavelet_inverse()`](haar_wavelet_inverse.md) : Inverse Haar
  wavelet transform
- [`heat_wavelet_latent()`](heat_wavelet_latent.md) : Heat wavelet
  latent constructor (explicit basis)
- [`heat_wavelet_loadings_handle()`](heat_wavelet_loadings_handle.md) :
  Construct a shared LoadingsHandle via heat-wavelet lifting
- [`hrbf_generate_basis()`](hrbf_generate_basis.md)
  [`hrbf_project_matrix()`](hrbf_generate_basis.md)
  [`hrbf_reconstruct_matrix()`](hrbf_generate_basis.md) : Hierarchical
  radial basis functions (HRBF) for latent fMRI
- [`hrbf_latent()`](hrbf_latent.md) : Build a LatentNeuroVec using an
  HRBF basis
- [`hrbf_meta()`](hrbf_meta.md) : Retrieve HRBF metadata if present
- [`hrbf_reconstruct_partial()`](hrbf_reconstruct_partial.md) :
  Partially reconstruct selected voxels/timepoints
- [`implicit_latent()`](implicit_latent.md) : Construct an
  ImplicitLatent object
- [`implicit_meta()`](implicit_meta.md) : Get metadata from
  ImplicitLatent object
- [`is_explicit_latent()`](is_explicit_latent.md) : Test whether a
  latent object is explicit
- [`is_haar_latent()`](is_haar_latent.md) : Test if object is a Haar
  latent representation
- [`is_hierarchical_template()`](is_hierarchical_template.md) : Check
  whether an object is a HierarchicalBasisTemplate
- [`is_hrbf_latent()`](is_hrbf_latent.md) : Check if latent object
  carries HRBF metadata
- [`is_implicit_latent()`](is_implicit_latent.md) : Test if object is an
  ImplicitLatent
- [`is_shared_reference()`](is_shared_reference.md) : Test whether an
  object is a shared reference
- [`is_surface_template()`](is_surface_template.md) : Test whether an
  object is a surface basis template
- [`is_template()`](is_template.md) : Test whether an object is a
  supported template
- [`is_transport_latent()`](is_transport_latent.md) : Test whether an
  object is a transport-backed latent object
- [`latent_dct_heatwavelet()`](latent_dct_heatwavelet.md) : Create a
  template LatentNeuroVec with heat-wavelet spatial loadings
- [`latent_domain()`](latent_domain.md) : Extract the decoded domain
  associated with a latent object
- [`latent_factory()`](latent_factory.md) : Simple factory to build a
  spec and encode in one call
- [`latent_meta()`](latent_meta.md) : Get lightweight metadata from a
  latent object
- [`latent_searchlight()`](latent_searchlight.md) : Apply a user-defined
  function in latent space over neighborhoods
- [`latent_space_id()`](latent_space_id.md) : Identify latent
  coordinates and decoder domains
- [`latent_support()`](latent_support.md) : Extract the decoded support
  associated with a latent object
- [`latent_units()`](latent_units.md) : Describe Latent Response and
  Coordinate Units
- [`latent_units_record()`](latent_units_record.md) : Declare a Latent
  Units Contract
- [`lift(`*`<ClusterReduction>`*`,`*`<spec_diffusion_wavelet>`*`)`](lift-ClusterReduction-spec_diffusion_wavelet-method.md)
  : Lift diffusion wavelets for clustered reduction
- [`lift(`*`<ClusterReduction>`*`,`*`<spec_heat_wavelet>`*`)`](lift-ClusterReduction-spec_heat_wavelet-method.md)
  : Lift heat wavelets for clustered reduction
- [`lift(`*`<ClusterReduction>`*`,`*`<spec_pca>`*`)`](lift-ClusterReduction-spec_pca-method.md)
  : Lift parcel/cluster-local PCA bases for ClusterReduction
- [`lift(`*`<ClusterReduction>`*`,`*`<spec_slepian>`*`)`](lift-ClusterReduction-spec_slepian-method.md)
  : Lift spatial Slepians for clustered reduction
- [`lift(`*`<GraphReduction>`*`,`*`<ANY>`*`)`](lift-GraphReduction-ANY-method.md)
  : Default lift method (placeholder)
- [`lift()`](lift.md) : Lift reduced bases back to voxel space (abstract
  generic)
- [`linear_access(`*`<LatentNeuroVec>`*`,`*`<numeric>`*`)`](linear_access-methods.md)
  [`linear_access(`*`<LatentNeuroVec>`*`,`*`<integer>`*`)`](linear_access-methods.md)
  : Linear access to LatentNeuroVec elements
- [`list_encoders()`](list_encoders.md) : List registered encoders
- [`lna_hrbf_basis_from_params()`](lna_hrbf_basis_from_params.md) :
  Build an HRBF Basis with Optional Neuroarchive Compatibility Semantics
- [`load_hierarchical_template()`](load_hierarchical_template.md) : Load
  a hierarchical template from disk
- [`load_template()`](load_template.md) : Load a saved template from
  disk
- [`loadings()`](loadings-methods.md) : Get the loadings matrix (spatial
  components)
- [`make_cluster_reduction()`](make_cluster_reduction.md) : Create a
  ClusterReduction from a mask and voxel-to-cluster map
- [`make_coarsened_reduction()`](make_coarsened_reduction.md) : Create a
  coarsened graph reduction
- [`map(`*`<LatentNeuroVec>`*`)`](map-methods.md) : Get the map object
- [`mask(`*`<ImplicitLatent>`*`)`](mask-methods.md)
  [`mask(`*`<LatentNeuroVec>`*`)`](mask-methods.md) : Get the mask
- [`mask_to_array()`](mask_to_array.md) : Convert mask to array
- [`materialize_group_delta_loadings()`](materialize_group_delta_loadings.md)
  : Materialize group-plus-delta loadings
- [`materialize_shared_temporal_spec()`](materialize_shared_temporal_spec.md)
  : Materialize a shared temporal specification
- [`neuroarchive_handoff_contract()`](neuroarchive_handoff_contract.md)
  : Create the fmrilatent-to-neuroarchive handoff contract
- [`offset(`*`<ANY>`*`)`](offset-methods.md)
  [`offset(`*`<LatentNeuroSurfaceVector>`*`)`](offset-methods.md)
  [`offset(`*`<BilatLatentNeuroSurfaceVector>`*`)`](offset-methods.md)
  [`offset(`*`<BlockLatentNeuroVector>`*`)`](offset-methods.md)
  [`offset(`*`<LatentNeuroVec>`*`)`](offset-methods.md) : Get the offset
  vector
- [`parcel_basis_template()`](parcel_basis_template.md) : Build a shared
  parcel basis template
- [`parcel_similarity_matrix()`](parcel_similarity_matrix.md) : Build a
  similarity matrix for parcel clustering (Schaefer-like)
- [`parent_maps_from_levels()`](parent_maps_from_levels.md) : Derive
  parent maps for a nested set of parcellations
- [`plot_basis_gram()`](plot_basis_gram.md) : Plot Gram matrix of a
  basis (orthogonality check)
- [`plot_benchmark_roundtrip()`](plot_benchmark_roundtrip.md) : Plot
  benchmark results
- [`plot_slepian_temporal()`](plot_slepian_temporal.md) : Plot temporal
  Slepians (DPSS)
- [`plot_spatial_atom()`](plot_spatial_atom.md) : Plot a spatial atom
  (loading vector) on a mask
- [`portable_linear_map`](portable_linear_map.md) : Portable linear-map
  contract
- [`predict(`*`<BoldZipSR>`*`)`](predict.BoldZipSR.md) : Predict from a
  BOLDZip-SR codec payload
- [`predict(`*`<HaarLatent>`*`)`](predict.HaarLatent.md) : Predict
  method for HaarLatent
- [`predict(`*`<ImplicitLatent>`*`)`](predict.ImplicitLatent.md) :
  Predict method for ImplicitLatent
- [`print(`*`<ParcelBasisTemplate>`*`)`](print.ParcelBasisTemplate.md) :
  Print method for ParcelBasisTemplate
- [`print(`*`<SurfaceBasisTemplate>`*`)`](print.SurfaceBasisTemplate.md)
  : Print method for SurfaceBasisTemplate
- [`project_effect()`](project_effect.md) : Compatibility wrapper for
  decoder-based coefficient projection
- [`project_hierarchical()`](project_hierarchical.md) : Project
  coefficients only (no LatentNeuroVec wrapper)
- [`project_vcov()`](project_vcov.md) : Compatibility wrapper for
  decoder-based covariance pushforward
- [`reconstruct_array()`](reconstruct_array.md) : Reconstruct a latent
  object as a 4D array
- [`reconstruct_matrix()`](reconstruct_matrix.md) : Reconstruct a latent
  object as a matrix
- [`register_encoder()`](register_encoder.md) : Register an encoder
  family
- [`register_handle_kind()`](register_handle_kind.md) : Register a lazy
  handle materializer
- [`render_shared_events()`](render_shared_events.md) : Render sparse
  events with a shared event dictionary
- [`resolve_shared_reference()`](resolve_shared_reference.md) : Resolve
  an in-session shared reference
- [`roi_subset_columns()`](roi_subset_columns.md) : Subset
  reconstruction matrix columns by ROI mask
- [`save_hierarchical_template()`](save_hierarchical_template.md) : Save
  a hierarchical template to disk
- [`save_template()`](save_template.md) : Save a template object to disk
- [`series(`*`<LatentNeuroVec>`*`,`*`<integer>`*`)`](series-methods.md)
  [`series(`*`<LatentNeuroVec>`*`,`*`<numeric>`*`)`](series-methods.md)
  [`series(`*`<LatentNeuroVec>`*`,`*`<ANY>`*`)`](series-methods.md) :
  Extract time series from LatentNeuroVec
- [`shared_component_contract()`](shared_component_contract.md) : Build
  a method-neutral shared component contract
- [`shared_event_dictionary()`](shared_event_dictionary.md) : Construct
  a reusable event shape dictionary
- [`shared_reference()`](shared_reference.md) : Create an in-session
  shared object reference
- [`shared_reference_clear()`](shared_reference_clear.md) : Reset the
  in-session shared reference cache
- [`shared_reference_info()`](shared_reference_info.md) : Summarize a
  shared reference without resolving it
- [`shared_temporal_spec()`](shared_temporal_spec.md) : Construct a
  shared temporal basis descriptor
- [`show(`*`<LatentNeuroSurfaceVector>`*`)`](show-methods.md)
  [`show(`*`<BilatLatentNeuroSurfaceVector>`*`)`](show-methods.md)
  [`show(`*`<BlockLatentNeuroVector>`*`)`](show-methods.md)
  [`show(`*`<LatentNeuroVec>`*`)`](show-methods.md) : Display a
  LatentNeuroVec object
- [`slepian_spatial_latent()`](slepian_spatial_latent.md) : Slepian
  spatial latent constructor (explicit basis)
- [`slepian_spatial_loadings_handle()`](slepian_spatial_loadings_handle.md)
  : Create a LoadingsHandle for spatial Slepians (graph Laplacian)
- [`slepian_spatiotemporal_latent()`](slepian_spatiotemporal_latent.md)
  : Spatiotemporal Slepian latent (implicit, separable)
- [`slepian_temporal_handle()`](slepian_temporal_handle.md) : Create a
  BasisHandle for temporal Slepians (DPSS)
- [`slepian_temporal_latent()`](slepian_temporal_latent.md) :
  LatentNeuroVec using a temporal DPSS basis
- [`spec_hierarchical_template()`](spec_hierarchical_template.md) :
  Create hierarchical template spec
- [`spec_space_heat()`](spec_space_heat.md) : Spatial heat-wavelet spec
  (graph diffusion)
- [`spec_space_hrbf()`](spec_space_hrbf.md) : Spatial HRBF spec
- [`spec_space_parcel()`](spec_space_parcel.md) : Spatial parcel-basis
  spec (shared/template-based)
- [`spec_space_pca()`](spec_space_pca.md) : Spatial PCA spec
  (cluster-local)
- [`spec_space_slepian()`](spec_space_slepian.md) : Spatial Slepian spec
- [`spec_space_wavelet_active()`](spec_space_wavelet_active.md) :
  Spatial wavelet (active pencil) spec
- [`spec_st()`](spec_st.md) : Spatiotemporal spec (separable)
- [`spec_time_bspline()`](spec_time_bspline.md) : Temporal B-spline spec
- [`spec_time_dct()`](spec_time_dct.md) : Temporal DCT spec
- [`spec_time_slepian()`](spec_time_slepian.md) : Temporal Slepian/DPSS
  spec
- [`spectral_ward_hclust()`](spectral_ward_hclust.md) : Run
  spectral+Ward hierarchical clustering on a parcel graph
- [`surface_basis_template()`](surface_basis_template.md) : Build a
  shared surface basis template
- [`template_domain()`](template_domain.md) : Extract the domain
  associated with a template
- [`template_loadings()`](template_loadings.md) : Extract template
  loadings
- [`template_mask()`](template_mask.md) : Extract template mask
- [`template_measure()`](template_measure.md) : Extract optional measure
  or mass information for a template
- [`template_meta()`](template_meta.md) : Extract template metadata
- [`template_project()`](template_project.md) : Project data onto a
  template
- [`template_rank()`](template_rank.md) : Query the rank of a template
  basis
- [`template_roughness()`](template_roughness.md) : Extract the spatial
  roughness operator for a template asset
- [`template_support()`](template_support.md) : Extract the support
  associated with a template
- [`transport_latent()`](transport_latent.md) : Construct a
  transport-backed implicit latent object
- [`validate_nested_parcellations()`](validate_nested_parcellations.md)
  : Validate that parcellation levels are nested
- [`validate_neuroarchive_handoff_contract()`](validate_neuroarchive_handoff_contract.md)
  : Validate the fmrilatent-to-neuroarchive handoff contract
- [`validate_portable_linear_map()`](validate_portable_linear_map.md) :
  Validate an object against the portable linear-map contract
- [`validate_shared_component_contract()`](validate_shared_component_contract.md)
  : Validate a shared component contract
- [`validate_template_protocol()`](validate_template_protocol.md) :
  Validate the reusable template protocol
- [`voxel_subset_to_gsp()`](voxel_subset_to_gsp.md) : Convert voxel
  subset to an rgsp graph
- [`wavelet_active_latent()`](wavelet_active_latent.md) : Active-pencil
  wavelet latent (CDF 5/3)
- [`wrap_decoded()`](wrap_decoded.md) : Wrap flat decoded outputs into a
  domain-native representation
