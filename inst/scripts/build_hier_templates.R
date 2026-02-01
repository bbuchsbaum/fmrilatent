# Offline build script (not run on CRAN)
#
# Goal: generate reusable HierarchicalBasisTemplate objects for MNI space.
# Inputs: MNI gray-matter mask/prob map, Schaefer-400 cortical labels (with Yeo17 tags),
#         subcortical atlas (e.g., FreeSurfer aseg or Harvard-Oxford subcortical).
# Outputs: .rds templates under inst/extdata (small/medium/large SKUs) plus a manifest.

# Pseudocode sketch (fill in with actual file paths and atlas readers in your environment):

# library(fmrilatent)
# library(neuroim2)
# library(Matrix)

# load_mask <- function(path) {
#   vol <- neuroim2::read_vol(path)
#   LogicalNeuroVol(vol@.Data > 0.2, neuroim2::space(vol))
# }

# build_levels <- function(mask, schaefer400, aseg) {
#   # Align to mask voxel order; returns list of integer label vectors (mask order)
#   # Level3 (fine): Schaefer-400 cortex + split subcortical nuclei
#   # Level2: aggregate fine cortex to ~64 parcels (spectral or Ward on centroids), subcortex merged to nuclei
#   # Level1: aggregate to ~32 networks using Yeo17 mapping + hemispheres
#   # Level0: whole brain = 1
# }

# template_one <- function(mask, levels, k_per_level, sku, ridge = 1e-8) {
#   tmpl <- build_hierarchical_template(mask = mask, parcellations = levels,
#                                       k_per_level = k_per_level, ridge = ridge,
#                                       label = paste0("hier_", sku))
#   saveRDS(tmpl, file = file.path("inst", "extdata", paste0("hier_", sku, ".rds")))
#   invisible(tmpl)
# }

# main <- function() {
#   mask <- load_mask("/path/to/mni_gm_prob.nii.gz")
#   schaefer400 <- read_schaefer_labels(mask)  # user-defined helper
#   aseg <- read_subcortex_labels(mask)        # user-defined helper
#   levels <- build_levels(mask, schaefer400, aseg)

#   # Size presets
#   k_small  <- c(8, 3, 1)      # levels 0/1/2 (stop at ~32 parcels)
#   k_medium <- c(8, 5, 3, 1)   # levels 0/1/2/3 (up to ~64/400)
#   k_large  <- c(8, 5, 3, 1)   # same structure; adjust if adding an 800-tile leaf

#   template_one(mask, levels[1:3], k_small,  sku = "mni_small")
#   template_one(mask, levels[1:4], k_medium, sku = "mni_medium")
#   # template_one(mask, levels[1:5], k_large, sku = "mni_large")  # optional finer leaf
# }

# main()

# Notes
# - Keep voxel order consistent across mask, labels, and saved templates; record it in manifest.
# - Store a manifest (CSV/JSON) listing atoms with level/parcel/mode/parent for user interpretability.
# - Subcortex: allow nuclei splits for large structures; keep k_per_level low (1–3) to avoid unstable modes.
# - Spectral aggregation step should enforce nesting (each child belongs to exactly one parent).
# - Run offline; committed outputs are small (~tens of MB) once saved as sparse matrices.
