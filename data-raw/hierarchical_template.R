# Build a hierarchical Laplacian template (offline) and stash it under inst/extdata.
# Minimal knobs are hardcoded below; adjust as needed. Later we can loop over
# (space, parcels, networks) combinations for multiple SKUs.

suppressPackageStartupMessages(library(fmrilatent))
suppressPackageStartupMessages(library(neuroim2))
suppressPackageStartupMessages(library(neuroatlas))
suppressPackageStartupMessages(library(neurosurf))

# --- Constants ---------------------------------------------------------------

space <- "MNI152NLin2009cSym"   # TemplateFlow space for mask/atlas alignment
resolution <- 2                 # mm; must match Schaefer atlas resolution
mask_thresh <- 0.2              # probability threshold for GM mask
parcels <- 400L                 # Schaefer parcel count
networks <- 17L                 # Yeo network variant
k_levels <- c(16L, 32L, 64L)    # coarse -> fine clustering targets
k_per_level <- c(8L, 5L, 3L, 1L) # eigenmodes per parcel per level (includes finest level)
k_neighbors <- 6L
ridge <- 1e-8
solver <- "chol"
out_dir <- file.path("inst", "extdata")
out_file <- sprintf("hier_%s_schaefer%d_%dnet.rds", space, parcels, networks)
# TODO: add subcortical labels (aseg/Harvard-Oxford) and merge before building.

# --- Load mask/atlas ---------------------------------------------------------

# Gray-matter probability map in the chosen space; falls back to get_template()
# because get_template_gm is deprecated.
mask_img <- neuroatlas::get_template(
  space = space,
  variant = "probseg",
  label = "GM",
  resolution = resolution,
  extension = ".nii.gz"
)
mask_vol <- LogicalNeuroVol(mask_img@.Data > mask_thresh, neuroim2::space(mask_img))

# Schaefer volumetric atlas resampled into the mask space.
vol_atlas <- neuroatlas::get_schaefer_atlas(
  parcels = parcels,
  networks = networks,
  resolution = resolution,
  outspace = neuroim2::space(mask_vol),
  use_cache = TRUE
)$atlas

# --- Build template (Schaefer geodesic path) ---------------------------------

tmpl <- build_schaefer_hierarchical_template(
  mask = mask_vol,
  parcels = parcels,
  networks = networks,
  space = "fsaverage6",
  k_levels = k_levels,
  k_per_level = k_per_level,
  vol_atlas = vol_atlas,
  alpha = 0.5,
  beta = 0.3,
  gamma = 0.2,
  d0 = 30,
  component_policy = "largest",
  cache = TRUE,
  k_neighbors = k_neighbors,
  ridge = ridge,
  solver = solver
)

# --- Save outputs ------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, out_file)
save_hierarchical_template(tmpl, file = out_path, compress = "xz")
utils::write.csv(tmpl@atoms, file = file.path(out_dir, sub("\\.rds$", "_atoms.csv", out_file)), row.names = FALSE)

message("Saved hierarchical template to: ", normalizePath(out_path))
