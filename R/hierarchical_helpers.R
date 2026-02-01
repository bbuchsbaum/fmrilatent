# Hierarchical parcellation helpers (neuroimaging-agnostic)

# -----------------------------------------------------------------------------
# Schaefer surface-based hierarchical construction
# -----------------------------------------------------------------------------

#' Build hierarchical parcellation levels from Schaefer surface atlas
#'
#' @description
#' Constructs nested parcellation levels for the hierarchical template by:
#' 1. Loading Schaefer surface atlas from neuroatlas
#' 2. Computing geodesic distances and boundary contacts via neurosurf
#' 3. Building a parcel similarity matrix (geodesic + boundary + network)
#' 4. Performing spectral Ward clustering to create coarser levels
#' 5. Mapping surface parcels back to volumetric voxel labels
#'
#' Real-data glue to finish:
#' - Provide fsaverage Schaefer surfaces (via `neuroatlas::schaefer_surf`) so we can compute geodesic/boundary terms.
#' - Provide the volumetric Schaefer atlas aligned to the MNI mask (via `neuroatlas::get_schaefer_atlas`).
#' - Ensure mask and volumetric atlas share resolution/origin/spacing; this function assumes 2 mm MNI unless you pass a custom `vol_atlas`.
#' - Subcortex is **not** handled here; add it separately when assembling the full template.
#' 
#' @param mask LogicalNeuroVol defining the 3D brain mask (MNI space).
#' @param parcels Integer. Number of Schaefer parcels (100, 200, 300, 400, 500, 600, 800, 1000).
#' @param networks Integer. Yeo network variant (7 or 17).
#' @param space Character. Surface space ("fsaverage6", "fsaverage", "fsaverage5").
#' @param k_levels Integer vector of cluster counts for coarser levels, ordered coarse→fine.
#'   The finest level is always the original Schaefer parcellation. Example: c(16, 32, 64)
#'   produces 4 levels: L0=16, L1=32, L2=64, L3=Schaefer-parcels.
#' @param vol_atlas Optional. Pre-loaded volumetric Schaefer atlas (NeuroVol with parcel IDs).
#'   If NULL, attempts to load via neuroatlas::get_schaefer_atlas().
#' @param alpha Weight for boundary contact in similarity (default 0.5).
#' @param beta Weight for geodesic kernel in similarity (default 0.3).
#' @param gamma Weight for Yeo network agreement in similarity (default 0.2).
#' @param d0 Scale for geodesic exponential kernel (default 30 mm).
#' @param component_policy How to handle fragmented parcels ("error", "largest", "each", "merge").
#' @param cache Logical. Use cached geodesic distances if available (default TRUE).
#'
#' @return A list with components:
#'   \item{levels}{List of integer vectors (length = #voxels in mask), one per level, nested.}
#'   \item{network}{Character vector of Yeo network labels per finest-level parcel.}
#'   \item{hemi}{Character vector of hemisphere ("L"/"R") per finest-level parcel.}
#'   \item{geo_dist_lh, geo_dist_rh}{Parcel geodesic distance matrices (for diagnostics).}
#'   \item{boundary_lh, boundary_rh}{Boundary contact matrices (LH/RH).}
#' 
#' @details
#' Requires `neuroatlas` and `neurosurf` packages. The surface atlas provides geodesic
#' distances and boundary topology; the volumetric atlas maps parcels to mask voxels.
#'
#' @export
build_schaefer_levels <- function(mask,
                                   parcels = 400,
                                   networks = 17,
                                   space = "fsaverage6",
                                   k_levels = c(16, 32, 64),
                                   vol_atlas = NULL,
                                   alpha = 0.5,
                                   beta = 0.3,
                                   gamma = 0.2,
                                   d0 = 30,
                                   component_policy = c("largest", "error", "each", "merge"),
                                   cache = TRUE) {


  component_policy <- match.arg(component_policy)


  # --- Check dependencies ---

  if (!requireNamespace("neuroatlas", quietly = TRUE)) {
    stop("neuroatlas package required; install from bbuchsbaum/neuroatlas", call. = FALSE)
  }
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop("neurosurf package required; install from bbuchsbaum/neurosurf", call. = FALSE)
  }


  # --- Load surface atlas (LabeledNeuroSurface objects) ---
  surf_atlas <- neuroatlas::schaefer_surf(
    parcels = parcels,
    networks = networks,
    space = space,
    surf = "white"
  )


  # --- Load volumetric atlas if not provided ---
  if (is.null(vol_atlas)) {
    vol_atlas <- neuroatlas::get_schaefer_atlas(
      parcels = parcels,
      networks = networks,
      resolution = "2mm"
    )
  }


  # --- Compute geodesic distance matrices per hemisphere ---
  geo_dist_lh <- neurosurf::parcel_geodesic_distance_matrix(
    surf_atlas$lh_atlas,
    metric = "centroid",
    component_policy = component_policy,
    cache = cache
  )
  geo_dist_rh <- neurosurf::parcel_geodesic_distance_matrix(
    surf_atlas$rh_atlas,
    metric = "centroid",
    component_policy = component_policy,
    cache = cache
  )

  # --- Compute boundary contact matrices per hemisphere ---
  boundary_lh <- neurosurf::parcel_boundary_contact(
    surf_atlas$lh_atlas,
    component_policy = component_policy,
    counts = FALSE
  )
  boundary_rh <- neurosurf::parcel_boundary_contact(
    surf_atlas$rh_atlas,
    component_policy = component_policy,
    counts = FALSE
  )

  # --- Assemble full parcel similarity matrices per hemisphere ---
  # Get network labels per parcel (LH parcels are 1:n_lh, RH are n_lh+1:n_total)
  n_lh <- nrow(geo_dist_lh)
  n_rh <- nrow(geo_dist_rh)
  n_parc <- n_lh + n_rh

  # Extract network assignments from atlas metadata
  network_labels <- surf_atlas$network
  hemi_labels <- surf_atlas$hemi

  yeo_lh <- network_labels[seq_len(n_lh)]
  yeo_rh <- network_labels[seq(n_lh + 1, n_parc)]

  W_lh <- parcel_similarity_matrix(
    boundary_contact = boundary_lh * 1.0,  # convert logical to numeric
    geo_dist = geo_dist_lh,
    yeo17 = yeo_lh,
    alpha = alpha, beta = beta, gamma = gamma, d0 = d0
  )
  W_rh <- parcel_similarity_matrix(
    boundary_contact = boundary_rh * 1.0,
    geo_dist = geo_dist_rh,
    yeo17 = yeo_rh,
    alpha = alpha, beta = beta, gamma = gamma, d0 = d0
  )

  # --- Hierarchical clustering per hemisphere ---
  hc_lh <- spectral_ward_hclust(W_lh, k_embed = 3, hemi = NULL, network = yeo_lh)
  hc_rh <- spectral_ward_hclust(W_rh, k_embed = 3, hemi = NULL, network = yeo_rh)

  # --- Cut into nested levels ---
  # k_levels are for each hemisphere, so divide by 2 (roughly symmetric)
  k_per_hemi <- pmax(1L, as.integer(k_levels / 2))

  levels_lh <- cut_hclust_nested(hc_lh, k_per_hemi)
  levels_rh <- cut_hclust_nested(hc_rh, k_per_hemi)

  # --- Map parcel labels to voxel space ---
  mask_arr <- as.array(mask)
  mask_idx <- which(mask_arr)
  n_vox <- length(mask_idx)

  # Get volumetric parcel labels at mask voxels
  vol_labels <- vol_atlas$atlas[mask_idx]

  # Build voxel-level parcellations for each hierarchy level
  # Finest level = original Schaefer labels
  vox_levels <- vector("list", length(k_levels) + 1L)

  # For coarser levels, map: voxel → parcel → cluster

  for (lvl in seq_along(k_levels)) {
    vox_labels <- integer(n_vox)

    for (i in seq_len(n_vox)) {
      parc_id <- vol_labels[i]
      if (is.na(parc_id) || parc_id == 0L) {
        vox_labels[i] <- 0L
        next
      }

      # Determine hemisphere and local parcel index
      if (parc_id <= n_lh) {
        # Left hemisphere
        cluster_id <- levels_lh[[lvl]][parc_id]
        vox_labels[i] <- cluster_id
      } else {
        # Right hemisphere - offset cluster IDs to avoid collision
        local_id <- parc_id - n_lh
        cluster_id <- levels_rh[[lvl]][local_id]
        vox_labels[i] <- cluster_id + max(levels_lh[[lvl]])
      }
    }
    vox_levels[[lvl]] <- vox_labels
  }

  # Finest level = original parcel IDs
  vox_levels[[length(k_levels) + 1L]] <- as.integer(vol_labels)

  # Validate nesting
  validate_nested_parcellations(vox_levels)

  list(
    levels = vox_levels,
    network = network_labels,
    hemi = hemi_labels,
    geo_dist_lh = geo_dist_lh,
    geo_dist_rh = geo_dist_rh,
    boundary_lh = boundary_lh,
    boundary_rh = boundary_rh,
    k_levels = c(k_levels, parcels)
  )
}


#' Build hierarchical template from Schaefer surface atlas
#'
#' @description
#' Convenience wrapper that combines `build_schaefer_levels()` with
#' `build_hierarchical_template()` to produce a ready-to-use template
#' with geodesic-informed parcel clustering.
#'
#' @inheritParams build_schaefer_levels
#' @param k_per_level Integer vector of eigenmodes per parcel at each level.
#'   Length must match length(k_levels) + 1 (for the finest Schaefer level).
#' @param k_neighbors k for local graph construction inside parcels.
#' @param ridge Small diagonal ridge for Gram matrix stability.
#' @param solver Solver choice: "chol" or "qr".
#'
#' @return HierarchicalBasisTemplate object.
#'
#' @seealso \code{\link{build_hierarchical_template}}, \code{\link{build_schaefer_levels}}
#' @export
build_schaefer_hierarchical_template <- function(mask,
                                                  parcels = 400,
                                                  networks = 17,
                                                  space = "fsaverage6",
                                                  k_levels = c(16, 32, 64),
                                                  k_per_level = c(8, 5, 3, 1),
                                                  vol_atlas = NULL,
                                                  alpha = 0.5,
                                                  beta = 0.3,
                                                  gamma = 0.2,
                                                  d0 = 30,
                                                  component_policy = "largest",
                                                  cache = TRUE,
                                                  k_neighbors = 6L,
                                                  ridge = 1e-8,
                                                  solver = c("chol", "qr")) {

  solver <- match.arg(solver)

  # Build nested parcellation levels
  levels_info <- build_schaefer_levels(
    mask = mask,
    parcels = parcels,
    networks = networks,
    space = space,
    k_levels = k_levels,
    vol_atlas = vol_atlas,
    alpha = alpha, beta = beta, gamma = gamma, d0 = d0,
    component_policy = component_policy,
    cache = cache
  )

  if (length(k_per_level) != length(levels_info$levels)) {
    stop("k_per_level length (", length(k_per_level), ") must match number of levels (",
         length(levels_info$levels), ")", call. = FALSE)
  }

  # Build the template
  build_hierarchical_template(
    mask = mask,
    parcellations = levels_info$levels,
    k_per_level = k_per_level,
    k_neighbors = k_neighbors,
    ridge = ridge,
    solver = solver,
    label = paste0("schaefer", parcels, "_", networks, "net_geodesic")
  )
}


# -----------------------------------------------------------------------------
# Generic hierarchical helpers (neuroimaging-agnostic)
# -----------------------------------------------------------------------------

#' Cut an hclust into nested label vectors
#'
#' @description
#' Utility to turn a hierarchical clustering (on parcels) into a nested set of
#' parcel label vectors that downstream code can lift into voxel space. Real
#' data glue: you will pass the hclust built on a parcel-level similarity graph
#' (e.g., Schaefer-400 parcel graph constructed from surface geodesics and
#' boundary contacts); then map these parcel labels back to voxels using the
#' volumetric atlas.
#'
#' @param hc hclust object (e.g., from spectral_ward_hclust()).
#' @param k_levels Integer vector of cluster counts, ordered coarse→fine (e.g., c(16, 32, 64, 400)).
#' @return List of integer label vectors (same length as hc$labels), one per level, nested by construction.
#' @export
cut_hclust_nested <- function(hc, k_levels) {
  if (!inherits(hc, "hclust")) stop("hc must be an hclust object")
  if (is.unsorted(k_levels, strictly = FALSE)) {
    warning("k_levels not sorted; coercing to increasing")
    k_levels <- sort(unique(as.integer(k_levels)))
  }
  k_levels <- as.integer(k_levels)
  if (any(k_levels < 1L)) stop("k_levels must be >= 1")

  levels <- lapply(k_levels, function(k) stats::cutree(hc, k = k))
  validate_nested_parcellations(levels)
  levels
}

#' Build a similarity matrix for parcel clustering (Schaefer-like)
#'
#' @description
#' Combines geometric and functional priors into a parcel–parcel similarity
#' matrix W used to cluster Schaefer parcels into coarser levels. Real data glue:
#' - `boundary_contact`: from surface meshes (fsaverage), compute fraction of shared boundary between parcels.
#' - `geo_dist`: geodesic distance between parcel centroids on the surface (fallback: Euclidean in MNI).
#' - `yeo17`: network labels from the Schaefer/Yeo annotation (hemisphere-specific).
#'
#' @param boundary_contact Numeric matrix [n_parc x n_parc], symmetric, contact fraction (0..1).
#' @param geo_dist Numeric matrix [n_parc x n_parc], symmetric geodesic distances (mm).
#' @param yeo17 Integer or factor vector of length n_parc with Yeo17 network labels.
#' @param alpha Weight for boundary contact (default 0.5).
#' @param beta Weight for geodesic kernel (default 0.3).
#' @param gamma Weight for Yeo17 agreement (default 0.2).
#' @param d0 Scale for geodesic exponential (default 30 mm).
#' @return Similarity matrix W (symmetric, zero diagonal).
#' @export
parcel_similarity_matrix <- function(boundary_contact, geo_dist, yeo17,
                                     alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30) {
  if (!is.matrix(boundary_contact) || !is.matrix(geo_dist)) stop("boundary_contact and geo_dist must be matrices")
  if (!all(dim(boundary_contact) == dim(geo_dist))) stop("boundary_contact and geo_dist must have same dimensions")
  if (nrow(boundary_contact) != length(yeo17)) stop("yeo17 length must match parcel count")
  n <- nrow(boundary_contact)
  yeo17 <- as.integer(as.factor(yeo17))
  yeo_mat <- outer(yeo17, yeo17, FUN = "==") * 1.0
  geo_k <- exp(-geo_dist / d0)
  W <- alpha * boundary_contact + beta * geo_k + gamma * yeo_mat
  diag(W) <- 0
  # enforce symmetry
  W <- (W + t(W)) / 2
  W
}

#' Run spectral+Ward hierarchical clustering on a parcel graph
#'
#' @description
#' Produces an hclust over parcels using a spectral embedding of the similarity
#' matrix, then Ward linkage. Real data glue: supply W from `parcel_similarity_matrix()`,
#' `hemi` from parcel hemisphere labels, and `network` from Yeo17 labels. Penalties
#' bias merges to stay within hemisphere/network until higher in the tree.
#'
#' @param W Similarity matrix (symmetric).
#' @param k_embed Embedding dimension (2 or 3).
#' @param hemi Optional factor/character vector (n_parc) to enforce within-hemisphere merges first.
#' @param network Optional factor/character vector (n_parc) (e.g., Yeo17) to enforce within-network merges until relaxed.
#' @return hclust object.
#' @export
spectral_ward_hclust <- function(W, k_embed = 3, hemi = NULL, network = NULL) {
  if (!is.matrix(W)) stop("W must be a matrix")
  n <- nrow(W)
  if (n != ncol(W)) stop("W must be square")

  # Degree & Laplacian
  d <- rowSums(W)
  D_inv_sqrt <- 1 / sqrt(pmax(d, .Machine$double.eps))
  L <- diag(n) - (D_inv_sqrt * W) * D_inv_sqrt

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("RSpectra is required for spectral embedding", call. = FALSE)
  }

  k_use <- min(k_embed, n - 1L)
  eig <- RSpectra::eigs(L, k = k_use, which = "SM")
  embed <- eig$vectors

  # Optionally enforce hemisphere/network by scaling distances
  if (!is.null(hemi)) {
    hemi <- as.integer(as.factor(hemi))
    penalty <- outer(hemi, hemi, FUN = "!=") * 10  # big penalty across hemis
  } else {
    penalty <- matrix(0, n, n)
  }
  if (!is.null(network)) {
    network <- as.integer(as.factor(network))
    penalty <- penalty + outer(network, network, FUN = "!=") * 2  # softer penalty across networks
  }

  dist_mat <- stats::dist(embed)
  dist_mat <- as.matrix(dist_mat) + penalty
  hc <- stats::hclust(stats::as.dist(dist_mat), method = "ward.D2")
  hc
}

#' Derive parent maps for a nested set of parcellations
#'
#' @param levels List of integer label vectors (all same length), coarse→fine or fine→coarse; order agnostic.
#' @return List of integer named vectors: parents[[lvl]] maps child ids (names) at level lvl to parent ids at lvl-1; parents[[1]] is integer(0).
#' @export
parent_maps_from_levels <- function(levels) {
  validate_nested_parcellations(levels)
  # ensure coarse→fine order for parent mapping
  parents <- vector("list", length(levels))
  parents[[1L]] <- integer(0)
  if (length(levels) > 1L) {
    for (i in 2:length(levels)) {
      parents[[i]] <- .parent_map_one(child = levels[[i]], parent = levels[[i - 1L]], lvl = i)
    }
  }
  parents
}

#' Validate that parcellation levels are nested
#'
#' @param levels List of integer label vectors (same length) representing parcellations.
#' @return Invisibly TRUE; stops with an error if nesting fails or lengths differ.
#' @export
validate_nested_parcellations <- function(levels) {
  if (!is.list(levels) || length(levels) == 0L) stop("levels must be a non-empty list")
  lens <- vapply(levels, length, integer(1))
  if (length(unique(lens)) != 1L) stop("all parcellation levels must have the same length")
  if (length(levels) == 1L) return(invisible(TRUE))
  for (i in 2:length(levels)) {
    .parent_map_one(child = levels[[i]], parent = levels[[i - 1L]], lvl = i) # will stop on mismatch
  }
  invisible(TRUE)
}

# internal helper
.parent_map_one <- function(child, parent, lvl) {
  child_ids <- sort(unique(child))
  parent_map <- integer(length(child_ids))
  names(parent_map) <- as.character(child_ids)
  for (cid in child_ids) {
    idx <- which(child == cid)
    if (length(idx) == 0L) next
    p_ids <- unique(parent[idx])
    p_ids <- p_ids[!is.na(p_ids)]
    if (length(p_ids) != 1L) {
      stop("Levels are not nested at level ", lvl, ": child ", cid, " maps to ", paste(p_ids, collapse = ","))
    }
    parent_map[as.character(cid)] <- p_ids
  }
  parent_map
}
