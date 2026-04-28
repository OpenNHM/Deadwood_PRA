# Title: Assessing the Protective Effect of Lying Deadwood Against Snow Avalanches
#
#Please set input parameters to work
#
# Description:
# This script computes avalanche release membership and assesses the protective
# effect of lying deadwood dominated forest areas against snow avalanche release.
# It is designed for disturbed mountain forests where lying deadwood strongly
# influences surface roughness and therefore the likelihood of avalanche
# release. The workflow combines terrain information, deadwood
# structure, and canopy coverage to derive spatially explicit fuzzy membership maps.
#
# Authors: Leon Buehrle, Tommaso Baggio
# Contact: leon.buehrle@t-online.de; tbaggio93@gmail.com
# Year: 2026
#
# Notes:
# - Developed for high-resolution UAV and ULS deadwood structure analysis.
# - Code was sucessfully tested across 5 study sites including windthrow and snagfall dominated areas
# - Code was tested using different dense point clouds derived from low-cost DJI UAV systems to advanced DJI ULS Zenmuse L2 systems
# - The workflow integrates terrain, deadwood structure, canopy cover, and surface roughness.
#
# Disclaimer:
# THE SCRIPT IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ==================================================================================================
# Load required libraries
# ==================================================================================================

# LiDAR processing (ideally requires version 4.0.2)
library(lidR) # if higher version is used, do not change threshold_removevegetation > 3.9 (bug)

# Modern raster operations
library(terra)

# Spatial vector data handling
library(sf)

# Data manipulation
library(dplyr)

# Unit handling
library(units)

# Fast rasterization
library(fasterize)

# Additional spatial tools
library(spatialEco)

# Parallel computing
library(parallel)
library(doParallel)
library(foreach)

# Progress-bar apply functions
library(pbapply)

# Graph-based crown merging
library(igraph)

# ==================================================================================================
# Folder structure
# ==================================================================================================

# Set working directory containing the project data
# -> define your main project folder here
project_path <- "path/to/your/project"
setwd(project_path)

# Define output directory
# Define output directory (results will be stored here)
results_path <- file.path(project_path, "results")

# ==================================================================================================
# USER INPUT SETTINGS
# ==================================================================================================

# Configure key processing options before running the script.
# These settings control optional filtering steps and should be adjusted before
# launching the workflow. They are intended as the main user entry point.


# Define and create result subfolders
subfolders <- c(
  "zones", "VHM", "roughness", "slope_membership", "roughness_membership", "adapted_tree_parameters",
  "VHM_filled", "fuzzy_logic"
)

for (folder in subfolders) {
  full_path <- file.path(results_path, folder)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    message(paste("Subfolder", folder, "created."))
  } else {
    message(paste("Subfolder", folder, "already exists."))
  }
}

# ==================================================================================================
# Input data and parameter settings
# ==================================================================================================

# Load dense point cloud (X, Y, Z, R, G, B, n, r) and assign CRS
# The input point cloud should represent the snow-free surface and should
# contain RGB values if RGBVI filtering is intended to be used.
#nr is only used for LiDAR data
# Load LiDAR point cloud (LAS/LAZ file)
las_path <- file.path(project_path, "data", "pointcloud.laz")
las <- readLAS(las_path, select = "xyzRGBnr")
# Define CRS of point cloud (must match DTM)
st_crs(las) <- "EPSG:XXXX" 

# Load digital terrain model and assign CRS
# The DTM should describe the terrain beneath vegetation and deadwood and must
# match  coordinate reference system.
# Load Digital Terrain Model (GeoTIFF)
dtm_path <- file.path(project_path, "data", "DTM.tif")
DTM <- rast(dtm_path)
# Define coordinate reference system (use EPSG code or proj string)
crs(DTM) <- "EPSG:XXXX"

# Spatial resolution settings
# res_fin controls the base processing resolution for detailed raster products.
# res_fuzzy_logic defines the coarser resolution used for fuzzy aggregation.
# res_redistribution defines the zone size used for snow redistribution analysis.
res_fuzzy_logic <- 5      # (default: 5)
res_fin <- 0.2            # (default: 0.2)
res_redistribution <- 20  # (default: 20)

# Apply RGBVI-based vegetation filtering.
# TRUE  = additional spectral filtering after CSF classification
# FALSE = only CSF-based filtering is used
use_rgbvi_filter <- TRUE  # (default: TRUE)

# Threshold for RGBVI filtering.
# Adjust depending on sensor characteristics and illumination conditions.
# More negative values keep more points, while less negative or higher values
# remove a larger proportion of vegetation-like points.
rgbvi_threshold <- -0.015  # (default: -0.015)

# Height threshold used in CSF-based vegetation removal
threshold_removevegetation <- 3.5  # (default: 3.5) # if lidr version >4.0.2 is used, do not change threshold_removevegetation > 3.9 (bug)

# Quantile-based upfilling settings
# These parameters define the range and increment of quantiles used to estimate
# representative deadwood height within redistribution zones.
upfilling_quantile_step <- 0.003   # (default: 0.003)
upfilling_quantile_start <- 0.4    # (default: 0.4)
upfilling_quantile_end <- 1        # (default: 1)

# Snow-depth settings for roughness and fuzzy logic computation
# These parameters define the modeled snow depth range used to evaluate how the
# roughness effect of lying deadwood decreases as snow cover increases.
upfilling_SVH_start <- 0.0  # (default: 0.0)
upfilling_SVH_end <- 4      # (default: 4)
upfilling_SVH_step <- 0.1   # (default: 0.1)

# Tree detection and crown delineation parameters
# These parameters affect tree-top detection and crown delineation and therefore
# indirectly influence canopy-related membership estimates.
hmin_tree <- 5                    # (default: 5)
cmin_tree <- 2                    # (default: 2)
crown_min_height <- 4             # (default: 4)
crown_min_height_ratio <- 0.3     # (default: 0.3)
max_cr_factor <- 0.35             # (default: 0.35)

# Computational power
# Define number of cores (consider ~5 GB RAM per CPU)
num_cores <- 60  # (default: 60) cores (consider ~5 GB RAM per CPU)

# ==================================================================================================
# Start of the workflow
# ==================================================================================================

# Inspect input point cloud
print(las)

# Create digital surface model (DSM)
DSM <- rast(grid_metrics(las, res = res_fin, ~mean(Z)))

# Resample terrain model to the DSM grid
DTM_res <- resample(DTM, DSM, method = "bilinear")

# Normalize point cloud using terrain elevation

las_normalize <- normalize_height(las, DTM_res)

# Generate vegetation height model (VHM) including crowns
VHM_las <- rast(grid_metrics(las_normalize, res = res_fin, ~quantile(Z, 0.995)))

# Fill gaps in the VHM using median filters
# Median filtering reduces small gaps and noise in the canopy surface and helps
# produce more stable crown delineation results.
VHM_filled <- focal(VHM_las, w = matrix(1, 3, 3), fun = median, na.rm = TRUE, fillvalue = NA)
VHM_crowns_filled <- cover(
  VHM_filled,
  focal(VHM_filled, w = matrix(1, 5, 5), fun = median, na.rm = TRUE, fillvalue = NA)
)

# Free memory
gc()

# ==================================================================================================
# Tree-top detection and crown delineation
# ==================================================================================================

# Detect tree tops using a local maximum filter.
# Local maxima are interpreted as candidate tree tops. In stands with many tall
# trees, the lmf window size may need adjustment to avoid over- or under-detection.
ttops <- find_trees(VHM_crowns_filled, lmf(5, hmin = hmin_tree, shape = "circular"))

# Delineate crowns using the Silva et al. (2016) approach
# This crown segmentation step is used to derive canopy structure and to later
# estimate canopy coverage above deadwood-dominated surfaces.
crowns <- silva2016(
  chm = VHM_crowns_filled,
  treetops = ttops,
  exclusion = crown_min_height_ratio,
  max_cr_factor = max_cr_factor
)()

# Remove crown areas below the minimum crown height
crowns[VHM_crowns_filled <= crown_min_height] <- NA

# Convert crowns to polygons and calculate crown area
crown_polygons <- as.polygons(rast(crowns), dissolve = TRUE)
crown_polygons$area <- expanse(crown_polygons)

# Filter small crown polygons
tree_crown_shape <- crown_polygons[crown_polygons$area > cmin_tree, ]

# Convert crowns to sf and split multipart geometries
tree_sf <- st_as_sf(tree_crown_shape)

split_polys <- lapply(seq_len(nrow(tree_sf)), function(i) {
  geom <- tree_sf[i, ]
  st_cast(geom, "POLYGON")
})

tree_sf_single <- do.call(rbind, split_polys)

# Optional conversion back to terra
tree_polygons <- vect(tree_sf_single)

# Inspect output geometry types
tree_sf_single
st_geometry_type(tree_sf_single)

# Assign crown IDs to connected polygons
# Crown fragments belonging to the same segmented crown are grouped using graph-
# based connectivity, allowing nearby touching polygons to be merged.
tree_sf_single <- tree_sf_single %>%
  mutate(crown_id = NA_integer_)

ayers <- unique(tree_sf_single$focal_median)

for (l in layers) {
  sub <- tree_sf_single %>% filter(focal_median == l)
  
  if (nrow(sub) == 1) {
    tree_sf_single$crown_id[tree_sf_single$focal_median == l] <- l
    next
  }
  
  dmat <- st_distance(sub)
  dmat_num <- matrix(as.numeric(dmat), nrow = nrow(sub))
  g <- graph_from_adjacency_matrix(dmat_num <= 1, mode = "undirected", diag = FALSE)
  comps <- components(g)$membership
  
  current_max_id <- suppressWarnings(max(tree_sf_single$crown_id, na.rm = TRUE))
  if (is.infinite(current_max_id)) current_max_id <- 0
  
  new_ids <- current_max_id + comps
  tree_sf_single$crown_id[tree_sf_single$focal_median == l] <- new_ids
}

crown_merged <- tree_sf_single %>%
  group_by(crown_id) %>%
  summarise(
    focal_median = first(focal_median),
    geometry = st_union(geometry),
    .groups = "drop"
  )


# Recalculate crown area and remove very small crowns
crown_merged <- crown_merged %>%
  mutate(area = st_area(geometry)) %>%
  mutate(area = as.numeric(area)) %>%
  filter(area >= 2)

# Keep only tree tops located inside crown polygons
ttops <- st_as_sf(ttops)
ttops <- st_transform(ttops, st_crs(crown_merged))
inside <- st_within(ttops, crown_merged, sparse = TRUE)
ttops_filtered <- ttops[lengths(inside) > 0, ]

# Rasterize crown polygons
crown_raster <- rasterize(vect(crown_merged), VHM_crowns_filled, field = 1)

# ==================================================================================================
# Deadwood classification and vegetation filtering
# ==================================================================================================

# Classify deadwood points using the Cloth Simulation Filter (CSF)
# CSF is used here to remove elevated vegetation and retain a point set that
# better represents terrain and lying deadwood surfaces.
las_lastreturn <- classify_ground(
  las,
  csf(
    sloop_smooth = TRUE,
    class_threshold = threshold_removevegetation,
    cloth_resolution = 0.5,
    rigidness = 1L,
    iterations = 500L,
    time_step = 1
  ),
  last_returns = TRUE
)

# Retain points classified as deadwood/class 2
las_no_trees_CSF_canopy <- filter_poi(las_lastreturn, Classification == 2)
remove(las)

# Normalize filtered point cloud
# Converts absolute elevation (Z) into height above deadwood using the terrain model.
las_no_trees_CSF_canopy <- normalize_height(las_no_trees_CSF_canopy, DTM_res)

# Optional RGBVI-based vegetation filtering
# RGBVI helps distinguish vegetation from deadwood using RGB values.
if (use_rgbvi_filter) {
  # Detect RGB bit depth automatically from the maximum stored value
  rgb_max <- max(
    c(
      las_no_trees_CSF_canopy@data$R,
      las_no_trees_CSF_canopy@data$G,
      las_no_trees_CSF_canopy@data$B
    ),
    na.rm = TRUE
  )
  
  # Scale according to detected 8-bit or 16-bit storage
  rgb_scale <- if (rgb_max <= 255) 255 else 65535
  
  # Normalize RGB channels to the range 0 to 1
  R <- as.numeric(las_no_trees_CSF_canopy@data$R) / rgb_scale
  G <- as.numeric(las_no_trees_CSF_canopy@data$G) / rgb_scale
  B <- as.numeric(las_no_trees_CSF_canopy@data$B) / rgb_scale
  
  # Compute RGBVI
  las_no_trees_CSF_canopy@data$RGBVI <- ((R * B) - (G^2)) / ((R * B) + (G^2))
  
  # Remove remaining vegetation based on height and RGBVI
  las_filtered <- filter_poi(
    las_no_trees_CSF_canopy,
    !(
      (Z > 0.05 & RGBVI <= rgbvi_threshold)
    )
  )
} else {
  # Keep all CSF-classified points without additional RGBVI filtering
  las_filtered <- las_no_trees_CSF_canopy
}

# Rebuild LAS object with corrected elevation field
# Zref is renamed to Z for subsequent raster-based processing.
df <- las_filtered@data[, c("X", "Y", "Zref")]
colnames(df)[3] <- "Z"
las_fixed <- LAS(df, header = las_filtered@header)

# Restore RGB attributes
las_fixed@data$R <- las_filtered@data$R
las_fixed@data$G <- las_filtered@data$G
las_fixed@data$B <- las_filtered@data$B

# ==================================================================================================
# Terrain and deadwood surface models
# ==================================================================================================

# Create a surface model with trees removed
#Use of 80th percentile
DSM_no_trees <- rast(grid_metrics(las_fixed, res = res_fin, ~quantile(Z, 0.8)))
DTM_res_clipped <- crop(DTM_res, DSM_no_trees)
DSM_no_trees_filled <- cover(DSM_no_trees, DTM_res_clipped)

# Remove intermediate point cloud to free memory
remove(las_no_trees_CSF_canopy)

# Compute vegetation height model relative to terrain
VHM <- DSM_no_trees_filled - DTM_res
VHM <- VHM * (VHM > -1)

# Rasterize filtered tree tops to represent stem of standing trees in roughness calculation
ttops_vect <- vect(ttops_filtered)
DSM_terra <- DSM_no_trees
ttops_raster_buffer <- rasterize(ttops_vect, DSM_terra, field = "Z", fun = "max", background = NA)

# Separate VHM representations with and without tree tops
VHM_onlyttops <- VHM * 0
VHM_onlyttops <- cover(ttops_raster_buffer, VHM_onlyttops)
VHM_aftertreetop <- cover(ttops_raster_buffer, VHM)

# ==================================================================================================
# Snow redistribution zones
# ==================================================================================================

# Define analysis extent for the redistribution grid
# The redistribution grid partitions the study area into zones used to estimate
# representative snow redistribution effects.
VHM_extent <- ext(VHM_crowns_filled)

# Align extent to the redistribution grid size
xmin_new <- floor(VHM_extent$xmin / res_redistribution) * res_redistribution
ymin_new <- floor(VHM_extent$ymin / res_redistribution) * res_redistribution
xmax_new <- ceiling(VHM_extent$xmax / res_redistribution) * res_redistribution
ymax_new <- ceiling(VHM_extent$ymax / res_redistribution) * res_redistribution

# Compute grid dimensions
ncol_new <- (xmax_new - xmin_new) / res_redistribution
nrow_new <- (ymax_new - ymin_new) / res_redistribution

# Create redistribution raster
zones_redistribution <- rast(
  ncols = ncol_new, nrows = nrow_new,
  xmin = xmin_new, xmax = xmax_new,
  ymin = ymin_new, ymax = ymax_new,
  crs = crs(VHM_crowns_filled)
)

# Convert grid to polygons and assign zone IDs
polygony_redistribution <- as.polygons(zones_redistribution)
polygony_redistribution$zones <- seq_len(nrow(polygony_redistribution))

# Rasterize zone IDs
zones_redistribution <- rasterize(polygony_redistribution, zones_redistribution, field = "zones")

# Extend VHM to the redistribution grid
VHM_crowns_filled_extended <- extend(VHM_crowns_filled, zones_redistribution)
zones_redistribution_02 <- resample(zones_redistribution, VHM_crowns_filled_extended, method = "near")

# Identify valid zones with sufficient data coverage
# Zones with too little valid data are excluded to avoid unstable quantile and
# roughness estimates in poorly covered areas.
zone_cover <- zonal(VHM_crowns_filled_extended > -5, zones_redistribution_02, fun = "sum", na.rm = TRUE)
valid_zones <- zone_cover$zones[zone_cover[, 2] >= 6000]

# Keep only valid redistribution zones
zones_redistribution_renumbered <- ifel(zones_redistribution_02 %in% valid_zones, zones_redistribution_02, NA)
zones_redistribution_renumbered <- resample(zones_redistribution_renumbered, zones_redistribution, method = "near")
zone_vals <- values(zones_redistribution_renumbered)

# Mask invalid redistribution zones
zones_redistribution_renumbered <- ifel(
  zones_redistribution %in% valid_zones,
  zones_redistribution,
  NA
)

# ==================================================================================================
# Aggregation grids for fuzzy logic and canopy coverage
# ==================================================================================================

# Create 5 m analysis grid
# This grid is used for fuzzy membership aggregation and for storing most final
# membership outputs at a practitioner-friendly spatial resolution.
zones_5m <- rast(
  ext(zones_redistribution_renumbered),
  resolution = res_fuzzy_logic,
  crs = crs(zones_redistribution_renumbered)
)
values(zones_5m) <- 1:ncell(zones_5m)
zones_mask_5m <- resample(zones_redistribution_renumbered, zones_5m, method = "near")
zones_5m <- mask(zones_5m, zones_mask_5m)

# Create 10 m grid for canopy coverage analysis
# A coarser grid is used here because canopy coverage is treated as an aggregated
# zone-based variable rather than as a very high-resolution surface.
zones_10m <- rast(
  ext(zones_redistribution_renumbered),
  resolution = 10,
  crs = crs(zones_redistribution_renumbered)
)
values(zones_10m) <- 1:ncell(zones_10m)
zones_mask_10m <- resample(zones_redistribution_renumbered, zones_10m, method = "near")
zones_10m <- mask(zones_10m, zones_mask_10m)

# Free memory before parallel computation
gc()

# ==================================================================================================
# Parallel computation of zone-based SVH metrics
# ==================================================================================================

# Define number of cores and register cluster
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export variables required by parallel workers
clusterExport(
  cl,
  varlist = c(
    "upfilling_quantile_end", "upfilling_quantile_start", "upfilling_quantile_step",
    "results_path", "zones_redistribution"
  )
)

# Save intermediate rasters used by parallel workers
writeRaster(VHM_aftertreetop, filename = file.path(results_path, "VHM/deadwood_VHM.tif"), overwrite = TRUE)
writeRaster(VHM, filename = file.path(results_path, "VHM/beforettop_VHM.tif"), overwrite = TRUE)
writeRaster(VHM_onlyttops, filename = file.path(results_path, "VHM/VHM_onlyttops.tif"), overwrite = TRUE)

if (!file.exists(file.path(results_path, "zones/zones_redistribution.tif"))) {
  writeRaster(zones_redistribution_renumbered, file.path(results_path, "zones/zones_redistribution.tif"), overwrite = TRUE)
}

# Compute SVH metrics across quantiles and redistribution zones in parallel
# For each redistribution zone, representative deadwood height is estimated from
# quantiles and translated into stored volume height metrics under varying snow
# filling assumptions.
results_SVH <- foreach(
  quantile_deadwood_height = seq(upfilling_quantile_start, upfilling_quantile_end, by = upfilling_quantile_step),
  .combine = "rbind",
  .packages = c("terra", "spatialEco"),
  .export = c("results_path")
) %dopar% {
  
  # Load intermediate rasters in each worker
  deadwood_VHM <- rast(file.path(results_path, "VHM/deadwood_VHM.tif"))
  zones_redistribution_renumbered <- rast(file.path(results_path, "zones/zones_redistribution.tif"))
  
  #Function: Snow accumulation on protruding deadwood elements
  snow_above_deadwood <- function(x) {
    k <- 0.916
    1 - exp(-k * x)
  }
  
  # Function:Snow sliding below deadwood
  snow_below_deadwood <- function(x) {
    ifelse(
      x < 0.4,
      0,
      ifelse(x < 1.5, 0.5455 * (x - 0.4), 0.6)
    )
  }
  
  df_results_local <- list()
  zone_vals <- values(zones_redistribution_renumbered)
  
  # Group raster cells by zone ID
  cells_by_zone <- split(
    which(!is.na(zone_vals)),
    zone_vals[!is.na(zone_vals)]
  )
  
  # Build extents for each zone
  rx <- res(zones_redistribution_renumbered)[1] / 2
  ry <- res(zones_redistribution_renumbered)[2] / 2
  
  zone_extents <- lapply(cells_by_zone, function(cell_idx) {
    xy <- xyFromCell(zones_redistribution_renumbered, cell_idx[1])
    ext(
      xy[1] - rx, xy[1] + rx,
      xy[2] - ry, xy[2] + ry
    )
  })
  
  names(zone_extents) <- names(cells_by_zone)
  
  for (zone_id in names(zone_extents)) {
    extent <- zone_extents[[zone_id]]
    zone_buffer <- extend(extent, 1.2)
    VHM_crop <- crop(deadwood_VHM, zone_buffer)
    
    # Compute quantile-based deadwood height
    VHM_max_value <- global(
      VHM_crop,
      fun = function(x) quantile(x, probs = quantile_deadwood_height, na.rm = TRUE)
    )[[1]]
    
    VHM_max <- VHM_crop
    values(VHM_max) <- VHM_max_value
    
    # Compute stored volume height (SVH)
    SVH <- VHM_max - VHM_crop
    SVH[SVH < 0] <- 0
    
    mean_SVH_value <- global(SVH, fun = "mean", na.rm = TRUE)[[1]]
    
    # Adjust SVH using snow functions
    # These functions represent how snow accumulates on protruding deadwood and slides below deadwood 
    values(SVH)[values(SVH) == 0] <- snow_above_deadwood(mean_SVH_value)
    values(SVH) <- values(SVH) + snow_below_deadwood(values(VHM_crop))
    
    mean_SVH <- global(SVH, fun = "mean", na.rm = TRUE)[[1]]
    
    # Store results
    df_results_local[[length(df_results_local) + 1]] <- data.frame(
      Quantile = quantile_deadwood_height,
      Quantile_height = VHM_max_value,
      Zone = zone_id,
      mean = mean_SVH
    )
    
    rm(zone_buffer, VHM_crop, VHM_max, SVH)
    gc()
  }
  
  do.call(rbind, df_results_local)
}

# Combine all SVH results in a single data frame
df_results_SVH_zone_quantile <- do.call(rbind, results_SVH)

# Stop cluster and free memory
stopCluster(cl)
closeAllConnections()
rm(cl)
gc()

# ==================================================================================================
# Membership functions
# ==================================================================================================

# Generalized bell-shaped membership for roughness of the bed surface
# This function transforms roughness values into fuzzy membership values, where
# higher membership indicates conditions more favorable for avalanche release.
bell_membership_roughness_bed_surface <- function(x, a = 0.01, b = 2, c = -0.005) {
  1 / (1 + abs((x - c) / a)^(2 * b))
}

# Slope membership function
# This function assigns high membership to slope angles typically associated with
# avalanche release and suppresses values outside the critical slope range.
bell_membership_slope <- function(x, a = 7, b = 4, c = 43) {
  slope <- 1 / (1 + abs((x - c) / a)^(2 * b))
  slope[x < 32 | x > 55] <- 0
  slope
}

# Canopy-coverage membership function
# This function describes the moderating effect of canopy cover on avalanche
# release, with denser canopy generally reducing release likelihood.
bell_membership_canopycoverage <- function(x, a = 40, b = 3.5, c = -15) {
  1 / (1 + abs((x - c) / a)^(2 * b))
}

# Compute slope membership on the 5 m grid
dtm_resampled_5m <- resample(DTM, zones_5m, method = "bilinear")
DTM_5m_slope <- terrain(dtm_resampled_5m, v = "slope", unit = "degrees")
slope_membership <- bell_membership_slope(DTM_5m_slope, a = 7, b = 4, c = 43)

# Save slope membership if it does not yet exist
if (!file.exists(file.path(results_path, "slope_membership/slope_ALS2015_5m_membership.tif"))) {
  writeRaster(slope_membership, file.path(results_path, "slope_membership/slope_ALS2015_5m_membership.tif"), overwrite = TRUE)
}

# ==================================================================================================
# Canopy coverage membership
# ==================================================================================================

# Prepare crown-height raster relative to deadwood surface
VHM_crowns_filled_crop <- crop(VHM_crowns_filled_extended, ext(zones_redistribution_renumbered))
r_template <- rast(VHM_crowns_filled_crop)
res(r_template) <- 0.2

VHM_beforettop_extended <- extend(VHM, zones_10m)
crown_height_above_deadwood <- VHM_crowns_filled_crop - VHM_beforettop_extended

# Classify canopy cover presence above 4 m (2 times expected 30y-snow depth)
crown_height_binary <- ifel(crown_height_above_deadwood < 4, 0, 1)
crown_height_binary <- extend(crown_height_binary, zones_10m)
zones_10m_02res <- resample(zones_10m, crown_height_binary, method = "near")

# Compute canopy cover fraction per 10 m zone
# Canopy cover is estimated as the proportion of cells with crown height above a
# defined threshold (4m), summarized for each canopy analysis zone.
zone_stats <- zonal(crown_height_binary, zones_10m_02res, fun = "sum", na.rm = TRUE)
zones_10m_counts <- zonal(!is.na(VHM_crowns_filled_extended), zones_10m_02res, fun = "sum", na.rm = TRUE)
zone_stats[, 2] <- zone_stats[, 2] / zones_10m_counts[, 2] * 100

# Convert canopy cover to raster and membership values
raster_canopycoverage <- classify(
  zones_10m,
  zone_stats,
  others = NA
)

membership_canopycoverage <- bell_membership_canopycoverage(raster_canopycoverage)
membership_canopycoverage <- resample(membership_canopycoverage, zones_5m, method = "near")

# Save intermediate canopy outputs
writeRaster(zones_5m, file.path(results_path, "zones/zones_5m.tif"), overwrite = TRUE)
writeRaster(membership_canopycoverage, file.path(results_path, "adapted_tree_parameters/canopycoverage_membership.tif"), overwrite = TRUE)

# Free memory before fuzzy-logic processing
gc()

# ==================================================================================================
# Parallel fuzzy-logic computation of avalanche release membership
# This is the main snow-depth dependent modeling step. Roughness, slope, and
# canopy cover memberships are combined for each modeled snow depth.
# ==================================================================================================

# Define number of cores for fuzzy-logic modeling
max_cores <- 41
num_cores <- min(num_cores, max_cores)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export required objects to parallel workers
clusterExport(
  cl,
  varlist = c(
    "upfilling_SVH_start", "upfilling_SVH_end", "upfilling_SVH_step",
    "results_SVH", "bell_membership_roughness_bed_surface", "results_path"
  )
)

results_fuzzy_approach <- foreach(
  SVH_height = seq(upfilling_SVH_start, upfilling_SVH_end, by = upfilling_SVH_step),
  .packages = c("terra", "foreach", "doParallel", "spatialEco", "dplyr")
) %dopar% {
  
  # Load required rasters in each worker
  VHM_withttop <- rast(file.path(results_path, "VHM/deadwood_VHM.tif"))
  VHM_onlyttop <- rast(file.path(results_path, "VHM/VHM_onlyttops.tif"))
  VHM_beforettop <- rast(file.path(results_path, "VHM/beforettop_VHM.tif"))
  zones_redistribution_renumbered <- rast(file.path(results_path, "zones/zones_redistribution.tif"))
  slope_membership <- rast(file.path(results_path, "slope_membership/slope_ALS2015_5m_membership.tif"))
  zones_5m <- rast(file.path(results_path, "zones/zones_5m.tif"))
  membership_canopycoverage <- rast(file.path(results_path, "adapted_tree_parameters/canopycoverage_membership.tif"))
  
  # Define target raster for roughness computation
  target_raster <- zones_redistribution_renumbered
  res(target_raster) <- c(0.4, 0.4)
  
  # Select the closest quantile-based deadwood height per zone
  closest_quantile_per_zone <- results_SVH %>%
    group_by(Zone) %>%
    slice(which.min(abs(mean - SVH_height))) %>%
    ungroup() %>%
    dplyr::select(Zone, Quantile_height, mean, Quantile)
  
  roughness_list <- list()
  VHM_filled_list <- list()
  
  zone_vals <- values(zones_redistribution_renumbered)
  
  # Group cells by redistribution zone
  cells_by_zone <- split(
    which(!is.na(zone_vals)),
    zone_vals[!is.na(zone_vals)]
  )
  
  # Build extents for each zone
  rx <- res(zones_redistribution_renumbered)[1] / 2
  ry <- res(zones_redistribution_renumbered)[2] / 2
  
  zone_extents <- lapply(cells_by_zone, function(cell_idx) {
    xy <- xyFromCell(zones_redistribution_renumbered, cell_idx[1])
    ext(
      xy[1] - rx, xy[1] + rx,
      xy[2] - ry, xy[2] + ry
    )
  })
  
  names(zone_extents) <- names(cells_by_zone)
  
  for (zone_id in names(zone_extents)) {
    extent <- zone_extents[[zone_id]]
    zone_buffer <- extend(extent, 1.6)
    
    VHM_crop <- crop(VHM_withttop, zone_buffer)
    VHM_max <- VHM_crop
    
    # Assign quantile-based maximum deadwood height for the zone
    VHM_max_value <- closest_quantile_per_zone %>%
      filter(Zone == zone_id) %>%
      pull(Quantile_height)
    
    values(VHM_max) <- VHM_max_value
    
    # Prepare raster variants for roughness analysis
    VHM_onlyttop_crop <- crop(VHM_onlyttop, zone_buffer)
    VHM_beforettop_crop <- crop(VHM_beforettop, zone_buffer)
    VHM_withttop_crop <- crop(VHM_withttop, zone_buffer)
    
    
    # Combine the local VHM with the zone-specific snow cover.
    # This simulates a snow surface that fills depressions up to the
    # simulated snow height within each redistribution zone.
    VHM_onlyttop_filled <- app(c(VHM_max, VHM_onlyttop_crop), fun = max, na.rm = TRUE)
    VHM_filled_beforettop <- app(c(VHM_max, VHM_beforettop_crop), fun = max, na.rm = TRUE)
    VHM_filled_withttop <- app(c(VHM_max, VHM_withttop_crop), fun = max, na.rm = TRUE)
    
    
    # Calculate still outstanding deadwood or stems height
    # roughness-relevant residual structure.
    VHM_onlyttop_filled <- VHM_onlyttop_filled - VHM_max
    VHM_filled_beforettop <- VHM_filled_beforettop - VHM_max
    VHM_filled_withttop <- VHM_filled_withttop - VHM_max
    
    # Resample for roughness calculation
    VHM_onlyttop_filled_resampled <- resample(VHM_onlyttop_filled, target_raster, method = "bilinear")
    VHM_filled_resampled_beforettop <- resample(VHM_filled_beforettop, target_raster, method = "bilinear")
    
    VHM_onlyttops_filled_resampled <- crop(VHM_onlyttop_filled_resampled, zone_buffer)
    VHM_filled_resampled_beforettop <- crop(VHM_filled_resampled_beforettop, zone_buffer)
    
    # Compute roughness for winter terrain without standing trees and only for standing trees (different window size) 
	#as spatial roughness effect is smaller of standing trees
    roughness_onlyttop <- vrm(VHM_onlyttops_filled_resampled, 3)
    roughness_beforettop <- vrm(VHM_filled_resampled_beforettop, 7)
    
    #combine roughness values using higher value
    roughness_combined <- app(
      c(roughness_onlyttop, roughness_beforettop),
      fun = max,
      na.rm = TRUE
    )
    
    # Restrict outputs to the original zone extent
    roughness_clipped <- crop(roughness_combined, extent)
    VHM_crop_filled <- crop(VHM_filled_withttop, extent)
    
    roughness_list[[zone_id]] <- roughness_clipped
    VHM_filled_list[[zone_id]] <- VHM_crop_filled
  }
  
  # Merge zone-wise outputs
  roughness_list <- roughness_list[!sapply(roughness_list, is.null)]
  roughness_mosaic <- do.call(mosaic, c(unname(roughness_list), list(fun = "first")))
  roughness_mosaic <- extend(roughness_mosaic, ext(zones_5m))
  roughness_mosaic_membership <- bell_membership_roughness_bed_surface(roughness_mosaic)
  
  zones_05m <- resample(zones_5m, roughness_mosaic, method = "near")
  
  # Aggregate roughness membership by zone using different quantiles (5,10,20)
  zone_q5 <- zonal(
    roughness_mosaic_membership,
    zones_05m,
    fun = function(x) quantile(x, probs = 0.05, na.rm = TRUE)
  )
  roughness_aggregated_q5 <- zones_5m
  values(roughness_aggregated_q5) <- zone_q5[, 2][match(values(zones_5m), zone_q5[, 1])]
  
  zone_q10 <- zonal(
    roughness_mosaic_membership,
    zones_05m,
    fun = function(x) quantile(x, probs = 0.10, na.rm = TRUE)
  )
  roughness_aggregated_q10 <- zones_5m
  values(roughness_aggregated_q10) <- zone_q10[, 2][match(values(zones_5m), zone_q10[, 1])]
  
  zone_q20 <- zonal(
    roughness_mosaic_membership,
    zones_05m,
    fun = function(x) quantile(x, probs = 0.20, na.rm = TRUE)
  )
  roughness_aggregated_q20 <- zones_5m
  values(roughness_aggregated_q20) <- zone_q20[, 2][match(values(zones_5m), zone_q20[, 1])]
  
  # Merge filled VHM outputs
  VHM_filled_list <- VHM_filled_list[!sapply(VHM_filled_list, is.null)]
  VHM_filled_mosaic <- do.call(mosaic, c(unname(VHM_filled_list), list(fun = "first")))
  
  # Fuzzy aggregation of roughness, slope, and canopy membership
  # The aggregation combines limiting effects and average conditions, yielding a
  # balanced fuzzy estimate of avalanche release likelihood.
  mu_PRA <- function(roughness_input, slope_input, canopy_input) {
    mu_stack <- c(roughness_input, slope_input, canopy_input)
    mu_aggregated <- app(mu_stack, fun = function(x) {
      if (any(is.na(x))) {
        return(NA)
      } else {
        mu_min <- min(x)
        gamma <- 1 - mu_min
        gamma * mu_min + (1 - gamma) * mean(x)
      }
    })
    mu_aggregated
  }
  #calculate PRA for different aggregation percentiles (5,10,20)
  fuzzy_logic_result_q5roughness <- mu_PRA(roughness_aggregated_q5, slope_membership, membership_canopycoverage)
  fuzzy_logic_result_q10roughness <- mu_PRA(roughness_aggregated_q10, slope_membership, membership_canopycoverage)
  fuzzy_logic_result_q20roughness <- mu_PRA(roughness_aggregated_q20, slope_membership, membership_canopycoverage)
  
  # Define output filenames
  output_filename_roughness <- file.path(results_path, "roughness", paste0("roughness_SVH_", SVH_height, ".tif"))
  output_filename_roughness_membership <- file.path(results_path, "roughness_membership", paste0("roughness_membership_SVH_", SVH_height, ".tif"))
  output_filename_roughness_aggregated_q5 <- file.path(results_path, "roughness_membership", paste0("roughness_aggregated_q5_SVH_", SVH_height, ".tif"))
  output_filename_roughness_aggregated_q10 <- file.path(results_path, "roughness_membership", paste0("roughness_aggregated_q10_SVH_", SVH_height, ".tif"))
  output_filename_roughness_aggregated_q20 <- file.path(results_path, "roughness_membership", paste0("roughness_aggregated_q20_SVH_", SVH_height, ".tif"))
  output_filename_VHM_filled <- file.path(results_path, "VHM_filled", paste0("VHM_filled_", SVH_height, ".tif"))
  output_filename_fuzzy_logic_q5 <- file.path(results_path, "fuzzy_logic", paste0("fuzzy_logic_result_q5_SVH_", SVH_height, ".tif"))
  output_filename_fuzzy_logic_q10 <- file.path(results_path, "fuzzy_logic", paste0("fuzzy_logic_result_q10_SVH_", SVH_height, ".tif"))
  output_filename_fuzzy_logic_q20 <- file.path(results_path, "fuzzy_logic", paste0("fuzzy_logic_result_q20_SVH_", SVH_height, ".tif"))
  
  # Save outputs for the current snow depth
  writeRaster(roughness_mosaic_membership, output_filename_roughness_membership, overwrite = TRUE)
  writeRaster(roughness_mosaic, output_filename_roughness, overwrite = TRUE)
  writeRaster(roughness_aggregated_q20, output_filename_roughness_aggregated_q20, overwrite = TRUE)
  writeRaster(roughness_aggregated_q10, output_filename_roughness_aggregated_q10, overwrite = TRUE)
  writeRaster(roughness_aggregated_q5, output_filename_roughness_aggregated_q5, overwrite = TRUE)
  writeRaster(VHM_filled_mosaic, output_filename_VHM_filled, overwrite = TRUE)
  writeRaster(fuzzy_logic_result_q5roughness, output_filename_fuzzy_logic_q5, overwrite = TRUE)
  writeRaster(fuzzy_logic_result_q10roughness, output_filename_fuzzy_logic_q10, overwrite = TRUE)
  writeRaster(fuzzy_logic_result_q20roughness, output_filename_fuzzy_logic_q20, overwrite = TRUE)
  
  output_filename_roughness
}

# Stop cluster and free memory
stopCluster(cl)
closeAllConnections()
rm(cl)
gc()

# Save final outputs
# Please adapt file names if a clearer distinction between parameter settings is required.
# These final products include the deadwood-related surface model, crown raster
# outputs, and additional structural rasters needed for interpretation.
writeRaster(DSM_no_trees, file.path(results_path, "VHM", "DSM_notrees.tif"), overwrite = TRUE)
writeRaster(crowns, file.path(results_path, "adapted_tree_parameters", "crowns.tif"), overwrite = TRUE)

# Save filtered tree tops as shapefile
ttops <- st_as_sf(ttops_filtered)
st_write(ttops, file.path(results_path, "adapted_tree_parameters", "tree_tops.shp"), delete_layer = TRUE)

# Save additional VHM products
writeRaster(VHM, file.path(results_path, "VHM", "VHM_notrees.tif"), overwrite = TRUE)
writeRaster(VHM_crowns_filled, file.path(results_path, "VHM", "VHM_crowns_filled.tif"), overwrite = TRUE)
writeRaster(crown_height_above_deadwood, file.path(results_path, "VHM", "crown_height_above_deadwood.tif"), overwrite = TRUE)
