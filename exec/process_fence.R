# Processing fence files for use in NEON analysis
library(terra)
library(sf)
library(ggplot2)
import::from(magrittr, "%>%")
library(dplyr)

shp_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/spatial/NEON_Shapefiles'
fence_fl <- file.path(shp_dir, '/HI_completed_fencelines/HI_completed_fencelines.shp')

# get flow CRS
flow <- sf::st_read(file.path(shp_dir, 'Flow Map', 'Geological_Units.shp'))
flow_crs <- terra::crs(flow)
rm(flow); gc()

# get NAR overlay
outline <- terra::vect(file.path(shp_dir, 'PUUM_Outline', 'PUUM_OUTLINE.shp'))
outline <- terra::project(outline, flow_crs)

# get fence file
fence <- terra::vect(fence_fl)
fence <- terra::project(fence, flow_crs)
# crop fence to NAR
fence <- terra::crop(fence, ext(outline))
# drop outline secondary geometry
#outline <- outline[1]
#outline <- sf::st_as_sf(outline)
#outline <- outline[, 'geometry']

# save fence file
#terra::writeVector(fence, file.path(shp_dir, 'NAR_fencelines_raw.shp'), filetype = 'ESRI Shapefile')
rm(fence); gc()

# QGIS work:
#     1) save the lines as point vertices
#     2) re-order each polygon by hand, save to a .csv

p_df <- read.csv(file.path(shp_dir, "../fencelines_modified/point_notes.csv"))
#p_df <- p_df[which(p_df$polyID == 8), ]
# get fenceline file, now as vertices
fence <- terra::vect(file.path(shp_dir, '../fencelines_modified/fencelines_vertices.shp'))
fence <- fence[, -which(names(fence) %in% c('polyID', 'withinPoly'))]

# drop some for testing
fence <- fence[which(fence$fid %in% p_df$fid), ]
# random extras popped in
fence <- fence[-which(fence$fid == 1586), ]
fence <- fence[-which(fence$fid == 2567), ]
# add the metadata
fence_crs <- crs(fence)
fence <- as.data.frame(fence, geom = 'WKT')
fence <- dplyr::left_join(fence, p_df, by = 'fid')
fence_meta <- within(fence, rm(fid, Island, StatusFL, CompletedD, SecondaryM, TertiaryMa, BirdingTap, LineSource, SHAPE_Leng, vertex_ind, vertex_par, vertex_p_1, vertex_p_2, distance, angle, geometry, withinPolyOrder))
fence_meta <- fence_meta[!duplicated(fence_meta), ]

fence <- terra::vect(fence, geom = 'geometry', crs = fence_crs) %>%
  terra::project(flow_crs) %>%
  sf::st_as_sf() %>%
  dplyr::group_by(polyID) %>%
  dplyr::arrange(withinPolyOrder) %>%
  dplyr::summarize(do_union = F) %>%
  sf::st_cast('POLYGON')

# dropped meta, add it back in
fence <- dplyr::left_join(fence, fence_meta, by = 'polyID')

# save resulting polygon file
sf::st_write(fence, file.path(shp_dir, '../fencelines_modified/fenceline_polygons_NAR.shp'), driver = 'ESRI Shapefile', append = F)

