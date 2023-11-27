# BEM Nov 2023
library(NEON1)
library(terra)

shp_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/spatial/NEON_Shapefiles'
fence_fl <- file.path(shp_dir, '../fencelines_modified/fenceline_polygons_NAR.shp')

# get fence file
fence <- terra::vect(fence_fl)

# get plots
plots <- st_read(file.path(shp_dir, 'PUUM_Plots_fixed', 'AllPlot_PUUM.shp'))
# # 21-30 are mosquito, but 5 removed, 5 added as 50s
plots <- plots[-which(plots$layer == "MOS_PLOTS.SHP"), ]
# 1-20 is distributed, including 9, 17, 18
plots$layer[is.na(plots$layer)] <- "Distributed_Plots_SHP"
# 31-50 are tower plots, should be included
incl_plots <- sapply(strsplit(plots$Name_2, '_'), \(xx) as.numeric(xx[2]))
incl_plots[order(incl_plots)]
stopifnot(!any(c(21:30, 51:100) %in% incl_plots))
plots <- st_transform(plots, crs(fence))
plots <- st_crop(plots, fence)

# need more fence metadata from Mike
