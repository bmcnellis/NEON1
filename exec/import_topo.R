library(NEON1)
library(terra)
library(neonUtilities)
library(Rsagacmd)

### Get plot metadata
met <- NEON1::met
pts <- terra::vect(data.frame(lon = met$decimalLongitude, lat = met$decimalLatitude))
crs(pts) <- crs('EPSG:4326')
pts <- terra::project(pts, NEON1::NEON_aop_crs)
met_out <- met[, c('plotID', 'decimalLatitude', 'decimalLongitude')]

### Set up SAGA
saga <- Rsagacmd::saga_gis(raster_backend = 'terra', cores = 8)

# Variables to calculate will be:
# Topographic wetness index
TWI_FUN <- saga$ta_hydrology$topographic_wetness_index_twi
# Requires slope and catchment area
SLOPE_FUN <- saga$ta_morphometry$slope_aspect_curvature
SCA_FUN <- saga$ta_hydrology$flow_width_and_specific_catchment_area
# catchment area requires flow accumulation
FAT_FUN <- saga$ta_hydrology$flow_accumulation_top_down
# VDCN requires channel network
chan_FUN <- saga$ta_channels$channel_network
# first fill sinks
#sink_FUN <- saga$ta_preprocessor$fill_sinks_wang_liu
# then flow-accumulation top down
# = FAT_FUN

#e1_dir <- '../spatial/AOP'
e1_dir <- '/media/bem/data/NEON/spatial/AOP'
# DP3.30024.001 is 3.36 GB downloaded 2024-08-08 12:31 PM HST
#e1 <- neonUtilities::byFileAOP('DP3.30024.001', 'PUUM', 2020, T, F, e1_dir)
aop_dir <- 'DP3.30024.001/neon-aop-products/2020/FullSite/D20/2020_PUUM_2/L3/DiscreteLidar/DTMGtif'

### Elevation DTM

tif_fls <- list.files(file.path(e1_dir, aop_dir), full.names = T)
td <- tempdir()

res <- data.frame()

for (i in seq_along(tif_fls)) {

  # load a raster
  ii <- tif_fls[i]
  irast <- terra::rast(ii)
  stopifnot(identical(crs(irast), crs(pts)))

  # see if any of the plots occur within that raster
  ie <- terra::extract(irast, pts)
  if (length(table(ie[, 2])) == 0) {
    next
  } else {
    cat(paste0('\nfile: ', ii, '\n'))
  }

  # if so, calculate the SAGA covariates and extract points
  irast_fl <- file.path(td, 'irast.tif')
  ifat_fl <- file.path(td, 'ifat.tif')
  isca_fl <- file.path(td, 'isca.tif')
  islope_fl <- file.path(td, 'islope.tif')
  terra::writeRaster(irast, irast_fl, filetype = 'GTiff', overwrite = T)
  ifat <- FAT_FUN(irast_fl)$flow
  terra::writeRaster(ifat, ifat_fl, filetype = 'GTiff', overwrite = T)
  iCN <- chan_FUN(irast_fl, init_grid = ifat)$chnlntwrk
  isca <- SCA_FUN(irast_fl, ifat_fl)$sca
  terra::writeRaster(isca, isca_fl, filetype = 'GTiff', overwrite = T)
  islope <- SLOPE_FUN(irast_fl)$slope
  terra::writeRaster(islope, islope_fl, filetype = 'GTiff', overwrite = T)
  iTWI <- TWI_FUN(islope_fl, isca_fl)

  # return the values
  irast_val <- terra::extract(irast, pts)
  ifat_val <- terra::extract(ifat, pts)
  isca_val <- terra::extract(isca, pts)
  islope_val <- terra::extract(islope, pts)
  iTWI_val <- terra::extract(iTWI, pts)

  imet_out <- cbind(met_out, elev = irast_val[, 2], fat = ifat_val[, 2], sca = isca_val[, 2], slope = islope_val[, 2], TWI = iTWI_val[, 2])
  imet_out <- imet_out[complete.cases(imet_out), ]

  if (length(res) == 0) {
    res <- imet_out
  } else {
    res <- rbind(res, imet_out)
  }

}

res <- res[!duplicated(res), ]

topo <- res
usethis::use_data(topo)

file.remove(list.files(td, full.names = T))
