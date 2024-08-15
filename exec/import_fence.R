# BEM Nov 2023
# updated Aug 2024

library(NEON1)
library(terra)
library(viridis)

res_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/reports'

data('met')

td <- tempdir()
shp_zip <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/data/PuuMakaala_Fence_Units_20240517.zip'
unzip(shp_zip, exdir = td)
fence_fl <- file.path(td, 'PuuMakaala_Fence_Units_20240517.shp')

# get fence file
fence <- terra::vect(fence_fl)
crs(fence) <- 'EPSG:32605'
# save fence metadata
met_fence <- fence
# program was mixed DOFAW-NEPM and DOFAW-NEPM-TMA, all records by dlnr_lshizuma
# DateUngFre was unpopulated
met_fence$EstComp <- grepl('estimate|estimation', met_fence$Comments)
met_fence$Ingress <- c('2023', rep(NA, nrow(met_fence) - 1))
met_fence <- met_fence[, -which(names(met_fence) %in% c('Reserve', 'Program', 'Comments', 'DateUngFre', 'EditDate', 'Editor', 'CreationDa', 'GlobalID', 'Creator', 'Display', 'Shape__Are', 'Shape__Len'))]
met_fence$Acres <- round(met_fence$Acres * 0.404686, 2)
met_fence$DateComp <- as.numeric(strptime(met_fence$DateComp, '%Y/%m/%d')$year + 1900)
names(met_fence) <- c('enclosure_name', 'enclosure_ha', 'enclosure_completed', 'enclosure_ung_free', 'enclosure_maint_days', 'enclosure_completed_estimate', 'enclosure_ingress')

# get plots
plots <- terra::vect(data.frame(plotID = met$plotID, lon = met$decimalLongitude, lat = met$decimalLatitude), crs = 'EPSG:4326')
plots <- terra::project(plots, crs(fence))

# add fence data to plot data and re-save
met$enclosure_name <- terra::extract(met_fence, plots)$enclosure_name
#usethis::use_data(met, overwrite = T)
fence <- met[, c('plotID', 'enclosure_name')]
fence <- fence[!duplicated(fence), ]
usethis::use_data(fence, overwrite = T)
fence0 <- met_fence
met_fence <- as.data.frame(met_fence)
usethis::use_data(met_fence, overwrite = T)

# make plots

fence0$enclosure_completed <- as.factor(fence0$enclosure_completed)

# which plots are missing dates?
#png(file.path(res_dir, 'enclosures_missing_complete_dates.png'), units = 'in', res = 300, height = 8, width = 8)
plot(fence0, 'enclosure_completed', col = viridis::viridis(14, option = 'G', direction = -1)[1:length(unique(fence0$enclosure_name))], type = 'classes', main = 'PUUM Fence Units & Estimated Completion Dates')
plot(fence0[which(fence0$enclosure_name == 'Kata-Stein'), ], add = T, col = 'red')
plot(fence0[which(fence0$enclosure_name == 'Unencumbered Kulani Pasture'), ], add = T, col = 'red')
plot(plots, add = T)
text(terra::vect(data.frame(lon = c(256000, 258000), lat = c(2165000, 2163500)), crs = crs(fence0)), labels = c('???', '???'))
#dev.off()

#terra::writeVector(fence0, file.path(res_dir, 'PUUM_fencing.shp'), filetype = 'ESRI Shapefile')
#write.csv(met, file.path(res_dir, 'PUUM_fencing.csv'), row.names = FALSE)
