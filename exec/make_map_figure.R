library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(NEON1)

# also creates data object for flow age

shp_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/spatial/NEON_Shapefiles'
fig_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures'

# fix roads shapefile
#roads <- vect(file.path(shp_dir, 'PUUM_Roads', 'PUUM_roads_edited.shp'))
#writeVector(roads, file.path(shp_dir, 'PUUM_Roads', 'PUUM_roads_fixed.shp'))
#rm(roads); gc()

# This is the flow shapefile for all of Hawaii
flow <- st_read(file.path(shp_dir, 'Flow Map', 'Geological_Units.shp'))
# Drop the 1984 flow since it doesnt really intersect any of the plots or NAR
flow <- flow[-which(flow$age_range == 'A.D. 1984'), ]
# This is the outline for the Pu'u Maka'ala Natural Area Reserve
# The cut-out is the Kulani Correctional Facility
outline <- st_read(file.path(shp_dir, 'PUUM_Outline', 'PUUM_OUTLINE.shp'))
outline <- st_transform(outline, crs(flow))
outline <- outline[1, ]
roads <- st_read(file.path(shp_dir, 'PUUM_Roads', 'PUUM_roads_fixed.shp'))
roads <- st_transform(roads, crs(flow))
roads <- st_crop(roads, outline)
flow <- st_crop(flow, outline)
#flow$age_val <- key_age_range(flow$age_range)
plots <- st_read(file.path(shp_dir, 'PUUM_Plots_fixed', 'AllPlot_PUUM.shp'))
# # 21-30 are mosquito, but 5 removed, 5 added as 50s
plots <- plots[-which(plots$layer == "MOS_PLOTS.SHP"), ]
# 1-20 is distributed, including 9, 17, 18
plots$layer[is.na(plots$layer)] <- "Distributed_Plots_SHP"
# 31-50 are tower plots, should be included
incl_plots <- sapply(strsplit(plots$Name_2, '_'), \(xx) as.numeric(xx[2]))
incl_plots[order(incl_plots)]
stopifnot(!any(c(21:30, 51:100) %in% incl_plots))
plots <- st_transform(plots, crs(flow))
plots <- st_crop(plots, outline)

# what flow data is associated with what plots?
flow_meta <- flow[unlist(st_within(plots, flow)), c('strat_code', 'symbol', 'age_group', 'age_range')]
flow_meta <- st_drop_geometry(flow_meta)
flow_meta <- data.frame(plotID = st_drop_geometry(plots$Name_2), flow_meta)
row.names(flow_meta) <- NULL
usethis::use_data(flow_meta, overwrite = T)

tower <- st_read(file.path(shp_dir, 'PUUM_Tower', 'Tower.shp'))
tower <- st_transform(tower, crs(flow))
tower <- st_crop(tower, flow)

#moist <- st_read(file.path(shp_dir, 'PUUM_Moisture', 'Moist_Clip.shp'))
#moist <- st_transform(moist, crs(flow))
#moist <- st_crop(moist, flow)

rain <- st_read(file.path(shp_dir, 'Anual Rain', 'Annual_Rainfall_(in).shp'))
rain <- st_transform(rain, crs(flow))
rain <- st_crop(rain, flow)
rain$contour_mm <- round(rain$contour * 25.4)
rain$contour_m <- round(rain$contour_mm / 1000, 1)
rain$contour_m <- format(rain$contour_m, digits = 2)
# drop some layers to fit on plot better
rain <- rain[-c(2, 4, 6, 8), ]

KCF <- st_read(file.path(shp_dir, 'KCF_Prison_Outline', 'KCF_Outline.shp'))
KCF <- st_transform(KCF, crs(flow))
KCF <- st_crop(KCF, flow)

ggplot() +
  # Add flow layer - fill aesthetic
  geom_sf(data = flow, aes(fill = relevel_age(age_range)), show.legend = 'polygon') +
  # Fill in KCF plygon - no aesthetics
  geom_sf(data = KCF, aes(color = 'KCF'), fill = 'grey70') +
  # Add rainfall isolines with labels - label aesthetic
  geom_sf(data = rain, aes(color = 'rain'), linewidth = 0.5) +
  geom_sf_label(data = rain, aes(label = contour_m), size = 2.0, alpha = 0.8) +
  # Add outline, plots, and tower - no aesthetics
  geom_sf(data = outline, aes(color = 'NAR'), fill = NA, linewidth = 0.7) +
  geom_sf(data = plots, aes(color = 'plots'), fill = 'white', size = 3.2, shape = 21, alpha = 0.8) +
  geom_sf(data = tower, aes(color = 'tower'), fill = 'yellow', shape = 23, size = 3.8) +
  # add scale bar
  ggsn::scalebar(data = flow, location = 'bottomleft', anchor = c(x = -155.225, y = 19.485), dist = 1, dist_unit = 'km', st.size = 3, st.bottom = F, st.color = 'white', height = 0.025, transform = T, model = 'WGS84') +
  #ggsn::north(data = flow, location = 'topright', symbol = 16) +

  scale_fill_viridis(option = "A", discrete = T, direction = -1) +
  # Color scale - needs a white line, black line, black-outline circle, yellow square
  scale_color_manual(
    labels = c('KCF', 'NAR', 'Survey plots', 'MAP isolines', 'Canopy tower'),
    values = c('KCF' = 'grey70', 'NAR' = 'white', 'plots' = 'black', 'rain' = 'black', 'tower' = 'black'),
    guide = guide_legend(override.aes = list(
      linetype = c('solid', 'solid', 'blank', 'solid', 'blank'),
      shape = c(NA, NA, 1, NA, 23)
    ))
  ) +

  theme(
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black')
  ) +
  labs(fill = 'Flow Age', x = '', y = '', color = '')

ggsave(file.path(fig_dir, 'map_figure.pdf'), width = 11, height = 8.5, units = 'in')
ggsave(file.path(fig_dir, 'map_figure.png'), width = 11, height = 8.5, units = 'in')
