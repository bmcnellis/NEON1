# BEM
# updated August 2024
# updated October 2026: added Figure 2 (disturbance map); renamed from make_map_figure.R

library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(NEON1)

# makes Figure 1 (flow age map) and Figure 2 (disturbance map)
# also creates data object for flow age

#shp_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/spatial/NEON_Shapefiles'
shp_dir <- '../spatial/NEON_Shapefiles'
#fig_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures'
fig_dir <- '../results/figures'

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
#st_write(outline[1], file.path(shp_dir, 'PUUM_outline_fixed.shp'))
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
# set flow meta median based on provided age range
flow_meta$age_median <- rep(NA, nrow(flow_meta))
flow_meta$age_median <- ifelse(flow_meta$age_range == '750-1,500 yr', 1125, flow_meta$age_median)
flow_meta$age_median <- ifelse(flow_meta$age_range == '1,500-3,000 yr', 2250, flow_meta$age_median)
flow_meta$age_median <- ifelse(flow_meta$age_range == '200-750 yr', 475, flow_meta$age_median)
flow_meta$age_median <- ifelse(flow_meta$age_range == '3,000-5,000 yr', 4000, flow_meta$age_median)
flow_meta$age_median <- ifelse(flow_meta$age_range == '5,000-11,000 yr', 8000, flow_meta$age_median)
#usethis::use_data(flow_meta, overwrite = T)

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

plots <- plots[which(plots$Name_2 %in% NEON1::incl_plots), ]
lus <- NEON1::lus
colnames(lus)[1] <- 'Name_2'
plots <- dplyr::left_join(plots, lus, by = 'Name_2')
plots_cow <- plots[plots$cow, ]
plots_log <- plots[plots$log, ]
#file.remove(file.path(shp_dir, 'NEON_plots_reduced.shp'))
#st_write(plots, file.path(shp_dir, 'NEON_plots_reduced.shp'))
#st_write(plots_cow, file.path(shp_dir, 'NEON_plots_cow.shp'))
#st_write(plots_log, file.path(shp_dir, 'NEON_plots_log.shp'))

map_flow <- ggplot() +
  # Add flow layer - fill aesthetic
  # ggplot2 4.0 errors on geom_sf(show.legend = 'polygon') ("`x` must be a vector,
  # not `NULL`"); TRUE gives the same filled-polygon key
  geom_sf(data = flow, aes(fill = relevel_age(age_range)), show.legend = TRUE) +
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
  #ggsn::scalebar(data = flow, location = 'bottomleft', anchor = c(x = -155.225, y = 19.485), dist = 1, dist_unit = 'km', st.size = 3, st.bottom = F, st.color = 'white', height = 0.025, transform = T, model = 'WGS84') +
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

ggsave(file.path(fig_dir, 'map_figure.pdf'), plot = map_flow, width = 11, height = 8.5, units = 'in')
ggsave(file.path(fig_dir, 'map_figure.png'), plot = map_flow, width = 11, height = 8.5, units = 'in')

### Figure 2: disturbance map

# Fence units come from the same state layer that import_fence.R reads, so the
# years-ungulate-free shown here are the values that enter the models.
td <- tempdir()
unzip('../data/PuuMakaala_Fence_Units_20240517.zip', exdir = td)
fence_u <- st_read(file.path(td, 'PuuMakaala_Fence_Units_20240517.shp'))
# The shipped .prj declares NAD83 / UTM 5N (EPSG:26905), and transforming out of
# it returns NA coordinates here because the NAD83 datum grid is not installed.
# import_fence.R works around this by declaring EPSG:32605 (WGS 84 / UTM 5N)
# instead; the same override is used here so the figure and the model data agree.
st_crs(fence_u) <- 32605
fence_u <- st_transform(fence_u, st_crs(flow))
# terra/sf >= 1.8 return the shapefile date field already parsed as POSIXct
fence_u$completed <- as.integer(format(as.POSIXct(fence_u$DateComp, tz = 'UTC'), '%Y'))
# Kata-Stein and Unencumbered Kulani Pasture have no published completion date;
# both are set to 2014 here, the same assumption import_fence.R makes
fence_u$date_estimated <- is.na(fence_u$completed)
fence_u$completed <- ifelse(fence_u$date_estimated, 2014, fence_u$completed)
# match the model predictor: 2024 - completion, and 0 where not yet ungulate-free
fence_u$ung_free <- 2024 - fence_u$completed
fence_u$ung_free <- ifelse(fence_u$UngFree == 'No', 0, fence_u$ung_free)
fence_u$ung_free <- factor(fence_u$ung_free, levels = sort(unique(fence_u$ung_free)))
fence_u <- st_crop(fence_u, outline)

# plot-level land-use history, the two binary disturbance predictors combined
plots$land_use <- ifelse(
  plots$cow & plots$log, 'Grazed and logged',
  ifelse(plots$cow, 'Grazed', ifelse(plots$log, 'Logged', 'Neither'))
)
plots$land_use <- factor(
  plots$land_use, levels = c('Neither', 'Grazed', 'Logged', 'Grazed and logged')
)
# readable overstory labels, reference level first
cover_labels <- c(
  'ohia_tall' = 'Mature ohia',
  'ohia_woodland' = 'Ohia woodland',
  'ohia_koa_tall' = 'Mature ohia-koa',
  'koa_tall' = 'Mature koa',
  'ohia_ash_tall' = 'Mature ohia-ash'
)
plots$cover_label <- factor(
  unname(cover_labels[plots$cover_type]), levels = unname(cover_labels)
)

# Panel a: ungulate exclosures and the years-ungulate-free predictor
map_fence <- ggplot() +
  geom_sf(data = fence_u, aes(fill = ung_free), color = 'grey30', linewidth = 0.3, show.legend = TRUE) +
  geom_sf(data = fence_u[fence_u$date_estimated, ], fill = NA, color = 'red', linewidth = 0.6) +
  geom_sf(data = KCF, fill = 'grey70', color = 'grey70') +
  geom_sf(data = outline, fill = NA, color = 'black', linewidth = 0.7) +
  geom_sf(data = plots, fill = 'white', color = 'black', size = 3.2, shape = 21, alpha = 0.9) +
  scale_fill_viridis(option = 'D', discrete = T, direction = -1) +
  theme(
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black')
  ) +
  labs(fill = 'Years\nungulate-free', x = '', y = '', title = '(a) Ungulate exclosures')

# Panel b: overstory cover type and the historic grazing and logging predictors
map_lus <- ggplot() +
  geom_sf(data = fence_u, fill = NA, color = 'grey75', linewidth = 0.3) +
  geom_sf(data = KCF, fill = 'grey70', color = 'grey70') +
  geom_sf(data = outline, fill = NA, color = 'black', linewidth = 0.7) +
  geom_sf(data = plots, aes(fill = cover_label, shape = land_use), color = 'black', size = 3.6) +
  scale_fill_viridis(option = 'C', discrete = T, direction = -1) +
  scale_shape_manual(values = c('Neither' = 21, 'Grazed' = 22, 'Logged' = 24, 'Grazed and logged' = 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black')
  ) +
  labs(fill = 'Overstory', shape = 'Land-use history', x = '', y = '', title = '(b) Land-use history')

pdf(file.path(fig_dir, 'disturbance_figure.pdf'), width = 14, height = 6, pointsize = 14)
NEON1::Multiplot(map_fence, map_lus, cols = 2)
dev.off()
png(file.path(fig_dir, 'disturbance_figure.png'), width = 14, height = 6, units = 'in', res = 300, pointsize = 14)
NEON1::Multiplot(map_fence, map_lus, cols = 2)
dev.off()
