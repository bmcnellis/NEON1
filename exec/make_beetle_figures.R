# BEM Nov 2024

# Purpose: make an NMDS with the PUUM beetle data

library(NEON1)
library(vegan)
library(dplyr)

bl_reo <- read.csv('../results/tables/bet_counts_simple.csv')

elev <- bl_reo |>
  select(plotID, elevation) |>
  distinct() |>
  group_by(plotID) |>
  mutate(elevation = mean(elevation))

bl_grp <- bl_reo |>
  select(plotID, trapID, collectDate, individualCount, taxonID) |>
  group_by(plotID, collectDate, taxonID) |>
  summarize(individualCount = sum(individualCount), .groups = 'drop') |>
  filter(!is.na(taxonID)) |>
  left_join(elev, by = c('plotID')) |>
  mutate(date = as.POSIXlt(as.Date(collectDate))) |>
  mutate(year = date$year + 1900) |>
  mutate(year_fac = as.factor(year)) |>
  mutate(mon = sprintf('%02d', date$mon + 1)) |>
  mutate(mon_fac = as.factor(mon)) |>
  mutate(YM = paste0(as.character(year), as.character(mon))) |>
  mutate(YM_fac = as.factor(YM)) |>
  mutate(sample = paste(plotID, collectDate, sep = '_')) |>
  as.data.frame()

bl_mat <- bl_grp |>
  select(sample, taxonID, individualCount) |>
  tidyr::pivot_wider(names_from = 'taxonID', values_from = 'individualCount') |>
  tibble::column_to_rownames('sample') |>
  # drop samples with < 30 occurences
  select(-AGOMUE, -MECPUN1, -MECNEO, -MECNIT, -MECMAU, -MECVUL2, -MECGAG)
bl_mat <- bl_mat[rowSums(!is.na(bl_mat)) != 0, ]
bl_mat <- bl_mat |>
  as.matrix() |>
  apply(c(1, 2), \(xx) ifelse(is.na(xx), 0, xx))

bl_grp <- bl_grp[match(row.names(bl_mat), bl_grp$sample), ]
identical(bl_grp$sample, row.names(bl_mat))
grp_year <- bl_grp[, c('year')]
grp_plot <- bl_grp[, c('plotID')]
grp_YM <- bl_grp[, c('YM_fac')]

bl_MDS <- vegan::metaMDS(bl_mat, k = 3, trymax = 100)

m_s <- 0.0 # small, reduced margin size
m_x <- 0.5 # x axis margin size when labels are present
m_y <- 0.7 # y axis margin size when labels are present

c_y_p <- scales::alpha(viridis::magma((5 * 5))[c(rep(F, 4), T)], 0.8)
c_y_g <- scales::alpha(viridis::magma((8 * 5))[c(rep(F, 4), T)], 0.8)
c_y_a <- scales::alpha(viridis::magma((7 * 5))[c(rep(F, 4), T)], 0.8)

p_lty <- 2 # line type for plot panels
l_cex <- 1.0 # legend size scaler for color legend
lb_cex <- 1.6 # legend size scaler for site/species legend
c_p <- scales::alpha("grey50", 0.6)

# PUUM - Year
ordiplot(bl_MDS, type = 'none')
points(bl_MDS, display = "species", pch = 19, cex = 0.8, col = c_p)
points(bl_MDS, display = "sites", pch = 19, cex = 0.7, col = 'black')
#ordihull(bl_MDS, groups = grp_year, draw = 'lines', col = c_y_p, lty = 6, lwd = 2.5)
ordihull(bl_MDS, groups = grp_year, draw = 'lines', lwd = 2.5)
legend(x = "bottomright", legend = levels(bl_grp$year_fac), title = "Year", fill = c_y_p, cex = l_cex)

# PUUM - Plot
ordiplot(bl_MDS, type = 'none')
points(bl_MDS, display = "species", pch = 19, cex = 0.8, col = c_p)
points(bl_MDS, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(bl_MDS, groups = grp_plot, draw = 'lines', col = 'black', lwd = 2)
legend(x = "topright", legend = c('Plots', 'Species'), fill = c('black', 'grey50'), cex = lb_cex)

# PUUM - YM
ordiplot(bl_MDS, type = 'none')
points(bl_MDS, display = "species", pch = 19, cex = 0.8, col = c_p)
points(bl_MDS, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(bl_MDS, groups = grp_YM, draw = 'lines', col = 'black', lty = p_lty, lwd = 2)
legend(x = "topright", legend = c('Plots', 'Species'), fill = c('black', 'grey50'), cex = lb_cex)
