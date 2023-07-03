# BEM 12 May 2023

#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
library(NEON1)
library(vegan)
library(pairwiseAdonis)

set.seed(1)

plot_file <- "C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/NMDS_multipanel.tiff"
beet_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_count-beetles.zip'
beet_file <- "C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/NMDS_beet.tiff"

beet_tbl <- neonUtilities::stackByTable(beet_zip, savepath = 'envt')
data_div <- neonPlantEcology::download_plant_div(sites = c("PUUM", "ABBY", "GUAN"))
data_df <- neonPlantEcology::get_longform_cover(data_div, scale = 'plot')
data_mat <- neonPlantEcology::get_community_matrix(data_div)

# PUUM, just beetles
beet_p <- beet_tbl[['bet_parataxonomistID']] %>%
  #filter(taxonRank == 'species') %>%
  select('siteID', 'plotID', 'trapID', 'collectDate', 'taxonID', 'scientificName') %>%
  # Mecyclothorax = mostly Hawaiian
  # Trechus = nonnative
  # Blackburnia = native
  # Agonum = nonnative
  mutate(genus = sapply(strsplit(scientificName, ' '), \(xx) xx[1])) %>%
  mutate(native = ifelse(genus %in% c('Mecyclothorax', 'Blackburnia'), 'N', 'I')) %>%
  #mutate(mm_yyyy = paste(lubridate::month(collectDate), lubridate::year(collectDate), sep = '_')) %>%
  mutate(mm_yyyy = collectDate) %>%
  mutate(plot_date = paste(plotID, mm_yyyy, sep = '_')) %>%
  # PUUM_013_2020-06-23 is insane outlier??
  filter(plot_date != 'PUUM_013_2020-06-23')

# transform to community matrix
beet_mat <- beet_p %>%
  mutate(count = rep(1, nrow(.))) %>%
  select('plot_date', 'taxonID', 'count') %>%
  tidyr::pivot_wider(names_from = taxonID, values_from = count, values_fn = sum, values_fill = 0) %>%
  as.data.frame() %>%
  `rownames<-`(.[, 'plot_date']) %>%
  select(-plot_date)

#NMDS_screeplot(beet_mat, k = 15)
beet_mds <- vegan::metaMDS(beet_mat, distance = 'bray', k = 3, try = 50, trymax = 300)
# best solution wasnt repeated
beet_mds$stress # 0.039
# for plotting/stats:
beet_g_y <- as.factor(lubridate::year(strptime(sapply(strsplit(row.names(beet_mat), '_'), \(xx) xx[3]), '%Y-%m-%d')))
beet_g_p <- as.factor(sapply(strsplit(row.names(beet_mat), '_'), \(xx) xx[2]))
adonis2(beet_mat ~ beet_g_p + beet_g_y, by = 'margin')
# both significant
pairwise.adonis(beet_mat, factors = beet_g_y)
# all years but 19 vs 20 different from one another

# PUUM
PUUM_mat <- data_mat[grepl('PUUM', row.names(data_mat)), ]
PUUM_mat <- PUUM_mat[, colSums(PUUM_mat) > 0]
#NMDS_screeplot(PUUM_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions, can try to get away with 3
m_PUUM <- metaMDS(PUUM_mat, k = 3, try = 50, trymax = 200)
m_PUUM$stress # 0.1367
# for plotting/stats:
PUUM_g_y <- as.factor(sapply(strsplit(row.names(PUUM_mat), '_'), \(xx) xx[4]))
PUUM_g_p <- as.factor(sapply(strsplit(row.names(PUUM_mat), '_'), \(xx) xx[2]))
# significance:
adonis2(PUUM_mat ~ PUUM_g_p + PUUM_g_y, by = 'margin')
pairwise.adonis(PUUM_mat, factors = PUUM_g_y)

# GUAN
GUAN_mat <- data_mat[grepl('GUAN', row.names(data_mat)), ]
GUAN_mat <- GUAN_mat[, colSums(GUAN_mat) > 0]
#NMDS_screeplot(GUAN_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions, but still high
m_GUAN <- metaMDS(GUAN_mat, k = 4, try = 50, trymax = 200)
m_GUAN$stress # 0.1661
# for plotting:
GUAN_g_y <- as.factor(sapply(strsplit(row.names(GUAN_mat), '_'), \(xx) xx[4]))
GUAN_g_p <- as.factor(sapply(strsplit(row.names(GUAN_mat), '_'), \(xx) xx[2]))
# significance:
adonis2(GUAN_mat ~ GUAN_g_p + GUAN_g_y, by = 'margin')
pairwise.adonis(GUAN_mat, factors = GUAN_g_y)

# ABBY
ABBY_mat <- data_mat[grepl('ABBY', row.names(data_mat)), ]
ABBY_mat <- ABBY_mat[, colSums(ABBY_mat) > 0]
#NMDS_screeplot(ABBY_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions, but still high
m_ABBY <- metaMDS(ABBY_mat, k = 3, try = 50, trymax = 200)
m_ABBY$stress # 0.1638
# for plotting:
ABBY_g_y <- as.factor(sapply(strsplit(row.names(ABBY_mat), '_'), \(xx) xx[4]))
ABBY_g_p <- as.factor(sapply(strsplit(row.names(ABBY_mat), '_'), \(xx) xx[2]))
# significance:
adonis2(ABBY_mat ~ ABBY_g_p + ABBY_g_y, by = 'margin')
pairwise.adonis(ABBY_mat, factors = ABBY_g_y)

# Make 6-panel NMDS plot

#c_y <- scales::alpha(RColorBrewer::brewer.pal(9, "YlOrBr")[c(9, 8, 7, 5, 3)], 0.8)
c_y_p <- scales::alpha(viridis::magma((5 * 5))[c(rep(F, 4), T)], 0.8)
c_y_g <- scales::alpha(viridis::magma((8 * 5))[c(rep(F, 4), T)], 0.8)
c_y_a <- scales::alpha(viridis::magma((7 * 5))[c(rep(F, 4), T)], 0.8)
c_p <- scales::alpha("grey50", 0.6)

m_s <- 0.0 # small, reduced margin size
m_x <- 0.5 # x axis margin size when labels are present
m_y <- 0.7 # y axis margin size when labels are present

p_lty <- 2 # line type for plot panels
l_cex <- 1.0 # legend size scaler for color legend
lb_cex <- 1.6 # legend size scaler for site/species legend

y_1 <- "PUUM\nNMDS2" # Y axis label for top row
y_2 <- "GUAN\nNMDS2" # Y axis label for mid row
y_3 <- "ABBY\nNMDS2" # Y axis label for bot row

tiff(plot_file, width = 8.5, height = 11, units = 'in', res = 300)
par(mfrow = c(3, 2), mgp = c(0.5, 1, 0))

# PUUM - Year
par(mai = c(m_s, m_y, m_s, m_s))
ordiplot(m_PUUM, type = 'none', xaxt = 'n', yaxt = 'n', xlab = '', ylab = y_1)
points(m_PUUM, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_PUUM, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_PUUM, groups = PUUM_g_y, draw = 'lines', col = c_y_p, lty = 6, lwd = 2.5)
legend(x = "bottomright", legend = levels(PUUM_g_y), title = "Year", fill = c_y_p, cex = l_cex)

# PUUM - Plot
par(mai = c(m_s, m_s, m_s, m_s))
ordiplot(m_PUUM, type = 'none', xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
points(m_PUUM, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_PUUM, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_PUUM, groups = PUUM_g_p, draw = 'lines', col = 'black', lty = p_lty, lwd = 2)
legend(x = "bottomright", legend = c('Plots', 'Species'), fill = c('black', 'grey50'), cex = lb_cex)

# GUAN - Year
par(mai = c(m_s, m_y, m_s, m_s))
ordiplot(m_GUAN, type = 'none', xaxt = 'n', yaxt = 'n', xlab = '', ylab = y_2)
points(m_GUAN, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_GUAN, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_GUAN, groups = GUAN_g_y, draw = 'lines', col = c_y_g, lty = 6, lwd = 2.5)
legend(x = "bottomright", legend = levels(GUAN_g_y), title = "Year", fill = c_y_g, cex = l_cex)

# GUAN - Plot
par(mai = c(m_s, m_s, m_s, m_s))
ordiplot(m_GUAN, type = 'none', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
points(m_GUAN, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_GUAN, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_GUAN, groups = GUAN_g_p, draw = 'lines', col = 'black', lty = p_lty, lwd = 2)

# ABBY - Year
par(mai = c(m_x, m_y, m_s, m_s))
ordiplot(m_ABBY, type = 'none', xaxt = 'n', yaxt = 'n', ylab = y_3)
points(m_ABBY, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_ABBY, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_ABBY, groups = ABBY_g_y, draw = 'lines', col = c_y_a, lty = 6, lwd = 2.5)
legend(x = "bottomright", legend = levels(ABBY_g_y), title = "Year", fill = c_y_a, cex = l_cex)

# ABBY - Plot
par(mai = c(m_x, m_s, m_s, m_s))
ordiplot(m_ABBY, type = 'none', xaxt = 'n', yaxt = 'n', ylab = '')
points(m_ABBY, display = "species", pch = 19, cex = 0.8, col = c_p)
points(m_ABBY, display = "sites", pch = 19, cex = 0.7, col = 'black')
ordihull(m_ABBY, groups = ABBY_g_p, draw = 'lines', col = 'black', lty = p_lty, lwd = 2)

dev.off()

# PUUM - Beetles
b_cex <- 0.7
b2_cex <- 0.7

tiff(beet_file, width = 9, height = 6, units = 'in', res = 300)
par(mfrow = c(1, 2), mgp = c(0.5, 1, 0))

# Year
par(mai = c(m_x, m_y, m_s, m_s))
ordiplot(beet_mds, type = 'none', xaxt = 'n', yaxt = 'n')
points(beet_mds, display = "sites", pch = 19, cex = 0.7, col = 'black')
points(beet_mds, display = "species", pch = 19, cex = 0.8, col = c_p)
ordihull(beet_mds, groups = beet_g_y, draw = 'lines', col = c_y_a, lty = 6, lwd = 2.5)
legend(x = "bottomright", legend = levels(beet_g_y), title = "Year", fill = c_y_a, cex = b_cex)
# Plot
par(mai = c(m_x, m_s, m_s, m_s))
ordiplot(beet_mds, type = 'none', xaxt = 'n', yaxt = 'n')
points(beet_mds, display = "sites", pch = 19, cex = 0.7, col = 'black')
points(beet_mds, display = "species", pch = 19, cex = 0.8, col = c_p)
ordihull(beet_mds, groups = beet_g_p, draw = 'lines', col = c_y_a, lty = 6, lwd = 2.5)
legend(x = "bottomright", legend = c('Plots', 'Species'), fill = c('black', 'grey50'), cex = b2_cex)

dev.off()
