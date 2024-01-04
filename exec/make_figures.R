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
