# BEM

set.seed(1)

# TODO:
#      update this script to just do the figure making


# Libraries

# Load data
mod_dir <- '../results/model_results'
fig_dir <- '../results/figures'

load(file.path(mod_dir, 'm_p_mod0.rda'))
load(file.path(mod_dir, 'm_p_diag0.rda'))
# use 95% HPD instead of 90%
#post_beta_0 <- NEON1::posterior_from_coda(mc_p_0, 'Beta', 0.95, average = T, drop_ns = T)

# Figures 1 and 2 (maps) are made by exec/make_map_figures.R
# Figure numbers below follow the manuscript (v0_6): 3 variance partitioning,
# 4-6 parameter plots, 7 trait variance

# Figure 4: Parameter means and HPD (i.e. effect sizes) for elevation, age_median, pai
# whisker width is NEON1::param_plot()'s err_width (0.63, the ggplot2 default reduced by 30%)
# these are 90% HPD calculated using the `coda` package
pb_xx <- post_beta_0[which(post_beta_0$var %in% c('elevation', 'pai', 'age_median')), ]
pb_xx$var <- ifelse(pb_xx$var == 'age_median', 'Flow age', pb_xx$var)
pb_xx$var <- ifelse(pb_xx$var == 'elevation', 'Elevation', pb_xx$var)
pb_xx$var <- ifelse(pb_xx$var == 'pai', 'Plant area index', pb_xx$var)
pb_xx$spp <- gsub('_', ' ', pb_xx$spp)
#NEON1::param_plot(pb_xx, 'Parameter effect (Beta)')
pdf(file.path(fig_dir, 'param_plot_env.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::param_plot(pb_xx, 'Parameter effect (Beta)')
dev.off()

# Figure 5: Parameter means and HPD (i.e. effect size) for time_since_fence, logTRUE, cowTRUE
pb_yy <- post_beta_0[which(post_beta_0$var %in% c('time_since_fence', 'logTRUE', 'cowTRUE')), ]
pb_yy$var <- ifelse(pb_yy$var == 'cowTRUE', 'Historic grazing', pb_yy$var)
pb_yy$var <- ifelse(pb_yy$var == 'logTRUE', 'Historic logging', pb_yy$var)
pb_yy$var <- ifelse(pb_yy$var == 'time_since_fence', 'Years ungulate-free', pb_yy$var)
pb_yy$spp <- gsub('_', ' ', pb_yy$spp)
#NEON1::param_plot(pb_yy, 'Parameter effect (Beta)')
pdf(file.path(fig_dir, 'param_plot_disturb.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::param_plot(pb_yy, 'Parameter effect (Beta)')
dev.off()

# Figure 6: Parameter means and HPD (i.e. effect size) for cover types
pb_zz <- post_beta_0[which(post_beta_0$var %in% grep('cover_type', post_beta_0$var, value = T)), ]
# drop species never recorded in a cover type: with no presences the probit
# likelihood is flat toward -Inf, so the coefficient's magnitude is set by the
# prior (BEM, MS v0_5 review)
ct_obs <- sapply(split(as.data.frame(m_p_0$Y), m_p_0$XData$cover_type), colSums)
ct_obs <- data.frame(spp = rep(row.names(ct_obs), ncol(ct_obs)), var = paste0('cover_type', rep(colnames(ct_obs), each = nrow(ct_obs))), n_pres = as.vector(ct_obs))
pb_zz <- dplyr::left_join(pb_zz, ct_obs, by = c('spp', 'var'))
cat('Figure 6: credible cover type terms', nrow(pb_zz), '; never observed in that cover type', sum(pb_zz$n_pres == 0), '; plotted', sum(pb_zz$n_pres > 0), '\n')
print(table(pb_zz$var, pb_zz$n_pres > 0))
pb_zz <- pb_zz[pb_zz$n_pres > 0, ]
pb_zz$var <- ifelse(pb_zz$var == 'cover_typekoa_tall', 'Mature koa overstory', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_ash_tall', 'Mature ohia-ash overstory ', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_koa_tall', 'Mature ohia-koa overstory', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_woodland', 'Ohia woodland', pb_zz$var)
pb_zz$spp <- gsub('_', ' ', pb_zz$spp)
#NEON1::param_plot(pb_zz, 'Parameter effect (Beta)')
# one row with panel widths proportional to the number of species, so a panel
# with a single species does not stretch its whisker across the panel
pdf(file.path(fig_dir, 'param_plot_cover.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::param_plot(pb_zz, 'Parameter effect (Beta)') + ggplot2::facet_grid(. ~ var, scales = 'free_x', space = 'free_x', labeller = ggplot2::as_labeller(c('Mature koa overstory' = 'Koa', 'Mature ohia-ash overstory ' = 'Ohia-ash', 'Mature ohia-koa overstory' = 'Ohia-koa', 'Ohia woodland' = 'Ohia woodland'))) + ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
dev.off()

# Figure 3: variance partitioning
pdf(file.path(fig_dir, 'vp_plot.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::plot_vp(m_p_0, m_vp_0)
dev.off()

# Figure 7: variance explained by traits (i.e. invasiveness), with 90% HPD
# intervals computed per posterior draw (NEON1::r2t_draws); the mean equals
# m_vp_0$R2T$Beta, which is what NEON1::mve stores
mve <- NEON1::r2t_draws(m_p_0, prob = 0.9)
stopifnot(isTRUE(all.equal(mve$mean, as.numeric(m_vp_0$R2T$Beta), tolerance = 1e-6)))
write.csv(mve, file.path(mod_dir, 'r2t_trait_hpd0.csv'), row.names = F)
mve <- mve[-which(mve$var %in% c('(Intercept)')), ]
mve$var <- ifelse(mve$var == 'age_median', 'Flow age', mve$var)
mve$var <- ifelse(mve$var == 'elevation', 'Elevation', mve$var)
mve$var <- ifelse(mve$var == 'pai', 'Plant area index', mve$var)
mve$var <- ifelse(mve$var == 'time_since_fence', 'Years ungulate-free', mve$var)
mve$var <- ifelse(mve$var == 'logTRUE', 'Historic logging', mve$var)
mve$var <- ifelse(mve$var == 'cowTRUE', 'Historic grazing', mve$var)
mve$var <- ifelse(mve$var == 'cover_typekoa_tall', 'Mature koa overstory', mve$var)
mve$var <- ifelse(mve$var == 'cover_typeohia_ash_tall', 'Mature ohia-ash overstory', mve$var)
mve$var <- ifelse(mve$var == 'cover_typeohia_koa_tall', 'Mature ohia-koa overstory', mve$var)
mve$var <- ifelse(mve$var == 'cover_typeohia_woodland', 'Ohia woodland', mve$var)
mve_plot <- ggplot2::ggplot(data = mve, ggplot2::aes(x = var, y = mean, ymin = lower, ymax = upper)) +
  ggplot2::geom_bar(stat = 'identity', fill = 'grey60') +
  ggplot2::geom_errorbar(width = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = 'black', angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(color = 'black')
  ) +
  ggplot2::labs(x = 'Covariate interaction', y = 'Variance explained by interaction with invasive trait')
pdf(file.path(fig_dir, 've_invasive.pdf'), width = 8, height = 6, pointsize = 14)
print(mve_plot)
dev.off()
