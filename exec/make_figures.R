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

# Figure 1: Map figure

# Figure 2: Second map figure, showing disturbances

# Figure 4: Parameter means and HPD (i.e. effect sizes) for elevation, age_median, pai
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
pb_zz$var <- ifelse(pb_zz$var == 'cover_typekoa_tall', 'Mature koa overstory', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_ash_tall', 'Mature ohia-ash overstory ', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_koa_tall', 'Mature ohia-koa overstory', pb_zz$var)
pb_zz$var <- ifelse(pb_zz$var == 'cover_typeohia_woodland', 'Ohia woodland', pb_zz$var)
pb_zz$spp <- gsub('_', ' ', pb_zz$spp)
#NEON1::param_plot(pb_zz, 'Parameter effect (Beta)')
pdf(file.path(fig_dir, 'param_plot_cover.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::param_plot(pb_zz, 'Parameter effect (Beta)')
dev.off()

# Figure 7: variance partitioning
NEON1::plot_vp(m_p_0, m_vp_0)
pdf(file.path(fig_dir, 'vp_plot.pdf'), width = 8, height = 6, pointsize = 14)
NEON1::plot_vp(m_p_0, m_vp_0)
dev.off()

# Figure 8: variance explained by traits (i.e. invasiveness)
mve <- NEON1::mve
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
mve_plot <- ggplot(data = mve, aes(x = var, y = Beta)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),
    axis.text.y = element_text(color = 'black')
  ) +
  labs(x = 'Covariate interaction', y = 'Variance explained by interaction with invasive trait')
mve_plot
pdf(file.path(fig_dir, 've_invasive.pdf'), width = 8, height = 6, pointsize = 14)
mve_plot
dev.off()

