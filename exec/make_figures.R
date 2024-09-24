# BEM

set.seed(1)

# TODO:
#      update this script to just do the figure making


# Libraries

# Load data
load(file.path(mod_dir, 'm_p_mod.rda'))
load(file.path(mod_dir, 'm_p_diag.rda'))

# Figure 1: Map figure

# Figure 3/4?: variance partitioning
NEON1::plot_vp(m_p_0, m_vp_0)

## Parameter means and highest posterior density intervals (HPD)
# only plot elevation, pai, age_median because those are standardized
pb_xx <- post_beta_0[which(post_beta_0$var %in% c('elevation', 'pai', 'age_median')), ]
pb_xx$var <- ifelse(pb_xx$var == 'age_median', 'Flow age', pb_xx$var)
pb_xx$var <- ifelse(pb_xx$var == 'elevation', 'Elevation', pb_xx$var)
pb_xx$var <- ifelse(pb_xx$var == 'pai', 'Plant area index', pb_xx$var)
pb_xx$spp <- gsub('_', ' ', pb_xx$spp)

NEON1::param_plot(post_beta, 'Parameter effect (Beta)')

