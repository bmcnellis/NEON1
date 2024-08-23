# BEM 29 July 2024

set.seed(1)

# TODO:
#      implement cross-validation

# Caveats/Considerations:
#     * taxon not identified to species are excluded
#     * two enclosures (Kata-Stein & Unencumbered Kulani Pasture) had NA for completion date, this
#       was set to 2014 which is the age of the youngest bordering enclosure
#     * lava flow age is the median of the age range published on USGS maps
#     * pai, elev, age_median are centered and scaled
#     * time_since_fence NOT centered or scaled
#     * one enclosure (Army Road) was listed as NOT ungulate-free, time_since_fence set to 0

# "All variance partitioning values for each species were multiplied by the explanatory R2 value of that species to show amount
# of total variation in the response variable explained by each covariate."

### Libraries
library(NEON1)
library(dplyr)
library(Hmsc)
library(coda)
library(bayesplot)

dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'

### Directories
data_dir <- dir0
res_dir <- file.path(dir0, 'results')
mod_dir <- file.path(dir0, 'results/model_results')
fig_dir <- file.path(dir0, '/results/figures')
stopifnot(
  dir.exists(data_dir),
  dir.exists(res_dir),
  dir.exists(mod_dir),
  dir.exists(fig_dir),
  file.exists(file.path(mod_dir, 'm_p_mod.rda')),
  file.exists(file.path(mod_dir, 'm_p_diag.rda'))
)

### Load data
load(file.path(mod_dir, 'm_p_mod.rda'))
load(file.path(mod_dir, 'm_p_diag.rda'))

### Evaluate

## Cross-validation

## Posterior predictive checks
if (F) { # summarize the posterior predictive distribution into posterior mean and then extract and standardize the residuals.
  pm <- apply(mp_p, FUN = mean, MARGIN = 1)
  # y is response variable
  nres <- scale(y - pm)
  par(mfrow = c(1,2))
  hist(nres, las = 1)
  plot(pm, nres, las = 1)
  abline(a = 0, b = 0)
}
## pp_check using bayesplot::
if (F) {
  n_draws <- 1000 # the 12000 is way too many
  m0 <- matrix(sapply(mp_pp, as.numeric), ncol = length(as.numeric(m_p$Y)), nrow = length(mp_pp))
  m0 <- m0[sample.int(nrow(m0), n_draws), ]
  bayesplot::pp_check(as.numeric(m_p$Y), m0, bayesplot::ppc_dens_overlay)
}

## Effective sample size
hist(es_p_beta, breaks = 30, main = 'ess:Beta, model p')
hist(es_p_gamm, breaks = 30, main = 'ess:Gamma, model p')
hist(es_p_omeg, breaks = 30, main = 'ess:Omega, model p')
# which species have low ess?
es_p_beta <- NEON1::ess_as_df(es_p_beta, 'beta')
es_p_beta_spp <- summarise(group_by(es_p_beta, spp), ess_mean = mean(ess), ess_sd = sd(ess)/n())
es_p_beta_param <- summarise(group_by(es_p_beta, param), ess_mean = mean(ess), ess_sd = sd(ess)/n())
NEON1::ess_plot(es_p_beta_spp, spp, 'ESS (Beta)')
NEON1::ess_plot(es_p_beta_param, param, 'ESS (Beta)')

## RMSE/R2
hist(mf_p$RMSE, xlim = c(0,1), main = paste0("Mean = ", round(mean(mf_p$RMSE), 2)), breaks = 30)
# RMSE for probit, R2 for others

## Gelman's PSRF (Potential Scale Reduction Factor)
# Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations.
# Journal of Computational and Graphical Statistics, 7, 434-455.
hist(gd_p_beta, breaks = 30, main = 'psrf:Beta, model p')
table(abs(gd_p_beta[, 1] - 1) > 0.02)
hist(gd_p_gamm, breaks = 30, main = 'psrf:Gamma, model p')
hist(gd_p_omeg, breaks = 30, main = 'psrf:Omega, model p')
hist(gd_p_omeg[which(abs(gd_p_omeg[, 1] - 1) > 0.05), 1], breaks = 30)
table(abs(gd_p_omeg[, 1] - 1) > 0.05)

### Generating data/eval for figures

## Variance partitioning

# figure: variance partitioning
NEON1::plot_vp(m_p, m_vp)
# how much do traits explain environment-species relationships?
knitr::kable(round(m_vp$R2T$Beta * 100, 2))
round(m_vp$R2T$Y * 100, 2)
# not much, 0.47%

## Parameter means and highest posterior density intervals (HPD)
post_beta <- NEON1::posterior_from_coda(mc_p, 'Beta', 0.9, average = T, drop_ns = T)
# only plot elevation, pai, age_median because those are standardized
post_beta <- post_beta[which(post_beta$var %in% c('elevation', 'pai', 'age_median')), ]
post_beta$var <- ifelse(post_beta$var == 'age_median', 'Flow age', post_beta$var)
post_beta$var <- ifelse(post_beta$var == 'elevation', 'Elevation', post_beta$var)
post_beta$var <- ifelse(post_beta$var == 'pai', 'Plant area index', post_beta$var)
post_beta$spp <- gsub('_', ' ', post_beta$spp)

NEON1::param_plot(post_beta, 'Parameter effect (Beta)')

### old crap

plotBeta(m_p, post = ml_pb, param = 'Support', supportLevel = 0.90)
plotBeta(m_p, post = ml_pb, param = 'Mean', supportLevel = 0.90)
plotGamma(m_p, post = ml_pg, param = 'Support', supportLevel = 0.90)
knitr::kable(m_vp$R2T$Beta)

# i think this is now ma_p
sl <- 0.5
cp <- ((m_ca[[1]]$support > sl) + (m_ca[[1]]$support < (1 - sl)) > 0) * m_ca[[1]]$mean
corrplot::corrplot(
  cp, method = 'color', col = colorRampPalette(c('blue', 'white', 'red'))(200),
  title = paste('random effect level:', m_p$rLNames[1]), mar = c(0, 0, 1, 0)
)

