# BEM 29 July 2024

set.seed(1)

# TODO:
#      Get validation statistics for the two models:
#          Rhat, ESS, RMSE/R2, PSRF
#      Make plots for both models
#          Variance partitioning
#          Beta/Gamma parameter estimates with HPD
#          Species residual correlations

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
library(MCMCvis)

# model parameters
s0 <- c(13000, 3000, 6, 6) # iterations, burn, chains, parallel
s1 <- (s0[1] - s0[2]) * s0[3] # for ES

#dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'
dir0 <- '/media/bem/data/NEON'

### Directories
data_dir <- dir0
res_dir <- file.path(dir0, 'results')
mod_dir <- file.path(dir0, 'results/model_results')
fig_dir <- file.path(dir0, 'results/figures')
stopifnot(
  dir.exists(data_dir),
  dir.exists(res_dir),
  dir.exists(mod_dir),
  dir.exists(fig_dir),
  file.exists(file.path(mod_dir, 'm_p_mod.rda')),
  file.exists(file.path(mod_dir, 'm_p_diag.rda'))
)

### Load data for model 1
load(file.path(mod_dir, 'm_p_mod0.rda'))
load(file.path(mod_dir, 'm_p_diag0.rda'))

### Evaluate

## Cross-validation

## Rhat
# mc_s_beta_1, mc_s_beta_0, mc_s_gamma_1, etc.
# uses **$Rhat

## ESS
# mc_s_beta_1, mc_s_beta_0, mc_s_gamma_1, etc.
# uses **n.eff

## pp_check using bayesplot::
NEON1::pp_check(m_p = m_p, mp_pp = mp_pp, 2000)

## RMSE/R2
hist(mf_p$RMSE, xlim = c(0,1), main = paste0("Mean = ", round(mean(mf_p$RMSE), 2)), breaks = 30)
# RMSE for probit, R2 for others

## Gelman's PSRF (Potential Scale Reduction Factor)
# not sure what to do with this - report in table?
hist(gd_p_beta, breaks = 30, main = 'psrf:Beta, model p')
table(abs(gd_p_beta[, 1] - 1) > 0.01)
gd_p_beta[which(abs(gd_p_beta[, 1] - 1) > 0.01), ]
hist(gd_p_gamm, breaks = 30, main = 'psrf:Gamma, model p')
hist(gd_p_omeg, breaks = 30, main = 'psrf:Omega, model p')
hist(gd_p_omeg[which(abs(gd_p_omeg[, 1] - 1) > 0.01), 1], breaks = 30)
table(abs(gd_p_omeg[, 1] - 1) > 0.01)
gd_p_omeg[which(abs(gd_p_omeg[, 1] - 1) > 0.01), ]

### Generating data/eval for figures

## Variance partitioning

# figure: variance partitioning
NEON1::plot_vp(m_p, m_vp)
# how much do traits explain environment-species relationships?
knitr::kable(round(m_vp$R2T$Beta * 100, 2))
round(m_vp$R2T$Y * 100, 2)
# not much, 0.47%

## Parameter means and highest posterior density intervals (HPD)
# only plot elevation, pai, age_median because those are standardized
post_beta <- post_beta_0[which(post_beta$var %in% c('elevation', 'pai', 'age_median')), ]
post_beta$var <- ifelse(post_beta$var == 'age_median', 'Flow age', post_beta$var)
post_beta$var <- ifelse(post_beta$var == 'elevation', 'Elevation', post_beta$var)
post_beta$var <- ifelse(post_beta$var == 'pai', 'Plant area index', post_beta$var)
post_beta$spp <- gsub('_', ' ', post_beta$spp)

NEON1::param_plot(post_beta, 'Parameter effect (Beta)')
