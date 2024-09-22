# BEM 29 July 2024

set.seed(1)

# TODO:
#      Get validation statistics for the two models:
#          Rhat, ESS, RMSE/R2, PSRF
#      Make plots for both models
#          Variance partitioning
#          Beta/Gamma parameter estimates with HPD
#          Species residual correlations

# Motivation:
#
# Model diagnostics fall into two categories: model performance in a technical sense, and model performance in a scientific sense
#
# TECHNICAL: Rhat, ESS (in part), pp_check, PSRF
# SCIENTIFIC: ESS (in part), RMSE, xv-RMSE, variance partitioning

### Libraries
library(NEON1)
library(dplyr)
library(Hmsc)
library(coda)
library(bayesplot)
library(MCMCvis)

# model parameters
s0 <- c(6000, 3000, 4, 4) # iterations, burn, chains, parallel
s1 <- (s0[1] - s0[2]) * s0[3] # for ES

#dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'
#dir0 <- '..'
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
  file.exists(file.path(mod_dir, 'm_p_mod0.rda')),
  file.exists(file.path(mod_dir, 'm_p_diag0.rda')),
  file.exists(file.path(mod_dir, 'm_p_mod1.rda')),
  file.exists(file.path(mod_dir, 'm_p_diag1.rda'))
)

### Load data for models
load(file.path(mod_dir, 'm_p_mod0.rda'))
load(file.path(mod_dir, 'm_p_diag0.rda'))
load(file.path(mod_dir, 'm_p_mod1.rda'))
load(file.path(mod_dir, 'm_p_diag1.rda'))

### Evaluate

## Trace plots
# trace plot figures are in fig_dir

## Rhat
# Gill (2007) states failure to converge for one parameter is failure to converge for all parameters
rhat_beta_0 <- mc_s_beta_0[which(mc_s_beta_0$Rhat > 1.01), ]
rhat_beta_0
table(sapply(strsplit(row.names(rhat_beta_0), ' '), \(xx) xx[3]))
sum(table(sapply(strsplit(row.names(rhat_beta_0), ' '), \(xx) xx[3])))
rhat_gamm_0 <- mc_s_gamm_0[which(mc_s_gamm_0$Rhat > 1.01), ]
rhat_gamm_0
# no bad rhat
rhat_omeg_0 <- mc_s_omeg_0[which(mc_s_omeg_0$Rhat > 1.01), ]
rhat_omeg_0
# Adenophorus_tamariscinus, Coprosma_ernodeoides, Stenogyne_calaminthoides minor offendors
# Trichomanes_bauerianum, Myrsine_lessertiana major offenders
t0 <- table(sapply(strsplit(row.names(rhat_omeg_0), ' '), \(xx) xx[1]))
names(t0) <- gsub('Omega1\\[', '', names(t0))
t0 <- colSums(dplyr::bind_rows(t0, table(sapply(strsplit(row.names(rhat_omeg_0), ' '), \(xx) xx[3]))), na.rm = T)
sum(t0)

rhat_beta_1 <- mc_s_beta_1[which(mc_s_beta_1$Rhat > 1.05), ]
rhat_beta_1
table(sapply(strsplit(row.names(rhat_beta_1), ' '), \(xx) xx[3]))
sum(table(sapply(strsplit(row.names(rhat_beta_1), ' '), \(xx) xx[3])))
# 15 bad parameters, mostly cover_typeohia_woodland
# Elaphoglossum_alatum, Freycinetia_arborea, Pneumatopteris_sandwicensis x2 each
# Adenophorus_tamariscinus, Hymenophyllum_lanceolatum, Hymenophyllum_recurvum, Metrosideros_polymorpha,
# Psidium_cattleianum, Psilotum_complanatum, Stenogyne_calaminthoides, Trichomanes_bauerianum x1
# but all Rhats are = 1.02 and n.eff is relatively large
rhat_gamm_1 <- mc_s_gamm_1[which(mc_s_gamm_1$Rhat > 1.01), ]
rhat_gamm_1
# no bad parameters
rhat_omeg_1 <- mc_s_omeg_1[which(mc_s_omeg_1$Rhat > 1.01), ]
rhat_omeg_1
t1 <- table(sapply(strsplit(row.names(rhat_omeg_1), ' '), \(xx) xx[1]))
names(t1) <- gsub('Omega1\\[', '', names(t1))
t1 <- colSums(dplyr::bind_rows(t1, table(sapply(strsplit(row.names(rhat_omeg_1), ' '), \(xx) xx[3]))), na.rm = T)
sum(t1)
# 13 bad params
# Microlaena_stipoides is only potential issue, although some Rhats are large

## ESS
# mc_s_beta_1, mc_s_beta_0, mc_s_gamma_1, etc.
# uses **n.eff

## pp_check using bayesplot::
NEON1::pp_check(m_p = m_p, mp_pp = mp_pp, 2000)

## RMSE/R2
hist(mf_p$RMSE, xlim = c(0,1), main = paste0("Mean = ", round(mean(mf_p$RMSE), 2)), breaks = 30)
# RMSE for probit, R2 for others

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
