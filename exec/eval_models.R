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
library(ggplot2)
library(bayesplot)
library(MCMCvis)

# model parameters
s0 <- c(8000, 4000, 6, 6) # iterations, burn, chains, parallel
s1 <- (s0[1] - s0[2]) * s0[3] # for ES
# 10 dropped Hypochaeris_radicata, Hymenophyllum_recurvum, Axonopus_fissifolius' from 0, nothing from 1

#dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'
dir0 <- '..'
#dir0 <- '/media/bem/data/NEON'

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
#load(file.path(mod_dir, 'm_p_mod1.rda'))
#load(file.path(mod_dir, 'm_p_diag1.rda'))

### Evaluate

## Trace plots
# trace plot figures are in fig_dir

## Rhat
# Gill (2007) states failure to converge for one parameter is failure to converge for all parameters
# load in new rhat table
rhat <- read.csv(file.path(mod_dir, 'rhat_eval.csv'))
nrow(rhat[which(rhat$param == 'beta0'), ])
# 9 rows with bad Rhats for beta 0
nrow(rhat[which(rhat$param == 'gamm0'), ])
# no bad rhat for gamma 0
#nrow(rhat[which(rhat$param == 'omeg0'), ])
# 66 bad rows for omega 0 (but not present in rhat table for version 10)

# not worried about model 1 for now

## pp_check using bayesplot::
NEON1::pp_check(m_p = m_p_0, mp_pp = mp_pp_0, 2000)
# PP check looks good

# model fit diagnostics
# bind them all together with species names
mfs <- data.frame(spp = colnames(m_p_0$Y), RMSE = mf_p_0$RMSE, AUC = mf_p_0$AUC, TjurR2 = mf_p_0$TjurR2)
hist(mfs$RMSE, xlim = c(0,1), main = paste0("Mean = ", round(mean(mfs$RMSE), 2)), breaks = 30)
round(range(mfs$RMSE), 2)
# mean RMSE 0.24, range [0.08, 0.37], for 31 total species
# combine with species to show which are divergent
ggplot(mfs, aes(x = spp, y = RMSE)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# > 0.3 for Alyxia_stellata, Myrsine_lessertiana, Styphelia_tameiameiae, Uncinia_uncinata, Vaccinium_calycinum


# how much do traits explain environment-species relationships?
mve <- as.data.frame(m_vp_0$R2T)
mve <- data.frame(var = row.names(mve), mve, row.names = NULL)
# 22.02 % mean variance explained
# high of 41.45 for cowTRUE, with age_median at 11.30 and pai at 26.43

# check the HMSC book for other model fit diagnostics
