# BEM 12 May 2023
# updated 18 Dec 2023

# Setup
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
#library(vegan)
library(NEON1)
library(gllvm)

set.seed(1)

# Filenames

# (1) How much of plant composition is explained by time or disturbance at Puʻu Makaʻala?

# Create matrix of cover data
data_div <- neonPlantEcology::npe_download(sites = c("PUUM"))
# can add WREF and GUAN if need be
data_df <- neonPlantEcology::npe_longform(data_div, scale = 'plot')
data_df$cover <- data_df$cover / 100
data_df <- data_df[-which(data_df$taxonID == '2PLANT'), ]
colnames(data_df)[which(colnames(data_df) == 'eventID')] <- 'year'
data_mat <- neonPlantEcology::npe_community_matrix(data_div)
data_mat <- data_mat[, -which(colnames(data_mat) == '2PLANT')]
#data_mat <- data_mat[grepl('2018', row.names(data_mat)), ]
#data_mat <- data_mat[rowSums(data_mat) > 0, ]
data_mat <- data_mat / 100

# single zero-inflation parameter applied to all observations
# ziformula ~ 1
# hurdle model: ziformula ~ . ???

c0 <- glmmTMBControl(
  optCtrl = list(iter.max = 1e3, eval.max = 1e3),
  start_method = list(method = NULL, jitter.sd = 0)
)

# this is glmmTMB's version of a GLVM using a reduced rank covariance structure
# from: https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
mod0 <- glmmTMB(cover ~ age_group + year,
                data = data_df, family = beta_family(), ziformula = ~taxonID,
                control = c0)
# mod0 zi parameter correlates with plot occurence across all taxon, but
# still some poor fit
mod1 <- glmmTMB(cover ~ age_group + year + rr(taxonID + 0|plotID, d = 2),
                data = data_df, family = beta_family(), ziformula = ~taxonID,
                control = c0)

mod0_sim <- DHARMa::simulateResiduals(mod0)
plot(mod0_sim)

fit.rr <- mod0
ll <- fit.rr$obj$env$report(fit.rr$fit$parfull)$fact_load[[1]]
ll <- as.data.frame(ll)
ll <- setNames(ll, c('L1', 'L2'))


