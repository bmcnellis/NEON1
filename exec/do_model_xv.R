# BEM 29 July 2024

set.seed(1)

### Libraries
library(NEON1)
library(Hmsc)

# model parameters
s0 <- c(6000, 2000, 4, 4)# iterations, burn, chains, parallel

dir0 <- '/media/bem/data/NEON'
#dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'

### Directories
mod_dir <- file.path(dir0, 'results/model_results')
stopifnot(dir.exists(mod_dir))

# cross-validation
nf <- c(2, 5, 10)

load(file.path(mod_dir, 'm_p_mod0.rda'))
rm(list = ls()[-which(ls() %in% c('m_p_0', 'mod_dir', 'nf', 's0'))])
gc()

for (i in seq_along(nf)) {

  ip_0 <- Hmsc::createPartition(m_p_0, nfolds = nf[i])
  ipv_0 <- Hmsc::computePredictedValues(m_p_0, partition = ip_0, nParallel = s0[4])
  ipmf_0 <- Hmsc::evaluateModelFit(hM = m_p_0, predY = ipv_0)

  save(list = c('ipv_0', 'ipmf_0'), file = file.path(mod_dir, paste0('m_p_xv0_', i, '.rda')))

  rm(ip_0, ipv_0, ipmf_0)
  gc()
}

load(file.path(mod_dir, 'm_p_mod1.rda'))
rm(list = ls()[-which(ls() %in% c('m_p_1', 'mod_dir', 'nf', 's0'))])

for (i in seq_along(nf)) {
  ip_1 <- Hmsc::createPartition(m_p_1, nfolds = nf[i])
  ipv_1 <- Hmsc::computePredictedValues(m_p_1, partition = ip_1, nParallel = s0[4])
  ipmf_1 <- Hmsc::evaluateModelFit(hM = m_p_1, predY = ipv_1)

  save(list = c('ipv_1', 'ipmf_1'), file = file.path(mod_dir, paste0('m_p_xv1_', i, '.rda')))

  rm(ip_1, ipv_1, ipmf_1)
  gc()
}
