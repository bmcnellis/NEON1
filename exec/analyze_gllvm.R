# BEM 12 May 2023
# updated 18 Dec 2023

# Setup
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
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

df_X <- do.call('rbind', strsplit(row.names(data_mat), '_plot_'))
df_X <- setNames(as.data.frame(df_X), c('plotID', 'year'))
df_X <- dplyr::left_join(df_X, NEON1::flow_meta[, c('plotID', 'age_group')])
data_df <- dplyr::left_join(data_df, NEON1::flow_meta[, c('plotID', 'age_group')])
yr <- df_X$year
#df_X <- within(df_X, rm(year))

mod0 <- gllvm::gllvm(
  y = data_mat, X = df_X,
  family = 'beta', num.lv = 2L, link = 'logit', starting.val = 'zero'
)
#  # probably need to set a row effect by year
#  studyDesign = data.frame(year = yr), row.eff = ~(1|year)
#)
