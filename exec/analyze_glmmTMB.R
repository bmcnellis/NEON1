# BEM 12 May 2023
# updated 18 Dec 2023

# Setup
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
library(NEON1)
library(glmmTMB)

set.seed(1)

data_dir <- '/media/bem/data/NEON'

# Cover data
#data_div <- neonPlantEcology::npe_download(sites = c("WREF"))
data_div <- NEON1::PUUM_div
data_df <- neonPlantEcology::npe_longform(data_div, scale = 'plot')
data_df$cover <- data_df$cover / 100
data_df <- data_df[-which(data_df$taxonID == '2PLANT'), ]
colnames(data_df)[which(colnames(data_df) == 'eventID')] <- 'year'

# Other data
data_df <- dplyr::left_join(data_df, NEON1::flow_meta, by = 'plotID')
hemi <- neonUtilities::stackByTable(filepath = file.path(data_dir, 'NEON_hemispheric-photos-veg.zip'), savepath = 'envt')
traits <- neonUtilities::stackByTable(filepath = file.path(data_dir, 'NEON_traits-foliar.zip'), savepath = 'envt')

# Create data matrix
data_mat <- neonPlantEcology::npe_community_matrix(data_div)
data_mat <- data_mat[, -which(colnames(data_mat) == '2PLANT')]
data_mat <- data_mat / 100

c0 <- glmmTMBControl(
  optCtrl = list(iter.max = 1e3, eval.max = 1e3),
  start_method = list(method = NULL, jitter.sd = 0)
)

# this is glmmTMB's version of a GLVM using a reduced rank covariance structure
# from: https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
fit <- glmmTMB(cover ~ age_group + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                data = data_df, family = beta_family(), ziformula = ~taxonID,
                control = c0, start = list())

fit_sim <- DHARMa::simulateResiduals(fit)
plot(fit_sim)

ll <- fit$obj$env$report()$fact_load[[2]] |>
  as.data.frame() |>
  cbind(unique(data_df$taxonID)) |>
  cbind(unique(data_df$taxonID)[order(unique(data_df$taxonID))]) |>
  setNames(c("L1", "L2", "taxonID", "taxonID_ordered"))

results_LV_df <- ll
usethis::use_data(results_LV_df, overwrite = T)

