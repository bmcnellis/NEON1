# BEM 12 May 2023
# updated 6 Feb 2023

# Caveats:
#     soil age is an ordered factor
#     'Rhynchospora chinensis' assumed native, b/c native ssp, unsure on other ssp presence

# Setup
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
library(NEON1)
library(glmmTMB)
library(DHARMa)

set.seed(1)

#data_dir <- '/media/bem/data/NEON'
data_dir <- 'C:/Users/BrandonMcNellis/Documents/NEON_data'
resid_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/reports'

# Cover data
#PUUM_div <- neonPlantEcology::npe_download(sites = c("PUUM"))
#usethis::use_data(PUUM_div, overwrite = T)  # "2024-05-13 12:32:03 HST"
data_div <- NEON1::PUUM_div
data_df <- neonPlantEcology::npe_longform(data_div, scale = 'plot')
data_df$cover <- data_df$cover / 100
data_df <- data_df[-which(data_df$taxonID == '2PLANT'), ]
colnames(data_df)[which(colnames(data_df) == 'eventID')] <- 'year'

# Other data/data modification
data_df <- dplyr::left_join(data_df, NEON1::flow_meta, by = 'plotID')
#data_df$age_fac <- as.factor(data_df$age_group - 1)
data_df$age_fac <- as.factor(paste0('G:', data_df$age_group - 1, '_', data_df$age_range))
pmet <- NEON1::site_meta_PUUM_DP110058001[, c('plotID', 'decimalLatitude', 'decimalLongitude', 'elevation')]
pmet <- pmet[!duplicated(pmet), ]
data_df <- dplyr::left_join(data_df, pmet, by = 'plotID')
data_df$cover_logit <- NEON1::logit_Warton(data_df$cover)
data_df$nativeStatusCode[which(data_df$scientificName == 'Rhynchospora chinensis Nees & Meyen')] <- 'N'
#hemi <- NEON1::hemisphere_processed_or_whatever
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
#
# Models tested:
#   beta_family()
#     Sometimes a beta family fixes the grouped residuals, but it generally makes the QQ residuals bow-shaped.
#     A logit_Warton() transform of the response with a gaussian() family improves QQ residuals much better.
#   gaussian()
#     A gaussian link without a response transform has poor convergence and very poor fit. But, a logit_Warton()
#     transform seems to be a major improvement.
#   ziformula = ~taxonID
#     Zero-inflation terms don't seem to help the residuals very much. Perhaps, in this dataset, there isn't
#     much zero-inflation happening at the taxa level, but rather more at the plot level.





# Normal-family, logit-xform, no-ZI
# QQ much better for normality, but still within-group deviance
fit_1 <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = gaussian(),
                 control = c0, start = list())

# Normal-family, logit-xform, ZI
# doesnt appear to change much vs fit_4
fit_2 <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = gaussian(), ziformula = ~taxonID,
                 control = c0, start = list())

# Normal-family, logit-xform, non-ZI, nativeStatusCode as RE
# By-group fit improved, but QQ residuals look worse
fit_3 <- glmmTMB(cover_logit ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = gaussian(),
                 control = c0, start = list())

# Normal-family, logit-xform, ZI, nativeStatusCode as RE
# By-group fit improved, but QQ residuals look worse
fit_4 <- glmmTMB(cover_logit ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = gaussian(), ziformula = ~taxonID,
                 control = c0, start = list())

# Beta-family, ZI, nativeStatusCode as RE
# Same problems as with beta family
fit_5 <- glmmTMB(cover ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = beta_family(), ziformula = ~taxonID,
                 control = c0, start = list())

# Check residuals
cr_grp <- list(data_df$age_fac, data_df$nativeStatusCode, data_df$year, data_df$plotID)
NEON1::check_residuals(fit_1, resid_dir, cr_grp)
NEON1::check_residuals(fit_2, resid_dir, cr_grp)
NEON1::check_residuals(fit_3, resid_dir, cr_grp)
NEON1::check_residuals(fit_4, resid_dir, cr_grp)
NEON1::check_residuals(fit_5, resid_dir, cr_grp)
# 10, 11, 18, 19 are consistent problems


ll <- fit$obj$env$report()$fact_load[[2]] |>
  as.data.frame() |>
  cbind(unique(data_df$taxonID)) |>
  cbind(unique(data_df$taxonID)[order(unique(data_df$taxonID))]) |>
  setNames(c("L1", "L2", "taxonID", "taxonID_ordered"))

results_LV_df <- ll
usethis::use_data(results_LV_df, overwrite = T)

