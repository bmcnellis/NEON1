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
# Beta-family, ZI
fit_0 <- glmmTMB(cover ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = beta_family(), ziformula = ~taxonID,
                 control = c0, start = list())
#
# Beta-family, no-ZI
# omitting zero-inflation term causes strong grouping within the residuals
fit_1 <- glmmTMB(cover ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = beta_family(),
                 control = c0, start = list())
#
# Normal-family, ZI
# poor convergence, bad fit
#fit_2 <- glmmTMB(cover ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = gaussian(), ziformula = ~taxonID,
#                 control = c0, start = list())
#
# Normal-family, no-ZI
# poor convergence, bad fit
#fit_3 <- glmmTMB(cover ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = gaussian(),
#                 control = c0, start = list())
#
# Normal-family, logit-xform, no-ZI
# QQ much better for normality, but still within-group deviance
fit_4 <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df, family = gaussian(),
                 control = c0, start = list())
#
# Normal-family, logit-xform, ZI
# doesnt appear to change much vs fit_4
#fit_5 <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = gaussian(), ziformula = ~taxonID,
#                 control = c0, start = list())
#
# Normal-family, logit-xform, non-ZI, nativeStatusCode as RE
# By-group fit improved, but QQ residuals look worse
#fit_6 <- glmmTMB(cover_logit ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = gaussian(),
#                 control = c0, start = list())
#
# Normal-family, logit-xform, ZI, nativeStatusCode as RE
# By-group fit improved, but QQ residuals look worse
#fit_7 <- glmmTMB(cover_logit ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = gaussian(), ziformula = ~taxonID,
#                 control = c0, start = list())
#
# Beta-family, ZI, nativeStatusCode as RE
# Same problems as with beta family
#fit_8 <- glmmTMB(cover ~ age_fac + (1 | nativeStatusCode) + (1 | year) + rr(taxonID + 0|plotID, d = 2),
#                 data = data_df, family = beta_family(), ziformula = ~taxonID,
#                 control = c0, start = list())

r_0 <- resid(fit_0)
s_0 <- DHARMa::simulateResiduals(fit_0)
plot(s_0, rank = T)
plot(s_0, asFactor = T)

r_1 <- resid(fit_1)
s_1 <- DHARMa::simulateResiduals(fit_1)
plot(s_1, rank = T)
plot(s_1, asFactor = T)

#r_2 <- resid(fit_2)
#s_2 <- DHARMa::simulateResiduals(fit_2)
#plot(s_2, rank = T)
#plot(s_2, asFactor = T)

#r_3 <- resid(fit_3)
#s_3 <- DHARMa::simulateResiduals(fit_3)
#plot(s_3, rank = T)
#plot(s_3, asFactor = T)

r_4 <- resid(fit_4)
s_4 <- DHARMa::simulateResiduals(fit_4)
plot(s_4, rank = T)
plot(s_4, asFactor = T)
# plotResiduals should look like a no-sig-difference boxplot, with mean at 0.5 and box from 0.25-0.75
plotResiduals(s_4, form = data_df$age_fac)
# looks good, first 3 age groups deviate within-group
plotResiduals(s_4, form = data_df$nativeStatusCode)
# looks terrible, invasives are a big problem
plotResiduals(s_4, form = data_df$year)
# could be worse, all deviate within-group
plotResiduals(s_4, form = data_df$plotID)
# 2, 10, 11, 18, 19 are bad
# 10, 11, 18, 19 are consistent problems

#r_5 <- resid(fit_5)
#s_5 <- DHARMa::simulateResiduals(fit_5)
#plot(s_5, rank = T)
#plot(s_5, asFactor = T)

#r_6 <- resid(fit_6)
#s_6 <- DHARMa::simulateResiduals(fit_6)
#plot(s_6, rank = T)
#plot(s_6, asFactor = T)
#plotResiduals(s_6, form = data_df$age_fac)
#plotResiduals(s_6, form = data_df$nativeStatusCode)
#plotResiduals(s_6, form = data_df$year)
#plotResiduals(s_6, form = data_df$plotID)

#r_7 <- resid(fit_7)
#s_7 <- DHARMa::simulateResiduals(fit_7)
#plot(s_7, rank = T)
#plot(s_7, asFactor = T)
#plotResiduals(s_7, form = data_df$age_fac)
#plotResiduals(s_7, form = data_df$nativeStatusCode)
#plotResiduals(s_7, form = data_df$year)
#plotResiduals(s_7, form = data_df$plotID)

#r_8 <- resid(fit_8)
#s_8 <- DHARMa::simulateResiduals(fit_8)
#plot(s_8, rank = T)
#plot(s_8, asFactor = T)
#plotResiduals(s_8, form = data_df$age_fac)
#plotResiduals(s_8, form = data_df$nativeStatusCode)
#plotResiduals(s_8, form = data_df$year)
#plotResiduals(s_8, form = data_df$plotID)

ll <- fit$obj$env$report()$fact_load[[2]] |>
  as.data.frame() |>
  cbind(unique(data_df$taxonID)) |>
  cbind(unique(data_df$taxonID)[order(unique(data_df$taxonID))]) |>
  setNames(c("L1", "L2", "taxonID", "taxonID_ordered"))

results_LV_df <- ll
usethis::use_data(results_LV_df, overwrite = T)

