# BEM 12 May 2023
# updated 6 Feb 2023

set.seed(1)

# Caveats:
#     soil age is an ordered factor
#     'Rhynchospora chinensis' assumed native, b/c native ssp, unsure on other ssp presence

### Libraries
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
library(NEON1)
library(glmmTMB)
library(DHARMa)

### Directories
#data_dir <- '/media/bem/data/NEON'
data_dir <- 'C:/Users/BrandonMcNellis/Documents/NEON_data'
resid_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/reports'

### Constants

omit_codes_1 <- c('2PLANT')
omit_codes_2 <- c(
  "2PLANT", "ASPLE", "CAREX", "FRAGA", "MELIC3", "PEPER", "RHYNC3", "DRYOP", "ERAGR",
  "CYANE", "CYRTA", "SANTA", "SELAG", "STENO", "CLERM", "PRITC", "FRAGASPP", "MYRSI",
  "DIPLA2", "CAREXSPP", "MELIC3SPP", "MICRO3", "ASPLESPP", "DRYOPSPP", "RUBUS"
)

c0 <- glmmTMB::glmmTMBControl(
  optCtrl = list(iter.max = 1e3, eval.max = 1e3),
  start_method = list(method = NULL, jitter.sd = 0)
)

### Import & modification
#PUUM_div <- neonPlantEcology::npe_download(sites = c("PUUM"))
#usethis::use_data(PUUM_div, overwrite = T)  # "2024-05-13 12:32:03 HST"

data_div <- NEON1::PUUM_div
data_df <- NEON1::div_to_long(data_div, add_zeros = T)
data_df <- dplyr::left_join(data_df, NEON1::flow_meta, by = 'plotID')
# add hemisphere data
#hemi <- NEON1::hemisphere_processed_or_whatever
# add traits
#traits <- neonUtilities::stackByTable(filepath = file.path(data_dir, 'NEON_traits-foliar.zip'), savepath = 'envt')
colnames(data_df)[which(colnames(data_df) == 'percentCover')] <- 'cover'
data_df$cover <- data_df$cover / 100
data_df$cover_logit <- NEON1::logit_Warton(data_df$cover)
data_df$age_fac <- as.factor(paste0('G:', data_df$age_group - 1, '_', data_df$age_range))
data_df$nativeStatusCode[which(data_df$scientificName == 'Rhynchospora chinensis Nees & Meyen')] <- 'N'
data_df$nativeStatusCode <- ifelse(data_df$nativeStatusCode == 'N', '0_N', data_df$nativeStatusCode)
data_df$nativeStatusCode <- ifelse(data_df$nativeStatusCode == 'I', '1_I', data_df$nativeStatusCode)
data_df$nativeStatusCode <- ifelse(data_df$nativeStatusCode == 'UNK', '2_U', data_df$nativeStatusCode)
data_df$nativeStatusCode <- as.factor(data_df$nativeStatusCode)
data_df$endDate <- as.character(data_df$endDate)

# set up multiple input dataframes to check for omission changes
data_df_0 <- data_df
data_df_1 <- data_df[-which(data_df$taxonID %in% omit_codes_1), ]
data_df_2 <- data_df[-which(data_df$taxonID %in% omit_codes_2), ]
data_df_2$nativeStatusCode <- as.character(data_df_2$nativeStatusCode)
data_df_2$nativeStatusCode <- as.factor(data_df_2$nativeStatusCode)

# set up multiple input dataframes to use for hurdles
data_df_0_Z <- within(data_df_0, cov_0 <- cover < 0.0001)
data_df_0_NZ <- data_df_0[which(data_df_0$cover > 0.0001), ]

### Modelling

# this is glmmTMB's version of a GLVM using a reduced rank covariance structure
# from: https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
#
# Models tested:
#   beta_family()
#     Sometimes a beta family fixes the grouped residuals, but it generally makes the QQ residuals bow-shaped.
#     A logit_Warton() transform of the response with a gaussian() family improves QQ residuals much better.
#     Dropping the weird taxons didnt help.
#   gaussian()
#     A gaussian link without a response transform has poor convergence and very poor fit. But, a logit_Warton()
#     transform seems to be a major improvement. Dropping the weird taxons didnt help.
#   ziformula = ~taxonID
#     Zero-inflation terms don't seem to help the residuals very much. Perhaps, in this dataset, there isn't
#     much zero-inflation happening at the taxa level, but rather more at the plot level. Dropping the weird
#     taxa didnt help.
#     Brooks et al. (2017) says that the zi-inflation terms in glmmTMB are *structural* zeroes, whereas zeroes
#     in this data is probably *sampling* zeroes, i.e. the sampling units are too small relative to the overall
#     spatial percent cover of the species across the plots.
#   (1 | nativeStatusCode)
#     Native status as a random effect seemed to worsen the QQ normals, simliar to changing the family to
#     beta. It makes sense just fine as a fixed effect, with natives = 0 and value changes for invasive and unk, respectively.
#   nativeStatusCode
#     Native residuals are good, but the invasives are super heteroscedastic, with a tail towards low cover and outliers
#     at high cover. It looks like most invasives behave as natives, but some are behaving as invasives. Dropping the weird
#     taxa didnt help, and the invasives still appear to be a problem.
#   dispformula
#     The overall dispersion was made worse using dispformula = ~nativeStatusCode or ~nativeStatusCode + year. This
#     was the same for gaussian or beta families.

# Normal-family, logit-xform, no-ZI
# QQ much better for normality, but still within-group deviance
# only real problem groups are nativeStatusCode, because invasives are weird
# Specifically, high predicted cover values have lower scaled residuals than they should
fit_0 <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | endDate) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df_0_NZ, family = gaussian(),
                 control = c0, start = list())

fit_1 <- glmmTMB(cov_0 ~ age_fac + nativeStatusCode + (1 | endDate) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df_0_Z, family = binomial(),
                 control = c0, start = list())


# hurdle model fit to presence/absence from other plot sizes???

fit_XX <- glmmTMB(cover ~ age_fac + nativeStatusCode + (1 | endDate) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df_0, family = Gamma(),
                 control = c0, start = list())

fit_XX <- glmmTMB(cover_logit ~ age_fac + nativeStatusCode + (1 | endDate) + rr(taxonID + 0|plotID, d = 2),
                 data = data_df_0, family = gaussian(), ziformula = ~ 1
                 control = c0, start = list())
# Gamma can have dispersion parameter, gaussian can not, explore this

# Check residuals
cr_grp_0 <- list(data_df_0$age_fac, data_df_0$nativeStatusCode, data_df_0$endDate, data_df_0$plotID)
cr_grp_1 <- list(data_df_1$age_fac, data_df_1$nativeStatusCode, data_df_1$endDate, data_df_1$plotID)
cr_grp_2 <- list(data_df_2$age_fac, data_df_2$nativeStatusCode, data_df_2$endDate, data_df_2$plotID)

NEON1::check_residuals(fit_0, resid_dir, cr_grp_0)
NEON1::check_residuals(fit_1, resid_dir, cr_grp_0)
NEON1::check_residuals(fit_2, resid_dir, cr_grp_0)
# 10, 11, 18, 19 are consistent problems


ll <- fit_0$obj$env$report()$fact_load[[2]] |>
  as.data.frame() |>
  cbind(unique(data_df_0$taxonID)) |>
  cbind(unique(data_df_0$taxonID)[order(unique(data_df_0$taxonID))]) |>
  setNames(c("L1", "L2", "taxonID", "taxonID_ordered"))

#results_LV_df <- ll
#usethis::use_data(results_LV_df, overwrite = T)

