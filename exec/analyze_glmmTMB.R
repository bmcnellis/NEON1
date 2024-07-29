# BEM 12 May 2023
# updated Feb 2024
# updated July 2024

set.seed(1)

# Caveats:
#     soil age is an ordered factor

### Libraries
library(NEON1)
library(dplyr)
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

### Data import

df0 <- NEON1::div_one |>
  left_join(NEON1::flow_meta, by = 'plotID') |>
  left_join(NEON1::dhp, by = 'plotID') |>
  left_join(NEON1::met, by = 'plotID') |>
  left_join(NEON1::spp, by = 'scientificName') |>
  select(c(plotID, subplotID, endDate, binomialName, age_group, age_range, percentCover, pai, fcover, fipar)) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = percentCover, values_fn = mean) |>
  tidyr::pivot_longer(cols = -c(1:8), names_to = 'binomialName', values_to = 'percentCover') |>
  mutate(
    cov_log = !is.na(percentCover),
    percentCover = ifelse(is.na(percentCover), 0, percentCover),
    endDate = as.Date(endDate),
    date = as.character(endDate),
    plotDate = paste(plotID, date, sep = '_'),
    cover = percentCover / 100,
    cover_logit = NEON1::logit_Warton(cover),
    cover_beta = ifelse(cover > 0.99, 0.99, cover),
    cover_beta = ifelse(cover_beta < 0.001, 0.001, cover_beta),
    age_fac = as.factor(paste0('G:', age_group - 1, '_', age_range))
  ) |>
  select(-c(age_group, age_range, percentCover)) |>
  distinct()
# get native status, this is breaking the pivots for some reason
inv_nat <- NEON1::spp |>
  select(c(binomialName, nativeStatusCode)) |>
  filter(!is.na(nativeStatusCode)) |>
  mutate(nativeBinary = ifelse(nativeStatusCode == 'I', 'I', 'N')) |>
  select(-nativeStatusCode) |>
  mutate(nat_fac = as.factor(nativeBinary)) |>
  distinct()
df0 <- df0 |>
  left_join(inv_nat, by = 'binomialName')

### Modelling

c0 <- glmmTMB::glmmTMBControl(
  optCtrl = list(iter.max = 1e3, eval.max = 1e3),
  start_method = list(method = NULL, jitter.sd = 0)
)

df1 <- df0

# this is glmmTMB's version of a GLVM using a reduced rank covariance structure
# from: https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
#
# Null model: cover is only explained by space, time, or species covariances
# the best fit is gained by nesting date within plotID and an ordbeta() family
# can't figure out how to nest plotID within date for a RR covariance structure
fit_0_1 <- glmmTMB(cover ~ (1|date/plotID),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_1, resid_dir, list(df1$date, df1$plotID))
fit_0_2 <- glmmTMB(cover ~ 1,
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_2, resid_dir, list(df1$date, df1$plotID))
fit_0_3 <- glmmTMB(cover ~ 1 + rr(binomialName + 0|plotID, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_3, resid_dir, list(df1$date, df1$plotID))
fit_0_4 <- glmmTMB(cover ~ (1|date) + rr(binomialName + 0|plotID, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_4, resid_dir, list(df1$date, df1$plotID))
fit_0_5 <- glmmTMB(cover ~ (1|date/plotID) + rr(binomialName + 0|plotID, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_5, resid_dir, list(df1$date, df1$plotID))
fit_0_6 <- glmmTMB(cover ~ (1|plotDate),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_6, resid_dir, list(df1$date, df1$plotID))
fit_0_7 <- glmmTMB(cover ~ 1 + rr(binomialName + 0|plotDate, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_7, resid_dir, list(df1$date, df1$plotID))
fit_0_8 <- glmmTMB(cover ~ (1|date/plotID) + rr(binomialName + 0|plotID, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_8, resid_dir, list(df1$date, df1$plotID))
fit_0_9 <- glmmTMB(cover ~ (1|date/plotID) + rr(binomialName + 0|date/plotID, d = 2),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_0_9, resid_dir, list(df1$date, df1$plotID))

# Model set 1: this explores the data withiout the RR covariance structure and instead uses
# the best-fitting null models
fit_1_1 <- glmmTMB(cover ~ (1|date/plotID),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_1_1, resid_dir, list(df1$date, df1$plotID))
fit_1_2 <- glmmTMB(cover ~ pai*binomialName + (1|date/plotID),
                   data = df1, family = ordbeta(),
                   control = c0, start = list())
NEON1::check_residuals(fit_1_2, resid_dir, list(df1$date, df1$plotID, df1$binomialName))


