# BEM 29 July 2024

set.seed(1)

# Caveats:
#     soil age is an ordered factor

### Libraries
library(NEON1)
library(dplyr)
library(Hmsc)

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


