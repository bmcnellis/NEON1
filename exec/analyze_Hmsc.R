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

### Data import

inv_nat <- NEON1::spp |>
  select(c(binomialName, nativeStatusCode)) |>
  filter(!is.na(nativeStatusCode)) |>
  mutate(nativeBinary = ifelse(nativeStatusCode == 'I', 'I', 'N')) |>
  select(-nativeStatusCode) |>
  mutate(nat_fac = as.factor(nativeBinary)) |>
  distinct()
df0 <- NEON1::div_one |>
  left_join(NEON1::spp, by = 'scientificName') |>
  select(c(plotID, subplotID, endDate, binomialName, percentCover)) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = percentCover, values_fn = mean) |>
  tidyr::pivot_longer(cols = -c(1:3), names_to = 'binomialName', values_to = 'percentCover') |>
  mutate(percentCover = ifelse(is.na(percentCover), 0, percentCover)) |>
  select(-subplotID) |>
  group_by(plotID, endDate, binomialName) |>
  summarise(percentCover = mean(percentCover, na.rm = T), .groups = 'drop') |>
  left_join(inv_nat, by = 'binomialName') |>
  left_join(NEON1::flow_meta, by = 'plotID') |>
  left_join(NEON1::dhp, by = 'plotID') |>
  left_join(NEON1::met, by = 'plotID') |>
  mutate(plotDate = paste0(gsub('PUUM_', '', plotID), as.character(endDate), sep = '_')) |>
  mutate(ageFactor = as.factor(paste0('G:', age_group - 1, '_', age_range))) |>
  rename('nativeFactor' = nat_fac) |>
  select(plotDate, plotID, endDate, binomialName, percentCover, ageFactor, nativeFactor, pai, fcover, fipar) |>
  mutate(logicalCover = ifelse(percentCover < 0.0001, FALSE, TRUE)) |>
  mutate(fractionCover = percentCover / 100) |>
  mutate(betaCover = ifelse(fractionCover > 0.99, 0.99, fractionCover)) |>
  mutate(betaCover = ifelse(betaCover < 0.001, 0.001, betaCover)) |>
  mutate(logitCover = NEON1::logit_Warton(fractionCover)) |>
  distinct()
rm(inv_nat)

# create species matrix, plot-date as row with species as col
mat0 <- df0 |>
  select(c(plotDate, binomialName, fractionCover)) |>
  group_by(plotDate, binomialName) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = fractionCover) |>
  arrange(plotDate) |>
  tibble::column_to_rownames('plotDate') |>
  as.matrix()
# create environmental matrix
env1 <- df0 |>
  select(c(plotDate, ageFactor, pai)) |>
  arrange(plotDate) |>
  distinct()
# create study-design matrix
stu1 <- df0 |>
  select(c(plotDate, plotID, endDate)) |>
  arrange(plotDate) |>
  select(plotDate, plotID) |>
  distinct()

stopifnot(nrow(env1) == nrow(mat0), nrow(env1) == nrow(stu1))
stopifnot(identical(row.names(mat0), env1$plotDate))
stopifnot(identical(row.names(mat0), stu1$plotDate))
row.names(mat0) <- NULL
env1 <- env1[, -1]
stu1 <- stu1[, -1]

# modelling

xf1 <- as.formula(~ ageFactor)

mm <- Hmsc::Hmsc(Y = mat0, XData = env1, XFormula = xf1)
mm <- Hmsc::sampleMcmc(m,thin = 10, samples = 1000, transient = 500, nChains = 2, nParallel = 2)

mp <- Hmsc::convertToCodaObject(mm)
summary(mp$Beta)
