# BEM 29 July 2024

set.seed(1)

# TODO:

# Caveats/Considerations:
#     * taxon not identified to species are excluded
#     * two enclosures (Kata-Stein & Unencumbered Kulani Pasture) had NA for completion date, this
#       was set to 2014 which is the age of the youngest bordering enclosure
#     * lava flow age is the median of the age range published on USGS maps
#     * pai, elev, age_median are centered and scaled
#     * time_since_fence NOT centered or scaled
#     * one enclosure (Army Road) was listed as NOT ungulate-free, time_since_fence set to 0

# "All variance partitioning values for each species were multiplied by the explanatory R2 value of that species to show amount
# of total variation in the response variable explained by each covariate."

### Libraries
library(NEON1)
library(dplyr)
library(Hmsc)

#dir0 <- '/media/bem/data/NEON'
dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'

### Directories
data_dir <- dir0
res_dir <- file.path(dir0, 'results')
mod_dir <- file.path(dir0, 'results/model_results')
stopifnot(all(dir.exists(data_dir), dir.exists(res_dir), dir.exists(mod_dir)))

### Data import

inv_nat <- NEON1::spp |>
  select(c(binomialName, nativeStatusCode)) |>
  filter(!is.na(nativeStatusCode)) |>
  distinct()
fence <- NEON1::fence |>
  mutate(enclosure_completed = ifelse(enclosure_name == 'Kata-Stein', 2014, enclosure_completed)) |>
  mutate(time_since_fence = 2024 - enclosure_completed) |>
  mutate(time_since_fence = ifelse(enclosure_name == 'Army Road', 0, time_since_fence)) |>
  mutate(time_since_fence = ifelse(is.na(time_since_fence), 0, time_since_fence)) |>
  select(c(plotID, time_since_fence)) |>
  distinct()
df0 <- NEON1::div_one |>
  # convert scientificName to binomialName
  left_join(NEON1::spp, by = 'scientificName') |>
  # pivot the dataframe to fill in 0's, then aggregate across subplots
  select(c(plotID, subplotID, endDate, binomialName, percentCover)) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = percentCover, values_fn = mean) |>
  tidyr::pivot_longer(cols = -c(1:3), names_to = 'binomialName', values_to = 'percentCover') |>
  mutate(percentCover = ifelse(is.na(percentCover), 0, percentCover)) |>
  select(-subplotID) |>
  group_by(plotID, endDate, binomialName) |>
  summarise(percentCover = mean(percentCover, na.rm = T), .groups = 'drop') |>
  # join the other relevant data
  left_join(inv_nat, by = 'binomialName') |>
  left_join(fence, by = 'plotID') |>
  left_join(NEON1::met[, c('plotID', 'elevation', 'decimalLatitude', 'decimalLongitude')], by = 'plotID') |>
  left_join(NEON1::flow_meta[, c('plotID','age_median')], by = 'plotID') |>
  left_join(NEON1::dhp[, c('plotID', 'pai')], by = 'plotID') |>
  left_join(NEON1::lus[, c('plotID', 'cover_type', 'log', 'cow')], by = 'plotID') |>
  # add data from mutates
  mutate(plotDate = paste(plotID, as.character(endDate), sep = '_')) |>
  mutate(plotDate = gsub('PUUM_', '', plotDate)) |>
  mutate(date_day_center = as.integer(round(scale(as.numeric(endDate) / 86400, scale = F)))) |>
  # drop species which are not identified to species
  filter(!grepl('_sp|_spp', binomialName)) |>
  # assume all NI? are N
  mutate(nativeStatusCode = ifelse(nativeStatusCode == 'NI?', 'N', 'I')) |>
  mutate(nativeStatusCode = ifelse(nativeStatusCode == 'I', 'z_I', nativeStatusCode)) |>
  # center/scale and fix for Hmsc inputs
  mutate(
    pai = as.numeric(scale(pai)),
    elevation = as.numeric(scale(elevation)),
    decimalLatitude = as.numeric(scale(decimalLatitude, scale = F)),
    decimalLongitude = as.numeric(scale(decimalLatitude, scale = F)),
    age_median = as.numeric(scale(age_median)),
    cover_type = ifelse(cover_type == 'ohia_tall', 'AAA_ohia_tall', cover_type),
    cover_type = as.factor(cover_type)
  ) |>
  distinct()
rm(inv_nat, fence)

# create species matrix for presence-absence, plot-date as row with species as col
mat_p <- df0 |>
  mutate(fractionCover = percentCover / 100) |>
  select(c(plotDate, binomialName, fractionCover)) |>
  group_by(plotDate, binomialName) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = fractionCover) |>
  arrange(plotDate) |>
  tibble::column_to_rownames('plotDate') |>
  as.matrix() |>
  apply(c(1, 2), \(xx) ifelse(xx > 0.0001, 1, 0))
# create environmental matrix
env1 <- df0 |>
  select(c(plotDate, age_median, pai, time_since_fence, elevation, cover_type, log, cow)) |>
  arrange(plotDate) |>
  distinct() |>
  tibble::column_to_rownames('plotDate') |>
  as.data.frame()
# create trait df
tr1 <- df0 |>
  select(binomialName, nativeStatusCode) |>
  distinct() |>
  filter(binomialName %in% colnames(mat_p)) |>
  arrange(binomialName) |>
  tibble::column_to_rownames('binomialName')
# create study-design matrix
stu0 <- df0 |>
  select(c(plotDate, plotID, date_day_center)) |>
  distinct() |>
  arrange(plotDate) |>
  mutate(plotDate = as.factor(plotDate)) |>
  select(plotDate) |>
  mutate(p0 = plotDate) |>
  tibble::column_to_rownames('p0') |>
  as.data.frame()

# check matrix/df alignment and modify for modelling
stopifnot(
  identical(row.names(mat_p), row.names(env1)),
  identical(row.names(mat_p), row.names(stu0)),
  identical(colnames(mat_p), row.names(tr1))
)

# model terms
xf1 <- as.formula(~ age_median + pai + time_since_fence + elevation + log + cow + cover_type)
rL1 <- Hmsc::HmscRandomLevel(units = stu0$plotDate)
# "Ultimately, the random effect structure one uses in an analysis encodes the assumptions that one makes about
# how sampling units (subjects and items) vary, and the structure of dependency that this variation creates in one’s data."
# - Barr et al. 2014

# fit
# presence/absence probit model
if (!file.exists(file.path(mod_dir, 'm_p_mod.rda'))) {
  m_p <- Hmsc::Hmsc(
    Y = mat_p, XData = env1, studyDesign = stu0, XFormula = xf1, TrData = tr1,
    ranLevels = list('plotDate' = rL1), distr = 'probit', TrFormula = ~nativeStatusCode,
    XScale = F
  )
  #m_p <- Hmsc::sampleMcmc(m_p, thin = 10, samples = 12000, transient = 2000, nChains = 6, nParallel = 6)
  m_p <- Hmsc::sampleMcmc(m_p, thin = 10, samples = 4000, transient = 1000, nChains = 4, nParallel = 2)
  mc_p <- Hmsc::convertToCodaObject(m_p)
  ma_p <- Hmsc::computeAssociations(m_p, thin = 10)
  # requires normality of explanatory variables i think
  mp_p <- Hmsc::computePredictedValues(m_p)
  mp_pp <- predict(m_p)
  save(list = c('m_p', 'mc_p', 'ma_p', 'mp_p', 'mp_pp'), file = file.path(mod_dir, 'm_p_mod.rda'))
}

# presence/absence probit model
if (!file.exists(file.path(mod_dir, 'm_p_diag.rda'))) {

  mf_p <- Hmsc::evaluateModelFit(hM = m_p, predY = mp_p)
  m_vp <- Hmsc::computeVariancePartitioning(m_p)
  m_pv <- Hmsc::computePredictedValues(m_p)
  m_ca <- Hmsc::computeAssociations(m_p)

  es_p_beta <- coda::effectiveSize(mc_p$Beta)
  es_p_gamm <- coda::effectiveSize(mc_p$Gamma)
  es_p_omeg <- coda::effectiveSize(mc_p$Omega[[1]])

  gd_p_beta <- coda::gelman.diag(mc_p$Beta, multivariate = F)$psrf
  gd_p_gamm <- coda::gelman.diag(mc_p$Gamma, multivariate = F)$psrf
  gd_p_omeg <- coda::gelman.diag(mc_p$Omega[[1]], multivariate = F)$psrf

  mpe_p_beta <- Hmsc::getPostEstimate(m_p, parName = 'Beta')
  mpe_p_gamm <- Hmsc::getPostEstimate(m_p, parName = 'Gamma')
  mpe_p_omeg <- Hmsc::getPostEstimate(m_p, parName = 'Omega')

  #ml_pr <- getPostEstimate(m_p, parName = 'Lambda')
  ml_pb <- Hmsc::getPostEstimate(m_p, parName = 'Beta')
  ml_po <- Hmsc::getPostEstimate(m_p, parName = 'Omega')
  ml_pg <- Hmsc::getPostEstimate(m_p, parName = 'Gamma')


  save(list = c(
    'mf_p', 'm_ca', 'm_vp',# 'ml_pr'
    'es_p_beta', 'es_p_omeg', 'es_p_gamm',
    'gd_p_beta', 'gd_p_omeg', 'gd_p_gamm',
    'mpe_p_beta', 'mpe_p_omeg', 'mpe_p_gamm',
    'ml_pb', 'ml_po', 'ml_pg'
  ), file = file.path(mod_dir, 'm_p_diag.rda'))

}
