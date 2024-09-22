# BEM
# created  29 July 2024
# last run 21 Sep 2024

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
library(coda)
library(MCMCvis)

# model parameters
#s0 <- c(6000, 2000, 4, 4)# iterations, burn, chains, parallel
# need at least 4k burn-in and more than 4 chains for problematic species (e.g. CHITRI)
s0 <- c(8000, 4000, 6, 6)
# add more when final model is figured out

dir0 <- '/media/bem/data/NEON'
#dir0 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1'

### Directories
data_dir <- dir0
res_dir <- file.path(dir0, 'results')
mod_dir <- file.path(dir0, 'results/model_results')
stopifnot(all(dir.exists(data_dir), dir.exists(res_dir), dir.exists(mod_dir)))

# DONT REMOVE THE COMMENTED ROUNDS
# round 1: only CHETRI removed, oc >= 8, 7
#bad_spp_0 <- c('Cheirodendron_trigynum')
#bad_spp_1 <- c('Cheirodendron_trigynum')
# round 2: nothing removed, oc >= 15, 10
#bad_spp_0 <- character()
#bad_spp_1 <- character()
# beta 0 visu: Anthoxanthum, Axonopus, Cheirodendron, C. glaucum, Eragrostis, Hedyotis, Hymenophyllum_recurvum, Hypochaeris_radicata
# beta 0 rhat: Hypochaeris_radicata, Axonopus_fissifolius, and Eragrostis_cumingiiare are the worst offenders. They are weedy non-natives.
#              can take out Hypochaeris_radicata first because it probably isn't a species of serious concern
# gamm 0 visu: OK
# gamm 0 rhat: because all the beta0 problems are invasives, this will likely be improved
# beta 1 visu: Adenophorus, Athyrium, C. menzeisii, Cyrtandra spp., Elaphoglossum alatum, Freycinetia, Hedyotis, Hymenophyllum spp., Stenogyne
# beta 1 rhat: many species, can remove from df 0 first and see how that changes things
# gamm 0 visu: OK
# gamm 0 rhat: two problems with intercept, otherwise OK
# round 3: H radicata removed, oc >= 15, 10
#bad_spp_0 <- c('Hypochaeris_radicata')
#bad_spp_1 <- c('Hypochaeris_radicata')
# removing H radicata made beta 0 a lot better and fixed the gamm issues
# beta 0 rhat: Hymenophyllum_recurvum has highest Rhat
# beta 1 rhat: Hymenophyllum_recurvum also present in bad beta 1's
# round 4: H radicata, H recurvum removed, oc >= 15, 10
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum')
#bad_spp_1 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum')
# beta 0 rhat: everything mostly got worse? try dropping Axonopus
# beta 1 rhat: lots of stuff still bad, try dropping Sticherus_owhyensis
# round 5: below removed, oc >=15, 10
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius')
#bad_spp_1 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Sticherus_owhyensis')
# only 8 bad terms for beta 0, with 7 species, best yet. can try removing Eragrostis_cumingii
# round 6: removing Axonopus helped, removing Sticherus broke
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius', 'Eragrostis_cumingii')
#bad_spp_1 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius')
# round 7: misspelled Axonopus, starting over with beta 0
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius', 'Eragrostis_cumingii')
#bad_spp_1 <- c('')
# round 8: removing eragrostis made it worse, trying again with beta 1
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius')
#bad_spp_1 <- c('Hymenophyllum_lanceolatum')
# round 9: maybe this is as good as 0 is going to get, try 1 with Coprosma because of Omega1 rhats
#bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius')
#bad_spp_1 <- c('Coprosma_granadensis')
# round 10: 0 left unchanged, 1 oc = 7, no drops, s0 raised to [8000, 4000, 6, 6]
bad_spp_0 <- c('Hypochaeris_radicata', 'Hymenophyllum_recurvum', 'Axonopus_fissifolius')
bad_spp_1 <- c('')

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
# combine div_one and div_ten?
df0 <- NEON1::div_one |>
  # 125 species to start
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
  mutate(year = as.integer(as.POSIXlt(endDate)$year + 1900)) |>
  # drop species which are not identified to species
  # 122 species
  filter(!grepl('_sp|_spp', binomialName))
  # 93 species
# drop species which have less than 8 occurences in the entire dataset
# remaining species represent > 5% of the total possible occurences (142)
drop_spp <- df0 |>
  group_by(binomialName) |>
  summarise(oc = sum(ifelse(percentCover > 0.001, 1, 0))) |>
  # 8 / 142 means remaining represent > 5% of total possible occurences
  #filter(oc >= 8) |>
  # 15 / 142 means remaining represent > 10% of total possible occurences
  filter(oc >= 15) |>
  pull(binomialName)
df0 <- df0 |>
  # 93 species
  filter(binomialName %in% drop_spp)
# 45 species if oc >= 8
# 34 species if oc >= 15
if (length(bad_spp_0) > 0) {
  df0 <- df0 |>
    filter(!(binomialName %in% bad_spp_0))
}
df0 <- df0 |>
  # assume all NI? are N
  mutate(nativeStatusCode = ifelse(nativeStatusCode %in% c('N', 'NI?'), 'N', 'I')) |>
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
df1 <- NEON1::div_one |>
  left_join(NEON1::spp, by = 'scientificName') |>
  select(plotID, subplotID, endDate, taxonID) |>
  filter(endDate %in% NEON1::div_ten$endDate) |>
  rbind(NEON1::div_ten) |>
  distinct() |>
  group_by(plotID, endDate, taxonID) |>
  summarise(count = n(), .groups = 'drop') |>
  left_join(NEON1::spp, by = 'taxonID') |>
  # pivot the dataframe to fill in 0's, then aggregate across subplots
  select(c(plotID, endDate, binomialName, count)) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = count, values_fn = sum) |>
  tidyr::pivot_longer(cols = -c(1:2), names_to = 'binomialName', values_to = 'count') |>
  mutate(count = ifelse(is.na(count), 0, count)) |>
  group_by(plotID, endDate, binomialName) |>
  summarise(count = sum(count), .groups = 'drop') |>
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
  mutate(year = as.integer(as.POSIXlt(endDate)$year + 1900)) |>
  # drop species which are not identified to species
  # 175 species
  filter(!grepl('_sp|_spp', binomialName))
  # 134 species
drop_spp <- df1 |>
  group_by(binomialName) |>
  summarise(oc = sum(ifelse(count > 0, 1, 0))) |>
  # 5 / 93 means remaining represent > 5% of total possible occurences
  filter(oc >= 5) |>
  # 10 / 93 means remaining represent > 10% of total possible occurences
  #filter(oc >= 10) |>
  pull(binomialName)
df1 <- df1 |>
  # 134 species
  filter(binomialName %in% drop_spp)
# 88 species if oc >= 5
# 73 species if oc >= 10
if (length(bad_spp_1) > 0) {
  df1 <- df1 |>
    filter(!(binomialName %in% bad_spp_1))
}
df1 <- df1 |>
  mutate(nativeStatusCode = ifelse(nativeStatusCode %in% c('N', 'NI?'), 'N', 'I')) |>
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

rm(inv_nat, fence, drop_spp)

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
mat_t <- df1 |>
  mutate(count = as.integer(count)) |>
  select(c(plotDate, binomialName, count)) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = count) |>
  arrange(plotDate) |>
  tibble::column_to_rownames('plotDate') |>
  as.matrix() |>
  apply(c(1, 2), \(xx) ifelse(xx > 0, 1, 0))
# create environmental matrix
env0 <- df0 |>
  select(c(plotDate, age_median, pai, time_since_fence, elevation, cover_type, log, cow)) |>
  arrange(plotDate) |>
  distinct() |>
  tibble::column_to_rownames('plotDate') |>
  as.data.frame()
env1 <- df1 |>
  select(c(plotDate, age_median, pai, time_since_fence, elevation, cover_type, log, cow)) |>
  arrange(plotDate) |>
  distinct() |>
  tibble::column_to_rownames('plotDate') |>
  as.data.frame()
# create trait df
tr0 <- df0 |>
  select(binomialName, nativeStatusCode) |>
  distinct() |>
  filter(binomialName %in% colnames(mat_p)) |>
  arrange(binomialName) |>
  tibble::column_to_rownames('binomialName')
tr1 <- df1 |>
  select(binomialName, nativeStatusCode) |>
  distinct() |>
  filter(binomialName %in% colnames(mat_t)) |>
  arrange(binomialName) |>
  tibble::column_to_rownames('binomialName')
# create study-design matrix
stu0 <- df0 |>
  #select(c(plotDate, plotID, year)) |>
  select(plotDate) |>
  distinct() |>
  arrange(plotDate) |>
  mutate(across(everything(), as.factor)) |>
  mutate(p0 = plotDate) |>
  tibble::column_to_rownames('p0') |>
  as.data.frame()
stu1 <- df1 |>
  #select(c(plotDate, plotID, year)) |>
  select(plotDate) |>
  distinct() |>
  arrange(plotDate) |>
  mutate(across(everything(), as.factor)) |>
  mutate(p0 = plotDate) |>
  tibble::column_to_rownames('p0') |>
  as.data.frame()

# check matrix/df alignment and modify for modelling
stopifnot(
  identical(row.names(mat_p), row.names(env0)),
  identical(row.names(mat_p), row.names(stu0)),
  identical(colnames(mat_p), row.names(tr0)),
  identical(row.names(mat_t), row.names(env1)),
  identical(row.names(mat_t), row.names(stu1)),
  identical(colnames(mat_t), row.names(tr1))
)

# model terms
xf1 <- as.formula(~ age_median + pai + time_since_fence + elevation + log + cow + cover_type)
rL0 <- list('plotDate' = Hmsc::HmscRandomLevel(units = stu0$plotDate))
#rL0 <- list("plotDate" = Hmsc::HmscRandomLevel(units = stu0$plotDate), "plotID" = Hmsc::HmscRandomLevel(units = unique(stu0$plotID)), "year" = Hmsc::HmscRandomLevel(units = unique(stu0$year)))
#rL1 <- list("plotDate" = Hmsc::HmscRandomLevel(units = stu1$plotDate), "plotID" = Hmsc::HmscRandomLevel(units = unique(stu1$plotID)), "year" = Hmsc::HmscRandomLevel(units = unique(stu1$year)))
rL1 <- list('plotDate' = Hmsc::HmscRandomLevel(units = stu1$plotDate))

# "Ultimately, the random effect structure one uses in an analysis encodes the assumptions that one makes about
# how sampling units (subjects and items) vary, and the structure of dependency that this variation creates in one’s data."
# - Barr et al. 2014

# example script code:
#study_design = str_split(rownames(X_250), pattern = "_", simplify = T) %>%
#  as_tibble() %>%
#  mutate_each(funs = as.integer) %>%
#  set_names(c("samplearea_id_250", "year")) %>%
#  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_250")]), by = "samplearea_id_250") %>%
#  dplyr::select(samplearea_id, samplearea_id_250, year) %>%
#  mutate_each(funs = as_factor)
#random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
#                     "samplearea_id_250" = HmscRandomLevel(units = unique(study_design$samplearea_id_250)),
#                     "year" = HmscRandomLevel(units = unique(study_design$year)))
#model_occ_250 = Hmsc(Y = as.matrix(Y_250_occ), XData = X_250, XFormula = ~., distr = "probit",
#                     studyDesign = as.data.frame(study_design), ranLevels = random_levels)
#model_abund_250 = Hmsc(Y = as.matrix(Y_250_abund), XData = X_250, XFormula = ~., distr = "lognormal poisson",
#                       studyDesign = as.data.frame(study_design), ranLevels = random_levels)

# fit
# presence/absence probit model
# using only 1-m data
m_p_0 <- Hmsc::Hmsc(
  Y = mat_p, XData = env0, studyDesign = stu0, XFormula = xf1, TrData = tr0, ranLevels = rL0,
  distr = 'probit', TrFormula = ~nativeStatusCode, XScale = F
)
m_p_0 <- Hmsc::sampleMcmc(m_p_0, thin = 10, samples = s0[1], transient = s0[2], nChains = s0[3], nParallel = s0[4])
mc_p_0 <- Hmsc::convertToCodaObject(m_p_0)
ma_p_0 <- Hmsc::computeAssociations(m_p_0, thin = 10)
mp_p_0 <- Hmsc::computePredictedValues(m_p_0)
mp_pp_0 <- predict(m_p_0)

mf_p_0 <- Hmsc::evaluateModelFit(hM = m_p_0, predY = mp_p_0)
m_vp_0 <- Hmsc::computeVariancePartitioning(m_p_0)
m_pv_0 <- Hmsc::computePredictedValues(m_p_0)
m_ca_0 <- Hmsc::computeAssociations(m_p_0)

mc_s_beta_0 <- MCMCvis::MCMCsummary(mc_p_0$Beta)
mc_s_gamm_0 <- MCMCvis::MCMCsummary(mc_p_0$Gamma)
mc_s_omeg_0 <- MCMCvis::MCMCsummary(mc_p_0$Omega[[1]])

post_beta_0 <- NEON1::posterior_from_coda(mc_p_0, 'Beta', 0.9, average = T, drop_ns = T)
post_gamm_0 <- NEON1::posterior_from_coda(mc_p_0, 'Gamma', 0.9, average = T, drop_ns = T)
post_omeg_0 <- NEON1::posterior_from_coda(mc_p_0, 'Omega', 0.9, average = T, drop_ns = T)

## Gelman's PSRF (Potential Scale Reduction Factor)
# Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations.
# Journal of Computational and Graphical Statistics, 7, 434-455.
#gd_p_beta_0 <- coda::gelman.diag(mc_p_0$Beta, multivariate = T)$psrf
#gd_p_gamm_0 <- coda::gelman.diag(mc_p_0$Gamma, multivariate = T)$psrf
# multivariate = T fails here, see:
# https://stackoverflow.com/questions/57501259/coda-gelman-diag-error-in-chol-defaultw-the-leading-minor-of-order-nn-is
#gd_p_omeg_0 <- coda::gelman.diag(mc_p_0$Omega[[1]], multivariate = T)$psrf
#gd_p_omeg_0 <- "did not work"
# i'm not sure how much it matters if we're just gonna use rhat from MCMCvis

save(list = c('m_p_0', 'mc_p_0', 'ma_p_0', 'mp_p_0', 'mp_pp_0'), file = file.path(mod_dir, 'm_p_mod0.rda'))
save(list = c(
  'mf_p_0', 'm_ca_0', 'm_vp_0', 'mc_s_beta_0', 'mc_s_gamm_0', 'mc_s_omeg_0',
  #'gd_p_beta_0', 'gd_p_omeg_0', 'gd_p_gamm_0',
  'post_beta_0', 'post_gamm_0', 'post_omeg_0'), file = file.path(mod_dir, 'm_p_diag0.rda')
)

# fit
# presence/absence probit model
# using all 400-m data
m_p_1 <- Hmsc::Hmsc(
  Y = mat_t, XData = env1, studyDesign = stu1, XFormula = xf1, TrData = tr1, ranLevels = rL1,
  distr = 'probit', TrFormula = ~nativeStatusCode, XScale = F
)
m_p_1 <- Hmsc::sampleMcmc(m_p_1, thin = 10, samples = s0[1], transient = s0[2], nChains = s0[3], nParallel = s0[4])
mc_p_1 <- Hmsc::convertToCodaObject(m_p_1)
ma_p_1 <- Hmsc::computeAssociations(m_p_1, thin = 10)
mp_p_1 <- Hmsc::computePredictedValues(m_p_1)
mp_pp_1 <- predict(m_p_1)

mf_p_1 <- Hmsc::evaluateModelFit(hM = m_p_1, predY = mp_p_1)
m_vp_1 <- Hmsc::computeVariancePartitioning(m_p_1)
m_pv_1 <- Hmsc::computePredictedValues(m_p_1)
m_ca_1 <- Hmsc::computeAssociations(m_p_1)

mc_s_beta_1 <- MCMCvis::MCMCsummary(mc_p_1$Beta)
mc_s_gamm_1 <- MCMCvis::MCMCsummary(mc_p_1$Gamma)
mc_s_omeg_1 <- MCMCvis::MCMCsummary(mc_p_1$Omega[[1]])

post_beta_1 <- NEON1::posterior_from_coda(mc_p_1, 'Beta', 0.9, average = T, drop_ns = T)
post_gamm_1 <- NEON1::posterior_from_coda(mc_p_1, 'Gamma', 0.9, average = T, drop_ns = T)
post_omeg_1 <- NEON1::posterior_from_coda(mc_p_1, 'Omega', 0.9, average = T, drop_ns = T)

## Gelman's PSRF (Potential Scale Reduction Factor)
# Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations.
# Journal of Computational and Graphical Statistics, 7, 434-455.
#gd_p_beta_1 <- coda::gelman.diag(mc_p_1$Beta, multivariate = T)$psrf
#gd_p_gamm_1 <- coda::gelman.diag(mc_p_1$Gamma, multivariate = T)$psrf
#gd_p_omeg_1 <- coda::gelman.diag(mc_p_1$Omega[[1]], multivariate = T)$psrf
#gd_p_omeg_1 <- "did not work"

save(list = c('m_p_1', 'mc_p_1', 'ma_p_1', 'mp_p_1', 'mp_pp_1'), file = file.path(mod_dir, 'm_p_mod1.rda'))
save(list = c(
  'mf_p_1', 'm_ca_1', 'm_vp_1', 'mc_s_beta_1', 'mc_s_gamm_1', 'mc_s_omeg_1',
  #'gd_p_beta_1', 'gd_p_omeg_1', 'gd_p_gamm_1',
  'post_beta_1', 'post_gamm_1', 'post_omeg_1'), file = file.path(mod_dir, 'm_p_diag1.rda')
)

# Trace plots
MCMCvis::MCMCtrace(mc_p_0$Beta)
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_0_beta.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
MCMCvis::MCMCtrace(mc_p_0$Gamma)
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_0_gamma.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
#MCMCvis::MCMCtrace(mc_p_0$Omega[[1]])
#file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_0_omega.pdf'), overwrite = T)
#file.remove('MCMCtrace.pdf')

MCMCvis::MCMCtrace(mc_p_1$Beta)
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_1_beta.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
MCMCvis::MCMCtrace(mc_p_1$Gamma)
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_1_gamma.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
#MCMCvis::MCMCtrace(mc_p_1$Omega[[1]])
#file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_1_omega.pdf'), overwrite = T)
#file.remove('MCMCtrace.pdf')

# Rhat summaries
rhat_beta_0 <- mc_s_beta_0[which(mc_s_beta_0$Rhat > 1.01), ]
rhat_beta_0$term <- row.names(rhat_beta_0)
rhat_beta_0$species <- sapply(strsplit(row.names(rhat_beta_0), ' '), \(xx) xx[3])
rhat_beta_0$var <- sapply(strsplit(row.names(rhat_beta_0), ' '), \(xx) xx[1])
rhat_beta_0$var <- gsub('B\\[', '', rhat_beta_0$var)
rhat_beta_0$var <- gsub('cover_type', 'cover_type_', rhat_beta_0$var)
rhat_beta_0 <- data.frame(param = rep('beta0', nrow(rhat_beta_0)), term = rhat_beta_0$term, species = rhat_beta_0$species, var = rhat_beta_0$var, rhat_beta_0[, -which(colnames(rhat_beta_0) %in% c('term', 'var', 'species'))], row.names = NULL)
#rhat_beta_0 <- rhat_beta_0[order(rhat_beta_0$species, rhat_beta_0$var), ]
#write.csv(rhat_beta_0, file.path(mod_dir, 'rhat_beta_0.csv'), row.names = F)

rhat_beta_1 <- mc_s_beta_1[which(mc_s_beta_1$Rhat > 1.01), ]
rhat_beta_1$term <- row.names(rhat_beta_1)
rhat_beta_1$species <- sapply(strsplit(row.names(rhat_beta_1), ' '), \(xx) xx[3])
rhat_beta_1$var <- sapply(strsplit(row.names(rhat_beta_1), ' '), \(xx) xx[1])
rhat_beta_1$var <- gsub('B\\[', '', rhat_beta_1$var)
rhat_beta_1$var <- gsub('cover_type', 'cover_type_', rhat_beta_1$var)
rhat_beta_1 <- data.frame(param = rep('beta1', nrow(rhat_beta_1)), term = rhat_beta_1$term, species = rhat_beta_1$species, var = rhat_beta_1$var, rhat_beta_1[, -which(colnames(rhat_beta_1) %in% c('term', 'var', 'species'))], row.names = NULL)
#rhat_beta_1 <- rhat_beta_1[order(rhat_beta_1$species, rhat_beta_1$var), ]
#write.csv(rhat_beta_1, file.path(mod_dir, 'rhat_beta_1.csv'), row.names = F)

rhat_gamm_0 <- mc_s_gamm_0[which(mc_s_gamm_0$Rhat > 1.01), ]
rhat_gamm_0$term <- row.names(rhat_gamm_0)
rhat_gamm_0$species <- sapply(strsplit(row.names(rhat_gamm_0), ' '), \(xx) xx[3])
rhat_gamm_0$var <- sapply(strsplit(row.names(rhat_gamm_0), ' '), \(xx) xx[1])
rhat_gamm_0$var <- gsub('B\\[', '', rhat_gamm_0$var)
rhat_gamm_0$var <- gsub('cover_type', 'cover_type_', rhat_gamm_0$var)
rhat_gamm_0 <- data.frame(param = rep('gamm0', nrow(rhat_gamm_0)), term = rhat_gamm_0$term, species = rhat_gamm_0$species, var = rhat_gamm_0$var, rhat_gamm_0[, -which(colnames(rhat_gamm_0) %in% c('term', 'var', 'species'))], row.names = NULL)
#rhat_gamm_0 <- rhat_gamm_0[order(rhat_gamm_0$species, rhat_gamm_0$var), ]
#write.csv(rhat_gamm_0, file.path(mod_dir, 'rhat_gamm_0.csv'), row.names = F)

rhat_gamm_1 <- mc_s_gamm_1[which(mc_s_gamm_1$Rhat > 1.01), ]
rhat_gamm_1$term <- row.names(rhat_gamm_1)
rhat_gamm_1$species <- sapply(strsplit(row.names(rhat_gamm_1), ' '), \(xx) xx[3])
rhat_gamm_1$var <- sapply(strsplit(row.names(rhat_gamm_1), ' '), \(xx) xx[1])
rhat_gamm_1$var <- gsub('B\\[', '', rhat_gamm_1$var)
rhat_gamm_1$var <- gsub('cover_type', 'cover_type_', rhat_gamm_1$var)
rhat_gamm_1 <- data.frame(param = rep('gamm1', nrow(rhat_gamm_1)), term = rhat_gamm_1$term, species = rhat_gamm_1$species, var = rhat_gamm_1$var, rhat_gamm_1[, -which(colnames(rhat_gamm_1) %in% c('term', 'var', 'species'))], row.names = NULL)
#rhat_gamm_1 <- rhat_gamm_1[order(rhat_gamm_1$species, rhat_gamm_1$var), ]
#write.csv(rhat_gamm_1, file.path(mod_dir, 'rhat_gamm_1.csv'), row.names = F)

rhat_df <- rbind(rhat_beta_0, rhat_beta_1, rhat_gamm_0, rhat_gamm_1)
rhat_df <- rhat_df[order(rhat_df$param, rhat_df$species, rhat_df$var), ]
write.csv(rhat_df, file.path(mod_dir, 'rhat_eval.csv'), row.names = F)

rhat_omeg_0 <- mc_s_omeg_0[which(mc_s_omeg_0$Rhat > 1.01), ]
rhat_omeg_0$term <- row.names(rhat_omeg_0)
rhat_omeg_0$species <- sapply(strsplit(row.names(rhat_omeg_0), ' '), \(xx) xx[3])
rhat_omeg_0$var <- sapply(strsplit(row.names(rhat_omeg_0), ' '), \(xx) xx[1])
rhat_omeg_0$var <- gsub('Omega1\\[', '', rhat_omeg_0$var)
rhat_omeg_0 <- data.frame(param = rep('omeg0', nrow(rhat_omeg_0)), term = rhat_omeg_0$term, species = rhat_omeg_0$species, var = rhat_omeg_0$var, rhat_omeg_0[, -which(colnames(rhat_omeg_0) %in% c('term', 'var', 'species'))], row.names = NULL)

rhat_omeg_1 <- mc_s_omeg_1[which(mc_s_omeg_1$Rhat > 1.01), ]
rhat_omeg_1$term <- row.names(rhat_omeg_1)
rhat_omeg_1$species <- sapply(strsplit(row.names(rhat_omeg_1), ' '), \(xx) xx[3])
rhat_omeg_1$var <- sapply(strsplit(row.names(rhat_omeg_1), ' '), \(xx) xx[1])
rhat_omeg_1$var <- gsub('Omega1\\[', '', rhat_omeg_1$var)
rhat_omeg_1 <- data.frame(param = rep('omeg1', nrow(rhat_omeg_1)), term = rhat_omeg_1$term, species = rhat_omeg_1$species, var = rhat_omeg_1$var, rhat_omeg_1[, -which(colnames(rhat_omeg_1) %in% c('term', 'var', 'species'))], row.names = NULL)

rhat_omeg <- rbind(rhat_omeg_0, rhat_omeg_1)
rhat_omeg <- rhat_omeg[order(rhat_omeg$param, rhat_omeg$species, rhat_omeg$var), ]
write.csv(rhat_omeg, file.path(mod_dir, 'rhat_omeg.csv'), row.names = F)

# omega traces
MCMCvis::MCMCtrace(mc_p_0$Omega[[1]])
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_0_omega.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
MCMCvis::MCMCtrace(mc_p_1$Omega[[1]])
file.copy('MCMCtrace.pdf', file.path(mod_dir, 'trace_1_omega.pdf'), overwrite = T)
file.remove('MCMCtrace.pdf')
