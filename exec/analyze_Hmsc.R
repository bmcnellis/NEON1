# BEM 29 July 2024

set.seed(1)

# TODO: need to normalize/scale the cover values
#       break this into prepare_Hmsc and analyze_Hmsc after first model draft
#       add computeVariancePartioning to the workflow
#       should some species be excluded based on ESS?

# Caveats:
#     * soil age is an ordered factor
#     * taxon excluded: 1 family + 210 kingdom = 211 rows, mostly not-identified
#     * two enclosures (Kata-Stein & Unencumbered Kulani Pasture) had NA for completion date, this
#       was set to 2014 which is the age of the youngest bordering enclosure
#     * lava flow age is the median of the age range published on USGS maps
#     * cpai, elev, TWI, age_median are centered and scaled
#     * time_since_fence NOT centered or scaled

### Libraries
library(NEON1)
library(dplyr)
library(Hmsc)

### Directories
#data_dir <- '/media/bem/data/NEON'
data_dir <- 'C:/Users/BrandonMcNellis/Documents/NEON_data'
resid_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/reports'
#mod_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/model_results'
mod_dir <- '/media/bem/data/NEON/results/model_results'

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
  select(plotDate, plotID, endDate, binomialName, percentCover, ageFactor, age_median, nativeFactor, pai, fcover, fipar) |>
  mutate(logicalCover = ifelse(percentCover < 0.0001, FALSE, TRUE)) |>
  mutate(fractionCover = percentCover / 100) |>
  left_join(NEON1::fence, by = 'plotID') |>
  left_join(NEON1::met_fence[, c('enclosure_name', 'enclosure_completed', 'enclosure_ung_free')], by = 'enclosure_name') |>
  # set two NA enclosures to 2014, which is the youngest neighboring enclosure
  mutate(enclosure_completed = ifelse(is.na(enclosure_completed), 2014, enclosure_completed)) |>
  mutate(time_since_fence = 2024 - enclosure_completed) |>
  left_join(NEON1::topo[, c('plotID', 'elev', 'TWI')], by = 'plotID') |>
  # center and scale pai, elev, TWI, age_median
  # the as.numeric is because Hmsc doesnt like the scale() Value type
  mutate(pai = as.numeric(scale(pai)), elev = as.numeric(scale(elev)), TWI = as.numeric(scale(TWI)), age_median = as.numeric(scale(age_median))) |>
  distinct()
rm(inv_nat)

# create species matrix for fractional cover, plot-date as row with species as col
mat_0 <- df0 |>
  select(c(plotDate, binomialName, fractionCover)) |>
  group_by(plotDate, binomialName) |>
  tidyr::pivot_wider(names_from = binomialName, values_from = fractionCover) |>
  arrange(plotDate) |>
  tibble::column_to_rownames('plotDate') |>
  as.matrix()
# create species matrix for presence-absence, plot-date as row with species as col
mat_p <- mat_0 |>
  apply(c(1, 2), \(xx) ifelse(xx > 0.0001, 1, 0))
mat_f <- mat_0 |>
  apply(c(1, 2), \(xx) ifelse(xx < 0.001, NA, xx)) |>
  apply(c(1, 2), \(xx) ifelse(is.na(xx), NA, NEON1::logit_Warton(xx)))
# create environmental matrix
env1 <- df0 |>
  select(c(plotDate, age_median, pai, time_since_fence, elev, TWI)) |>
  arrange(plotDate) |>
  distinct() |>
  as.data.frame()
# create study-design matrix
stu1 <- df0 |>
  select(c(plotDate, plotID, endDate)) |>
  arrange(plotDate) |>
  mutate(dateFactor = as.factor(as.character(endDate))) |>
  mutate(plotID = as.factor(plotID)) |>
  mutate(plotDate = as.factor(plotDate)) |>
  select(plotDate, plotID, dateFactor) |>
  distinct() |>
  as.data.frame()

# check matrix/df alignment and modify for modelling
stopifnot(
  nrow(mat_f) == nrow(mat_p),
  nrow(mat_f) == nrow(env1),
  nrow(mat_f) == nrow(stu1),
  identical(row.names(mat_f), row.names(mat_p)),
  identical(row.names(mat_f), env1$plotDate),
  identical(row.names(mat_f), as.character(stu1$plotDate))
)
mat_f <- mat_f |>
  `rownames<-`(NULL)
env1 <- env1 |>
  select(-plotDate)

# check distribution of effects and covariates
hist(env1$pai, breaks = 30)
# pai looks bimodal
#hist(as.numeric(env1$ageFactor), breaks = 30)
hist(env1$age_median, breaks = 30)
# soil age looks lognormal
hist(env1$time_since_fence, breaks = 30)
# not very normal
hist(env1$elev, breaks = 30)
# super bimodal
hist(env1$TWI, breaks = 30)
# mostly normal

# modelling
# kable(VP$R2T$Beta), where VP is the result of computeVariancePartitioning

# fixed-effects
xf1 <- as.formula(~ age_median + pai + time_since_fence + elev*TWI)
# random effects
#rL <- Hmsc::HmscRandomLevel(units = stu1$plotID)
# random effect at sampling unit estimates species-species associations
rL <- Hmsc::HmscRandomLevel(units = stu1$plotDate)
# constrain latent factors to 2 levels
rL$nfMin <- 2
rL$nfMax <- 2

# fit
# presence/absence probit model
if (!file.exists(file.path(mod_dir, 'm_p_mod.rda'))) {
  m_p <- Hmsc::Hmsc(Y = mat_p, XData = env1, XFormula = xf1, studyDesign = stu1, ranLevels = list('plotDate' = rL), distr = 'probit')
  m_p <- Hmsc::sampleMcmc(m_p, thin = 10, samples = 5000, transient = 1000, nChains = 4, nParallel = 4)
  mc_p <- Hmsc::convertToCodaObject(m_p)
  ma_p <- Hmsc::computeAssociations(m_p, thin = 10)
  # requires that the normality of the variables should be normally distributed
  mp_p <- Hmsc::computePredictedValues(m_p)
  save(list = c('m_p', 'mc_p', 'ma_p', 'mp_p'), file = file.path(mod_dir, 'm_p_mod.rda'))
} else {
  load(file.path(mod_dir, 'm_p_mod.rda'))
}

# logit proportion cover/NA normal model

# evaluate
# Beta = species niches, Gamma = traits on niches, Omega = residual species associations, rho = phylogenetic signal

# presence/absence probit model
# Beta and Omega are in this model, could maybe figure out traits later
if (!file.exists(file.path(mod_dir, 'm_p_diag.rda'))) {
  mf_p <- Hmsc::evaluateModelFit(hM = m_p, predY = mp_p)
  es_p_beta <- coda::effectiveSize(mc_p$Beta)
  #es_p_gamm <- coda::effectiveSize(mc_p$Gamma)
  es_p_omeg <- coda::effectiveSize(mc_p$Omega[[1]])
  gd_p_beta <- coda::gelman.diag(mc_p$Beta, multivariate = F)$psrf
  #gd_p_gamm <- coda::gelman.diag(mc_p$Gamma, multivariate = F)$psrf
  gd_p_omeg <- coda::gelman.diag(mc_p$Omega[[1]], multivariate = F)$psrf
  mpe_p_beta <- Hmsc::getPostEstimate(m_p, parName = 'Beta')
  #mpe_p_gamm <- Hmsc::getPostEstimate(m_p, parName = 'Gamma')
  mpe_p_omeg <- Hmsc::getPostEstimate(m_p, parName = 'Omega')
  save(list = c('mf_p', 'es_p_beta', 'es_p_omeg', 'gd_p_beta', 'gd_p_omeg', 'mpe_p_beta', 'mpe_p_omeg'), file = file.path(mod_dir, 'm_p_diag.rda'))
} else {
  load(file.path(mod_dir, 'm_p_diag.rda'))
}

hist(mf_p$R2, xlim = c(0,1), main = paste0("Mean = ", round(mean(mf_p$R2), 2)))
hist(es_p_beta, breaks = 30, main = 'ess:Beta, model p')
#hist(es_p_gamm, breaks = 30, main = 'ess:Gamma, model p')
hist(es_p_omeg, breaks = 30, main = 'ess:Omega, model p')
hist(gd_p_beta, breaks = 30, main = 'psrf:Beta, model p')
#hist(gd_p_gamm, breaks = 30, main = 'psrf:Gamma, model p')
hist(gd_p_omeg, breaks = 30, main = 'psrf:Omega, model p')

# cross-validation
#ml_pa <- Hmsc::createPartition(ml, nfolds = 2, column = 'plotDate')
# commenting out for now, should probably compare it somehow
#ml_xv <- Hmsc::computePredictedValues(ml, partition = ml_pa, nParallel = 2)
#Hmsc::evaluateModelFit(hM = ml, predY = ml_xv)

# can we add a spatially explicit model to incorporate subplots?

# figures


ml_pb <- getPostEstimate(ml, parName = 'Beta')
Hmsc::plotBeta(ml, post = ml_pb, param = 'Support', supportLevel = 0.95)


# i think this is now ma_p
sl <- 0.5
cp <- ((ml_c[[1]]$support > sl) + (ml_c[[1]]$support < (1 - sl)) > 0) * ml_c[[1]]$mean
corrplot::corrplot(
  cp, method = 'color', col = colorRampPalette(c('blue', 'white', 'red'))(200),
  title = paste('random effect level:', ml$rLNames[1]), mar = c(0, 0, 1, 0)
)
