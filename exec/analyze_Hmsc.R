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

stopifnot(nrow(env1) == nrow(mat0), nrow(env1) == nrow(stu1))
stopifnot(identical(row.names(mat0), env1$plotDate))
stopifnot(identical(row.names(mat0), as.character(stu1$plotDate)))
row.names(mat0) <- NULL
env1 <- env1[, -1]
stu1$plotID <- as.factor(stu1$plotID)

# modelling
# split up for hurdle
matL <- 1 * (mat0 > 0.0001)
matP <- apply(mat0, c(1, 2), \(xx) ifelse(xx < 0.000001, NA, xx))
xf1 <- as.formula(~pai)
rL <- Hmsc::HmscRandomLevel(units = stu1$plotID)

# testing
matL <- matL[, sample(ncol(matL), size = 10, replace = F)]
t0 <- rowSums(matL) > 0
matL <- matL[t0, ]
env1 <- env1[t0, ]
stu1 <- stu1[t0, ]
stu1$plotDate <- as.factor(as.character(stu1$plotDate))
# regular random effect for plot
#rL <- Hmsc::HmscRandomLevel(units = stu1$plotID)
# random effect at sampling unit estimates species-species associations
rL <- Hmsc::HmscRandomLevel(units = stu1$plotDate)
# constraints latent factor number
#rL$nfMin <- 2
#rL$nfMax <- 2

# presence/absence probit model
ml <- Hmsc::Hmsc(Y = matL, XData = env1, XFormula = xf1, studyDesign = stu1, ranLevels = list('plotDate' = rL), distr = 'probit')
ml <- Hmsc::sampleMcmc(ml, thin = 10, samples = 2000, transient = 500, nChains = 4, nParallel = 2)
ml_p <- Hmsc::convertToCodaObject(ml)
ml_c <- Hmsc::computeAssociations(ml, thin = 10)


coda::effectiveSize(ml_p$Beta)
coda::gelman.diag(ml_p$Beta, multivariate = F)$psrf
# 'Potential scale reduction factor', upper C.I. should be close to 1
# requires that the normality of the variables should be normally distributed
ml_pd <- Hmsc::computePredictedValues(ml, expected = F)
Hmsc::evaluateModelFit(hM = ml, predY = ml_pd)
ml_pa <- Hmsc::createPartition(ml, nfolds = 2, column = 'plotDate')
# commenting out for now, should probably compare it somehow
#ml_xv <- Hmsc::computePredictedValues(ml, partition = ml_pa, nParallel = 2)
#Hmsc::evaluateModelFit(hM = ml, predY = ml_xv)

# can we add a spatially explicit model to incorporate subplots?

# figures
hist(coda::effectiveSize(ml$Beta), main= 'ess(beta)')
hist(coda::gelman.diag(ml$Beta, multivariate = F)$psrf, main= 'psrf(beta)')

ml_pb <- getPostEstimate(ml, parName = 'Beta')
Hmsc::plotBeta(ml, post = ml_pb, param = 'Support', supportLevel = 0.95)


sl <- 0.5
cp <- ((ml_c[[1]]$support > sl) + (ml_c[[1]]$support < (1 - sl)) > 0) * ml_c[[1]]$mean
corrplot::corrplot(
  cp, method = 'color', col = colorRampPalette(c('blue', 'white', 'red'))(200),
  title = paste('random effect level:', ml$rLNames[1]), mar = c(0, 0, 1, 0)
)
