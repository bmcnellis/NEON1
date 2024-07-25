# BEM July 2024

# Purpose: Import NEON data relevant to all the current analyses

# TODO: add the beetle and structure data here

# Methods:
#     Only uses rows identified at least to genus (excludes 211 species)
#     Treats both sp. and spp. as individual species
#     Treats unknown/uncertain native status as native.
#         Peperomia, Rhynchospora, and Carex are the big offenders here.

library(NEON1)
library(neonUtilities)
library(dplyr)

### Plant diversity/percent cover data
d0 <- neonUtilities::loadByProduct('DP1.10058.001', 'PUUM', include.provisional = T, check.size = F)
# "2024-07-24 17:03:49 HST"

# Plant cover data - 1m
div_one <- d0$div_1m2Data |>
  filter(
    divDataType == 'plantSpecies',
    taxonRank %in% c('genus', 'species', 'subspecies', 'variety')
    # drops 1 family + 210 kingdom = 211 rows, mostly not-identified
  ) |>
  select(-c(
    uid, namedLocation, domainID, siteID, geodeticDatum, coordinateUncertainty,
    elevationUncertainty, nlcdClass, boutNumber, samplingProtocolVersion, eventID,
    divDataType, targetTaxaPresent, otherVariablesPresent, otherVariablesPresent,
    taxonRank, family, plotType, decimalLatitude, decimalLongitude, elevation, scientificName,
    identificationQualifier, taxonIDRemarks, morphospeciesID, morphospeciesIDRemarks,
    identificationReferences, identificationHistoryID, otherVariables, heightPlantSpecies,
    remarks, measuredBy, recordedBy, samplingImpractical, samplingImpracticalRemarks,
    biophysicalCriteria, publicationDate, release, nativeStatusCode
  )) |>
  distinct()

# Plant cover data - 10m
div_ten <- d0$div_10m2Data100m2Data |>
  filter(taxonRank %in% c('genus', 'species', 'subspecies', 'variety')) |>
  select(plotID, subplotID, endDate, taxonID) |>
  distinct()

# Plot meta-data
met <- d0$div_1m2Data |>
  select(plotID, decimalLatitude, decimalLongitude, elevation, plotType) |>
  distinct()

# Species meta-data
# Combines both 1m and 10m+ data
spp <- d0$div_1m2Data |>
  filter(divDataType == 'plantSpecies') |>
  select(taxonID, scientificName, taxonRank, family, nativeStatusCode) |>
  bind_rows(select(.data = d0$div_10m2Data100m2Data, taxonID, scientificName, taxonRank, family, nativeStatusCode)) |>
  mutate(
    binomialName = sapply(strsplit(scientificName, ' '), \(xx) paste(xx[1], xx[2])),
    nativeSimple = case_when(nativeStatusCode == 'I' ~ 'I', is.na(nativeStatusCode) ~ NA, .default = 'N')
  ) |>
  distinct()

### dhp data

# PAI    : Plant area index
# GAI    : Green area index
# FIPAR  : Fraction  of intercepted photosynthetically active radiation
# FCOVER : Fraction of vegetation cover
# PAIe   : Effective plant area index

# Upward-facing images represent PAI, as the image classification is sensitive to all canopy elements
# Downward-facing images represent GAI, as the image classification is sensitive to green elements
# PAI = GAI hereafter

# Brown et al. (2023) suggest Hinge-method is best. Using PAI versus PAIe for unknown reasons.

dhp <- read.csv('../data/neon_dhp_puum.csv') |>
  rename('plotID' = plot) |>
  mutate(
    date = as.Date(strptime(date, '%d/%m/%Y', tz = 'HST')),
    pai = NEON1::sph(pai_up)[, 1],
    pai_err = NEON1::sph(pai_up)[, 2],
    fcover = NEON1::sph(fcover_up)[, 1],
    fcover_err = NEON1::sph(fcover_up)[, 2],
    fipar = NEON1::sph(fipar_up)[, 1],
    fipar_err = NEON1::sph(fipar_up)[, 2],
  ) |>
  select(c(plotID, date, pai, pai_err, fcover, fcover_err, fipar, fipar_err))


# what do the regularly-censused plots look like?
plot(dhp$date[which(dhp$plotID == 'PUUM_031')], dhp$pai[which(dhp$plotID == 'PUUM_031')], xlab = 'Date', ylab = 'PAI - 31')
plot(dhp$date[which(dhp$plotID == 'PUUM_039')], dhp$pai[which(dhp$plotID == 'PUUM_039')], xlab = 'Date', ylab = 'PAI - 39')
plot(dhp$date[which(dhp$plotID == 'PUUM_041')], dhp$pai[which(dhp$plotID == 'PUUM_041')], xlab = 'Date', ylab = 'PAI - 41')


# reduce down to one measurement per plot
dhp <- dhp |>
  group_by(plotID) |>
  summarise(across(-date, ~ mean(.x, na.rm = T)), .groups = 'drop')

# add measurements for 032, 034, 042
dhp <- rbind(dhp, data.frame(plotID = 'PUUM_032', within(dhp[which(dhp$plotID == 'PUUM_039'), ], rm(plotID))))
dhp <- rbind(dhp, data.frame(plotID = 'PUUM_034', within(dhp[which(dhp$plotID == 'PUUM_031'), ], rm(plotID))))
dhp <- rbind(dhp, data.frame(plotID = 'PUUM_042', within(dhp[which(dhp$plotID == 'PUUM_041'), ], rm(plotID))))

### save it all
usethis::use_data(div_one, overwrite = T)
usethis::use_data(div_ten, overwrite = T)
usethis::use_data(met, overwrite = T)
usethis::use_data(spp, overwrite = T)
usethis::use_data(dhp, overwrite = T)
