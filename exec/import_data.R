# BEM July 2024
# updated August 2024

# Last run: "2024-08-08 10:24:52 HST"

# Purpose: Import NEON data relevant to all the current analyses

# TODO: add data:
#           DP3.30024.001 elevation-lidar
#           DP1.10026.001 traits
# Trait data seems complicated, unsure how to connect the samples to each other

# Methods:
#     Only uses rows identified at least to genus (excludes 211 species)
#     Treats both sp. and spp. as individual species
#     Treats unknown/uncertain native status as native.
#         Peperomia, Rhynchospora, and Carex are the big offenders here.

library(NEON1)
library(neonUtilities)
#library(neonOS)
library(dplyr)
library(terra)

#e1_dir <- '../spatial/AOP'
e1_dir <- '/media/bem/data/NEON/spatial/AOP'
aop_dir <- 'DP3.30024.001/neon-aop-products/2020/FullSite/D20/2020_PUUM_2/L3/DiscreteLidar/DTMGtif'
aop_fl <- '/media/bem/data/NEON/spatial/aop_elev_mosaic.tif'

### Load data products
d0 <- neonUtilities::loadByProduct('DP1.10058.001', 'PUUM', include.provisional = T, check.size = F)
d1 <- neonUtilities::loadByProduct('DP1.10098.001', 'PUUM', include.provisional = T, check.size = F)
d2 <- neonUtilities::loadByProduct('DP1.10047.001', 'PUUM', include.provisional = T, check.size = F)
d3 <- neonUtilities::loadByProduct('DP1.10086.001', 'PUUM', include.provisional = T, check.size = F)
d4 <- neonUtilities::loadByProduct('DP1.10026.001', 'PUUM', include.provisional = T, check.size = F)
# DP3.30024.001 is 3.36 GB downloaded 2024-08-08 12:31 PM HST
#e1 <- neonUtilities::byFileAOP('DP3.30024.001', 'PUUM', 2020, T, F, e1_dir)

### Elevation DTM
if (!file.exists(aop_fl)) {
  tif_fls <- list.files(file.path(e1_dir, aop_dir), full.names = T)
  # mosaic all the files
  dtm <- terra::rast()
  for (i in seq_along(tif_fls)) {
    ii <- tif_fls[i]
    cat(paste0('\rprogress: ', round(i / length(tif_fls) * 100), rep(' ', 10)))
    irast <- terra::rast(ii)

    if (i == 1) {
      dtm <- irast
    } else {
      dtm <- terra::mosaic(dtm, irast)
    }

    if (i %% 100 == 0) {
      terra::writeRaster(dtm, aop_fl, overwrite = T, filetype = 'GTiff')
    }

  }

  terra::writeRaster(dtm, aop_fl, overwrite = T, filetype = 'GTiff')
} else {
  dtm <- terra::rast(aop_fl)
}

### Meta-data

# Species meta-data
# Combines both 1m and 10m+ data, structure data
spp <- d0$div_1m2Data |>
  filter(divDataType == 'plantSpecies') |>
  select(scientificName, taxonID, taxonRank, family, nativeStatusCode) |>
  bind_rows(select(.data = d0$div_10m2Data100m2Data, taxonID, scientificName, taxonRank, family, nativeStatusCode)) |>
  bind_rows(select(.data = d1$vst_mappingandtagging, taxonID, scientificName, taxonRank)) |>
  mutate(binomialName = gsub('\\.', '', sapply(strsplit(scientificName, ' '), \(xx) paste(xx[1], xx[2], sep = '_')))) |>
  distinct()
# Plot meta-data
met <- d0$div_1m2Data |>
  select(plotID, decimalLatitude, decimalLongitude, elevation, plotType) |>
  distinct()

# neonOS::joinTableNEON()

### Plant diversity/percent cover data
# Plant cover data - 1m
div_one <- d0$div_1m2Data |>
  filter(
    divDataType == 'plantSpecies',
    taxonRank %in% c('genus', 'species', 'subspecies', 'variety')
    # drops 1 family + 210 kingdom = 211 rows, mostly not-identified
  ) |>
  select(plotID, subplotID, endDate, scientificName, percentCover) |>
  distinct()
# Plant cover data - 10m
div_ten <- d0$div_10m2Data100m2Data |>
  filter(taxonRank %in% c('genus', 'species', 'subspecies', 'variety')) |>
  select(plotID, subplotID, endDate, taxonID) |>
  distinct()

### Plant structure data
str <- d1$vst_apparentindividual |>
  select(date, plotID, subplotID, individualID, growthForm, plantStatus, stemDiameter, basalStemDiameter, height) |>
  distinct()
# To derive the mapped location of each individual tree in vst_mappingandtagging, download the R
# geoNEON package (https://github.com/NEONScience/NEON‐geolocation) and use the getLocTOS() function.
str <- d1$vst_mappingandtagging |>
  select(individualID, scientificName) |>
  distinct() |>
  # there are 7 individualIDs which have two values for 'species', just take the last entry
  group_by(individualID) |>
  summarise(scientificName = last(scientificName), .groups = 'drop') |>
  left_join(str, by = 'individualID')

### Soil data
soil_init <- d2$spc_biogeochem |>
  select(
    plotID, horizonName, biogeoCenterDepth, # unique identifiers + depth
    ececCecd33, # cation exchange capacity
    phCacl2, phH2o, acidity, ec12pre, # pH, pH, "acidity", EC
    carbonTot, nitrogenTot, ctonRatio, estimatedOC, # total C, total N, C:N, organic C
    pOxalate, MehlichIIITotP, OlsenPExtractable # phosphorus, 3 methods
  )
soil_peri <- d3$sls_soilChemistry |>
  # d3$ntr_externalLab and d3$sls_soilpH is detailed BGC probably not relevant to current study
  select(plotID, collectDate, d15N, organicd13C, nitrogenPercent, organicCPercent) |>
  group_by(plotID, collectDate) |>
  summarise(across(everything(), ~ mean(.x, na.rm = T)), .groups = 'drop')

### Trait data
# cfc_fieldData connects all the metadata
trait_meta <- d4$cfc_fieldData |>
  select(plotID, collectDate, subplotID, sampleID, chlorophyllSampleID, individualID, scientificName)

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
usethis::use_data(str, overwrite = T)
usethis::use_data(soil_init, overwrite = T)
usethis::use_data(soil_peri, overwrite = T)
