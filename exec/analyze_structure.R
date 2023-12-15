# BEM May 2023
library(NEON1)
library(neonUtilities)
library(dplyr)
import::from(magrittr, "%>%")
library(ggplot2)

# Ecosystem structure DP3.30015.001
# Non-herbaceous perennial vegetation structure DP1.10045.001
# Plant presence and percent cover DP1.0058.001
# Vegetation structure DP1.10098.001
struct_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_struct-plant.zip'

#struct_file <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/clim_1.tif'

# Structure

struct_tbl <- neonUtilities::stackByTable(struct_zip, savepath = 'envt')
# `vst_apparentindividuals` : Biomass and productivity measurements of woody individuals
# `vst_non-woody` : Biomass and productivity measurements of non-herbaceous perennial plants
# `vst_shrubgroup` : Biomass and productivity measurements of groups of shrubs

struct <- struct_tbl[['vst_apparentindividual']] %>%
  # drop some unnecessary columns
  select(-c(uid, namedLocation, eventID, domainID, siteID, dendrometerInstallationDate, initialBandStemDiameter, initialDendrometerGap, dendrometerHeight, dendrometerGap, dendrometerCondition, bandStemDiameter, remarks, recordedBy, measuredBy, dataEntryRecordID, maxBaseCrownDiameter, ninetyBaseCrownDiameter, initialGapMeasurementDate, dataQF, publicationDate, release)) %>%
  # drop some growth forms
  filter(!growthForm %in% c('liana', 'small shrub', 'single shrub')) %>%
  # add variable for date
  mutate(date = as.Date(startDateTime)) %>%
  filter(date >= min(date[which(siteID == 'PUUM')])) %>%
  # add variable for month/year ('m_y')
  #mutate(m_y = paste(format(date, '%m'), format(date, '%Y'), sep = '_')) %>%
  mutate(m_y = lubridate::floor_date(date, 'month')) %>%
  group_by(siteID, m_y) %>%
  summarize(ppt = sum(secPrecipBulk)) %>%
  ungroup()

map <- struct_tbl[['vst_mappingandtagging']] %>%
  # species/other information is in the mapping and tagging table
  # drop some unnecessary columns
  select(-c(uid, namedLocation, eventID, domainID, siteID, remarks, recordedBy, measuredBy, dataQF, publicationDate, publicationDate, release, otherTagID, otherTagOrg, samplingProtocolVersion, identificationReferences, morphospeciesID, morphospeciesIDRemarks, identificationQualifier)) %>%
  # only two stems had previously tagged IDS
  select(-c(previouslyTaggedAs, recordType)) %>%
  # not concerned about supporting individual currently, b/c 84/4057 stems have support
  select(-supportingStemIndividualID) %>%
  # dont need plot/plotIDs if we have individual ID
  select(-c(plotID, subplotID, nestedSubplotID, pointID)) %>%
  # 94% of stems ID'd to species
  select(-c(taxonRank, scientificName)) %>%
  # only 1 stem duplicated in data set: 00152 has entry with map, and one without
  subset(!(individualID == 'NEON.PLA.D20.PUUM.00152' & is.na(stemDistance))) %>%
  # range(table(map$date, map$individualID)) shows only 1 individual measured per date
  select(-date)

# convert azimuth map to x/y?
# in `geoNEON` package

# add map to struct df
struct <- dplyr::left_join(struct, map, by = 'individualID')

ggplot(data = struct, aes(x = taxonID, y = height)) +
  geom_boxplot() +

  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),
    axis.text.y = element_text(color = 'black')
  )
