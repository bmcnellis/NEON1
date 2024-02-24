# NEON Ground beetles sampled from pitfall traps
# DP1.10022.001 @ PUUM

# The protocol dictates that each trap is collected once per bout (one expected record per trapID per plotID per
# collectDate in bet_fielddata). A record from bet_fielddata may have zero (if no sample collected) or multiple
# child records in bet_sorting depending on number of taxa contained in the sampleID. A record from from bet_sorting may
# have zero (if no contents of the subsampleID pinned) or multiple child records in bet_parataxonomistID depending
# on the number of individuals selected for pinning from each subsampleID. A record in bet_archivepooling may correspond
# to one or more records in bet_sorting, because multiple subsampleIDs may have been pooled into a single archiveVial.
# Each record in bet_parataxonomistID should have zero or one corresponding records in bet_expertTaxonomistIDProcessed,
# depending on whether that individualID was selected for professional identification. Each record in bet_parataxonomistID
# should also have zero or one corresponding records in bet_expertTaxonomistIDRaw. All beetles must be sorted prior to pinning,
# so the total number of beetles collected can be calculated as the sum of individualCount in bet_sorting, though further
# identifications may be updated based on the downstream workflow. Taxonomic IDs of species of concern have been 'fuzzed';
# see data package readme files for more information. If taxonomic determinations have been updated for any records in the
# table bet_archivepooling, bet_expertTaxonomistIDProcessed, bet_expertTaxonomistIDRaw, bet_parataxonomistID, bet_sorting,
# past determinations are archived in the bet_identificationHistory table, where the archived determinations are linked to
# current records using identificationHistoryID.

# Citation:
#
# NEON (National Ecological Observatory Network). Ground beetles sampled from pitfall traps (DP1.10022.001),
# RELEASE-2024. https://doi.org/10.48443/rcxn-t544. Dataset accessed from
# https://data.neonscience.org/data-products/DP1.10022.001/RELEASE-2024 on February 8, 2024.

# One sample unit is a trap ID x plotID x collectDate, or plot_trap_date

# Libraries and setup
library(neonUtilities)
library(dplyr)

# Pull relevant data frames from downloaded NEON data
bl <- neonUtilities::loadByProduct('DP1.10022.001', 'PUUM', include.provisional = F)
Y
bl_f <- bl$bet_fielddata
bl_f <- bl_f[which(bl_f$sampleCollected == 'Y'), ]
bl_p <- bl$bet_parataxonomistID
bl_x <- bl$bet_expertTaxonomistIDProcessed
bl_s <- bl$bet_sorting

# 'individualID' connects `bet_paratxonomistID` and `bet_expertTaxonomistID`
key_cols_1 <- c('namedLocation', 'domainID', 'siteID', 'plotID', 'setDate', 'collectDate', 'individualID')
# 'subsampleID' connects the identification data with the sorted sample data
key_cols_2 <- c('namedLocation', 'domainID', 'siteID', 'plotID', 'setDate', 'collectDate', 'subsampleID')
# 'sampleID' connects the field data to the sorted sample data
key_cols_3 <- c('namedLocation', 'domainID', 'siteID', 'plotID', 'setDate', 'collectDate', 'sampleID')

# Combine the multiple data streams into one sparse dataframe
# Combine `bet_paratxonomistID` and `bet_expertTaxonomistID`
# All `individualID` in bet_expertTaxonomistID are in bet_parataxonomistID,
bl_id <- dplyr::left_join(bl_p, bl_x, by = key_cols_1, suffix = c('_para', '_expe'))

# Combine `bet_paratxonomistID`, `bet_expertTaxonomistID`, and `bet_sorting`
# All `subsampleCode` in bet_parataxonomistID are in bet_sorting, reverse is not true
bl_sd <- dplyr::left_join(bl_s, bl_id, by = key_cols_2, suffix = c('_sort', '_id'))

# Combine `bet_paratxonomistID`, `bet_expertTaxonomistID`, `bet_sorting`, and `bet_fielddata`
bl_all <- dplyr::left_join(bl_sd, bl_f, by = key_cols_3, suffix = c('_sortID', '_field'))
bl_all <- bl_all[, order(colnames(bl_all))]

# Remove extraneous columns
# Remove extraneous taxonomic information
bl_all <- bl_all[, -which(colnames(bl_all) %in% c(
  'class', 'family', 'genus', 'kingdom', 'order', 'phylum', 'specificEpithet', 'subfamily', 'subgenus', 'tribe',
  'infraspecificEpithet'
))]
# Remove extraneous sample identification
# This assumes that the expert identification supersedes the parataxonomist identification
bl_all$scientificName_coal <- dplyr::coalesce(bl_all$scientificName_expe, bl_all$scientificName_para)
bl_all$scientificName_coal <- ifelse(is.na(bl_all$scientificName_expe), bl_all$scientificName_coal, bl_all$scientificName_expe)
# 'scientificName' from `bet_sorting` is the same as from `scientificName_para`
bl_all$scientificName <- bl_all$scientificName_coal
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('scientificName_coal', 'scientificName_para', 'scientificName_expe'))]
# Simliar to identification, 'nativeStatusCode' is the same as 'nativeStatusCode_para'
bl_all$nativeStatusCode <- dplyr::coalesce(bl_all$nativeStatusCode_expe, bl_all$nativeStatusCode_para)
bl_all$nativeStatusCode <- ifelse(is.na(bl_all$nativeStatusCode_expe), bl_all$nativeStatusCode, bl_all$nativeStatusCode_expe)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('nativeStatusCode_para', 'nativeStatusCode_expe'))]
bl_all <- bl_all[, -which(colnames(bl_all) %in% c(
  'identifiedBy', 'identifiedBy_expe', 'identifiedBy_para', 'identifiedDate', 'identifiedDate_para', 'identifiedDate_expe',
  'identificationHistoryID', 'identificationHistoryID_expe', 'identificationHistoryID_para',
  'identificationQualifier', 'identificationQualifier_expe', 'identificationQualifier_para',
  'identificationReferences', 'identificationReferences_expe', 'identificationReferences_para'
))]
# Simliar to identification, 'morphospeciesID' is the same as 'morphospeciesID_para'
bl_all$morphospeciesID <- dplyr::coalesce(bl_all$morphospeciesID_expe, bl_all$morphospeciesID_para)
bl_all$morphospeciesID <- ifelse(is.na(bl_all$morphospeciesID_expe), bl_all$morphospeciesID, bl_all$morphospeciesID_expe)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('morphospeciesID_para', 'morphospeciesID_expe'))]
# Simliar to identification, 'scientificNameAuthorship' is the same as 'scientificNameAuthorship_para'
bl_all$scientificNameAuthorship <- dplyr::coalesce(bl_all$scientificNameAuthorship_expe, bl_all$scientificNameAuthorship_para)
bl_all$scientificNameAuthorship <- ifelse(is.na(bl_all$scientificNameAuthorship_expe), bl_all$scientificNameAuthorship, bl_all$scientificNameAuthorship_expe)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('scientificNameAuthorship_para', 'scientificNameAuthorship_expe'))]
# Simliar to identification, 'taxonID' is the same as 'taxonID_para'
bl_all$taxonID <- dplyr::coalesce(bl_all$taxonID_expe, bl_all$taxonID_para)
bl_all$taxonID <- ifelse(is.na(bl_all$taxonID_expe), bl_all$taxonID, bl_all$taxonID_expe)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('taxonID_para', 'taxonID_expe'))]
# Simliar to identification, 'taxonRank' is the same as 'taxonRank_para'
bl_all$taxonRank <- dplyr::coalesce(bl_all$taxonRank_expe, bl_all$taxonRank_para)
bl_all$taxonRank <- ifelse(is.na(bl_all$taxonRank_expe), bl_all$taxonRank, bl_all$taxonRank_expe)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('taxonRank_para', 'taxonRank_expe'))]

# trapID taken from sort table, id table sometimes missing
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('trapID_id', 'trapID_sort'))]

# sampleCondition was coalesced because all but 3 entries had only 1 non-NA entry
bl_all$sampleCondition_para <- ifelse(bl_all$sampleCondition_para == 'OK', NA, bl_all$sampleCondition_para)
bl_all$sampleCondition_field <- ifelse(bl_all$sampleCondition_field == 'OK', NA, bl_all$sampleCondition_field)
bl_all$sampleCondition_sortID <- ifelse(bl_all$sampleCondition_sortID == 'OK', NA, bl_all$sampleCondition_sortID)
bl_all$samplCondition_combined <- dplyr::coalesce(bl_all$sampleCondition_field, bl_all$sampleCondition_sortID, bl_all$sampleCondition_para)
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('sampleCondition_field', 'sampleCondition_sortID', 'sampleCondition_para', 'sampleCondition_expe'))]

# Remove other extraneous columns
bl_all <- bl_all[, -which(colnames(bl_all) %in% c(
  'uid_expe', 'uid_field', 'uid_para', 'uid_sortID', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 'laboratoryName',
  'coordinateUncertainty', 'elevationUncertainty', 'nlcdClass', 'etOHChangeDate', 'processingDate_id', 'processingDate_sort',
  'publicationDate_expe', 'publicationDate_field', 'publicationDate_para', 'publicationDate_sortID',
  'recordedBy', 'recordedBy_id', 'recordedBy_sort', 'release_expe', 'release_field', 'release_para', 'release_sortID',
  'sampleCode_field', 'sampleCode_sortID', 'totalLength', 'samplingProtocolVersion', 'subsampleCode_id', 'subsampleCode_sort',
  'setDate', 'siteID', 'plotType', 'namedLocation', 'eventID', 'domainID', 'trappingDays', 'cupStatus', 'fluidLevel', 'lidStatus',
  'sampleCollected', 'remarks_para', 'remarks_sortID', 'remarks_expe', 'remarks_field', 'samplingImpractical', 'identificationRemarks',
  'samplCondition_combined'
))]

bl_all <- bl_all[which(bl_all$sampleType == 'carabid'), ]
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('sampleType'))]

# Reorder columns for csv
bl_reo <- bl_all[, c(
  'plotID', 'elevation', 'trapID', 'collectDate', 'sampleID', 'subsampleID',
  'individualID', 'individualCount', 'taxonID', 'scientificName', 'sex',
  'scientificNameAuthorship', 'taxonRank', 'morphospeciesID', 'nativeStatusCode'
)]

write.csv(bl_reo, '../results/tables/bet_counts_simple.csv', row.names = F)
