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

key_cols_1 <- c('namedLocation', 'domainID', 'siteID', 'plotID', 'setDate', 'collectDate', 'individualID')

# Combine the multiple data streams into one sparse dataframe
# bet_parataxonomistID and bet_expertTaxonomistID need to be combined first
# All `individualID` in bet_expertTaxonomistID are in bet_parataxonomistID,
bl_id <- dplyr::left_join(bl_p, bl_x, by = key_cols_1, suffix = c('_para', '_expe'))

# All `subsampleCode` in bet_parataxonomistID are in bet_sorting, reverse is not true
bl_sd <- dplyr::left_join(bl_s, bl_d, by = c(key_cols, 'subsampleID'), suffix = c('_sort', '_para'))

bl_fs <- dplyr::left_join(bl_fs, bl_d, by = c('plotID', 'trapID', 'collectDate'), suffix = c('_fs', '_para'))
# Field data is last - the unique identifiers are plotID, trapID, collectDate

# Combine the multiple data streams into one sparse dataframe
# Parent dataframe is bet_parataxonomistID
# Modify expert data to be merged with para data
# These columns share names with para data, but have different values:
bl_rnx <- c(
  uid_ex = 'uid',
  identifiedDate_ex = 'identifiedDate',
  taxonID_ex = 'taxonID',
  taxonRank_Ex = 'taxonRank',
  scientificName_ex = 'scientificName',
  scientificNameAuthorship_ex = 'scientificNameAuthorship',
  identificationQualifier_ex = 'identificationQualifier',
  morphospeciesID_ex = 'morphospeciesID',
  identificationReferences_ex = 'identificationReferences',
  identifiedBy_ex = 'identifiedBy',
  nativeStatusCode_ex = 'nativeStatusCode',
  remarks_ex = 'remarks',
  identificationHistoryID_ex = 'identificationHistoryID'
)
bl_x <- dplyr::rename(bl_x, any_of(bl_rnx))
# The rest are key columns:
#     namedLocation, domainID, siteID, plotID, setDate, collectDate,
#     individualID, sampleCondition, publicationDate
bl_c <- dplyr::left_join(bl_d, bl_x)
# Modify field data to be merged with combined ID data
bl_rnf <- c(
  uid_f = 'uid',
  recordedBy_f = 'recordedBy',
  remarks_f = 'remarks'
)
bl_f <- dplyr::rename(bl_f, any_of(bl_rnf))
# The rest are key columns: namedLocation, domainID, siteID, plotID, trapID,
#     setDate, collectDate, sampleCondition, publicationDate, release
bl_0 <- dplyr::left_join(bl_c, bl_f)

bl_rns <- c(

)
# The rest are key columns:
#     namedLocation, domainID, siteID,
