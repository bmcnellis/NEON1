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
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('class', 'family', 'genus', 'kingdom', 'order', 'phylum'))]
# Remove extraneous sample identification
# This assumes that the expert identification supersedes the parataxonomist identification
bl_all$scientificName_coal <- dplyr::coalesce(bl_all$scientificName_expe, bl_all$scientificName_para)
bl_all$scientificName_coal <- ifelse(is.na(bl_all$scientificName_expe), bl_all$scientificName_coal, bl_all$scientificName_expe)
# What should supercede - the scientific name from the sorting, or the para/expe ID's?
bl_all <- bl_all[, -which(colnames(bl_all) %in% c('scientificName_coal', 'scientificName_para', 'scientificName_expe'))]
