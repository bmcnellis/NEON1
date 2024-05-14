# BEM May 2024

library(neonUtilities)

# get 'Plant presence and percent cover' for PUUM, DP1.10058.001

PUUM_DP110058001 <- neonUtilities::loadByProduct('DP1.10058.001', 'PUUM')

d0 <- PUUM_DP110058001[['div_1m2Data']]
d0 <- d0[, c('domainID', 'siteID', 'decimalLatitude', 'decimalLongitude', 'elevation', 'plotType', 'nlcdClass', 'plotID', 'subplotID', 'endDate', 'release')]
d0 <- d0[!duplicated(d0), ]

site_meta_PUUM_DP110058001 <- d0
usethis::use_data(site_meta_PUUM_DP110058001)
