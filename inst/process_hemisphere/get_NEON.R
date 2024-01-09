# BEM 1/9/24

library(neonUtilities)

#
xx <- neonUtilities::loadByProduct(dpID = "DP1.10017.001", site = "PUUM")
Sys.time() # [1] "2024-01-09 09:50:29 HST"
xx <- xx[['dhp_perimagefile']]

write.csv(xx, './inst/process_hemisphere/dhp_perimagefile.csv', row.names = F)
