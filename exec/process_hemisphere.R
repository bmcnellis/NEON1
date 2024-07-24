# BEM July 2024

library(NEON1)
library(dplyr)

sph <- function(vec) {

  v <- sapply(strsplit(vec, '\\+/-'), \(xx) xx[1])
  v <- ifelse(v == 'NULL', NA, v)
  v <- as.numeric(v)

  e <- sapply(strsplit(vec, '\\+/-'), \(xx) xx[2])
  e <- as.numeric(e)

  data.frame(val = v, err = e)

}

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
    pai = sph(pai_up)[, 1],
    pai_err = sph(pai_up)[, 2],
    fcover = sph(fcover_up)[, 1],
    fcover_err = sph(fcover_up)[, 2],
    fipar = sph(fipar_up)[, 1],
    fipar_err = sph(fipar_up)[, 2],
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

usethis::use_data(dhp, overwrite = T)
