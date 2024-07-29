# R 4.2.2.
library(neonUtilities) # 2.4.1
library(dplyr) # 1.1.4
library(tidyr) # 1.3.1
library(glmmTMB) # 1.1.8
library(DHARMa) # 0.4.6

# Download and format data
d0 <- neonUtilities::loadByProduct('DP1.10058.001', 'WREF', check.size = F) |>
  (`[[`)('div_1m2Data') |>
  dplyr::select(c(plotID, endDate, scientificName, percentCover)) |>
  dplyr::mutate(endDate = as.character(endDate), cover = percentCover / 100, percentCover = NULL) |>
  dplyr::filter(endDate > as.POSIXct(as.Date('2021-01-01'))) |>
  na.omit() |>
  tidyr::pivot_wider(names_from = scientificName, values_from = cover, values_fn = mean) |>
  tidyr::pivot_longer(cols = -c(1, 2), names_to = 'scientificName', values_to = 'cover') |>
  dplyr::mutate(cover = ifelse(is.na(cover), 0, cover))

# glmmTMB model with only random effects for date/plot fits well, but no information on species
c0 <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1e3, eval.max = 1e3), start_method = list(method = NULL, jitter.sd = 0))
m0 <- glmmTMB::glmmTMB(cover ~ (1|endDate/plotID), data = d0, family = ordbeta(), control = c0, start = list())
s0 <- DHARMa::simulateResiduals(m0)

DHARMa::plotResiduals(s0, form = d0$endDate)
DHARMa::plotResiduals(s0, form = d0$plotID)
DHARMa::plotResiduals(s0, form = d0$scientificName)

# glmmTMB multi-species model with RR grouped by plotID fit poorly
m1 <- glmmTMB(cover ~ 1 + rr(scientificName + 0|plotID, d = 2), data = d0, family = ordbeta(), control = c0, start = list())
s1 <- DHARMa::simulateResiduals(m1)
DHARMa::plotResiduals(s1, form = d0$endDate)
DHARMa::plotResiduals(s1, form = d0$plotID)
DHARMa::plotResiduals(s1, form = d0$scientificName)
