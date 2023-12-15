# BEM May 2023
library(NEON1)
library(neonUtilities)
library(dplyr)
import::from(magrittr, "%>%")
library(ggplot2)

ppt_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_precipitation.zip'
temp_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_temp-air-single.zip'
wind_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_wind-2d.zip'
par_zip <- 'C:/Users/BrandonMcNellis/Documents/NEON_data/NEON_par.zip'

clim_pfile_1 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/clim_1.tif'
wpar_pfile_1 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/wpar_1.tif'
ppt_pfile_2 <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/results/figures/ppt_2.tif'

# Precipitation DP1.00006.001
ppt_zip <- neonUtilities::zipsByProduct()
ppt_tbl <- neonUtilities::stackByTable(ppt_zip, savepath = 'envt')
# PRIPRE_30min: Primary Precipitation pooled over 30 minutes
# No ABBY in primary precipitation?
# SECPRE_30min: Secondary Precipitation pooled over 30 minutes
ppt <- ppt_tbl[['SECPRE_30min']] %>%
  # only use if quality flag == 0
  filter(secPrecipRangeQF == 0) %>%
  # drop unnecessary columns
  select(-c(domainID, horizontalPosition, verticalPosition, secPrecipRangeQF, secPrecipSciRvwQF, publicationDate, release)) %>%
  # add variable for date
  mutate(date = as.Date(startDateTime)) %>%
  filter(date >= min(date[which(siteID == 'PUUM')])) %>%
  # add variable for month/year ('m_y')
  #mutate(m_y = paste(format(date, '%m'), format(date, '%Y'), sep = '_')) %>%
  mutate(m_y = lubridate::floor_date(date, 'month')) %>%
  # summarize by site and date
  #group_by(siteID, date) %>%
  group_by(siteID, m_y) %>%
  summarize(ppt = sum(secPrecipBulk)) %>%
  # scale ppt to cm
  #mutate(ppt = ppt / 10) %>%
  ungroup()

# Single aspirated air temperature (DP1.00002.001)

temp_tbl <- neonUtilities::stackByTable(temp_zip, savepath = 'envt')

temp <- temp_tbl[['SAAT_30min']] %>%
  filter(finalQF == 0) %>%
  select(-c(domainID, horizontalPosition, verticalPosition, endDateTime, tempSingleNumPts, tempSingleExpUncert, tempSingleStdErMean, publicationDate, release)) %>%
  mutate(date = as.Date(startDateTime)) %>%
  filter(date >= min(date[which(siteID == 'PUUM')])) %>%
  mutate(m_y = lubridate::floor_date(date, 'month')) %>%
  group_by(siteID, m_y) %>%
  summarize(tmean = mean(tempSingleMean), tmin = mean(tempSingleMinimum), tmax = mean(tempSingleMaximum)) %>%
  ungroup()

# Join ppt/temp by siteID/m_y
clim <- dplyr::left_join(ppt, temp, by = c('siteID', 'm_y'))

# Summary statistics for climate variables
# min of ppt in PUUM
clim %>%
  filter(siteID == 'PUUM') %>%
  slice_min(ppt)
# 10.9 mm for June 2021
# max of ppt in PUUM
clim %>%
  filter(siteID == 'PUUM') %>%
  slice_max(ppt)
# 820 ppt for February 2023

# Relative humidity DP1.00098.001
# Shortwave and longwave radiation (net radiometer) DP1.00023.001

# correlation between sites in PPT
ppt_xx <- tidyr::pivot_wider(ppt, names_from = siteID, values_from = ppt)
cor(ppt_xx$PUUM, ppt_xx$ABBY, method = 'pearson', use = 'pairwise.complete.obs')
cor(ppt_xx$PUUM, ppt_xx$GUAN, method = 'pearson', use = 'pairwise.complete.obs')

# Climatograph

#t_offset <- 0
yp <- c(0, 850)   # in this example, precipitation
ys <- c(0, 31)    # in this example, temperature
b <- diff(yp) / diff(ys)
a <- yp[1] - (b * ys[1])

ppt_col <- 'grey35'

clim_plot <- ggplot(data = clim, aes(x = m_y, y = ppt)) +

  geom_col(fill = ppt_col) +
  #geom_line(mapping = aes(y = tmean + t_offset), color = "black", size = 0.8) +
  #geom_line(mapping = aes(y = a + tmax * b), color = "red", size = 0.5) +
  #geom_line(mapping = aes(y = a + tmin * b), color = "blue", size = 0.5) +
  geom_line(mapping = aes(y = a + tmean * b), color = "black", size = 0.8) +

  scale_y_continuous(
    "Precipitation (mm)",
    #sec.axis = sec_axis(~ . -t_offset, name = "Temperature (°C)"
    sec.axis = sec_axis(~ (. - a) / b, name = "Temperature (°C)"
    )) +


  facet_wrap(~ factor(siteID, levels = c('ABBY', 'PUUM', 'GUAN')), ncol = 1) +
  #facet_wrap(~ siteID, ncol = 1, scales = 'free_y') +

  theme_bw() +
  theme(
    axis.text.y.left = element_text(color = ppt_col),
    axis.ticks.y.left = element_line(color = ppt_col),
    axis.title.y.left = element_text(color = ppt_col),
    axis.text.y.right = element_text(color = 'black'),
    axis.ticks.y.right = element_line(color = 'black'),
    axis.text.x = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey95'),
    strip.text = element_text(color = 'black')
  ) +
  labs(x = '')

HighResTiff(clim_plot, clim_pfile_1, 6, 6, 300)
# dates are from "2019-06-01" to "2023-05-01"

# Plot of total precipitation for all 3 plots for 2020, 2021, and 2022

ppt_20_22 <- ppt %>%
  filter(lubridate::year(m_y) %in% c(2020:2022)) %>%
  group_by(siteID, year = lubridate::year(m_y)) %>%
  summarize(MSP = sum(ppt)) %>%
  ungroup() %>%
  mutate(MSP = round(MSP, 0)) %>%
  mutate(year = as.character(year))

ppt_20_22_plot <- ggplot(data = ppt_20_22, aes(x = year, y = MSP, shape = siteID)) +
  geom_point() +

  scale_shape_manual(values = c(19, 4, 8)) +
  scale_x_discrete(breaks = c(2020, 2021, 2022)) +
  scale_y_continuous(breaks = seq(0, 2500, 500), limits = c(0, 2500)) +

  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black')
  ) +
  labs(x = '', y = 'Sum Annual Precipitation (mm)', shape = 'Site')

HighResTiff(ppt_20_22_plot, ppt_pfile_2, 6, 4, 300)

