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

# PPT

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

# TEMP

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

# WIND

wind_tbl <- neonUtilities::stackByTable(wind_zip, savepath = 'envt')

wind <- wind_tbl[['twoDWSD_30min']] %>%
  filter(windSpeedFinalQF == 0) %>%
  select(c(siteID, startDateTime, windSpeedMean, windSpeedMinimum, windSpeedMaximum)) %>%
  mutate(date = as.Date(startDateTime)) %>%
  filter(date >= min(date[which(siteID == 'PUUM')])) %>%
  mutate(m_y = lubridate::floor_date(date, 'month')) %>%
  group_by(siteID, m_y) %>%
  # lets do "30-minute maximum wind speed, average per day"
  summarize(wmax = mean(windSpeedMaximum, na.rm = T)) %>%
  ungroup()

# PAR

# unzip the sensor_position files
unzip(par_zip, list = T) %>%
  filter(grepl('sensor_positions', Name)) %>%
  {unzip(par_zip, files = .$Name, exdir = tempdir())}

# create sensor position key
pos_key <- file.path(tempdir(), 'NEON_par') %>%
  list.files(recursive = T, full.names = T) %>%
  lapply(read.csv) %>%
  dplyr::bind_rows() %>%
  select(c(HOR.VER, referenceDescription, zOffset)) %>%
  distinct() %>%
  na.omit() %>%
  mutate(siteID = gsub('Puu Makaala Tower', 'PUUM', referenceDescription)) %>%
  mutate(siteID = gsub('Guanica Forest Tower', 'GUAN', siteID)) %>%
  mutate(siteID = gsub('Abby Road Tower', 'ABBY', siteID)) %>%
  mutate(verticalPosition = formatC(HOR.VER * 1E3, width = 3, format = 'd', flag = '0')) %>%
  select(c(siteID, verticalPosition, zOffset))

par_tbl <- neonUtilities::stackByTable(par_zip, savepath = 'envt')

par <- par_tbl[['PARPAR_30min']] %>%
  filter(PARFinalQF == 0) %>%
  select(c(siteID, startDateTime, verticalPosition, PARMean, PARMinimum, PARMaximum)) %>%
  mutate(date = as.Date(startDateTime)) %>%
  filter(date >= min(date[which(siteID == 'PUUM')])) %>%
  mutate(m_y = lubridate::floor_date(date, 'month')) %>%
  dplyr::left_join(pos_key) %>%
  select(-verticalPosition) %>%
  group_by(siteID) %>%
  filter(zOffset %in% c(min(zOffset), max(zOffset))) %>%
  # bot is 0.21, 0.83, 0.93 for GUAN, ABBY, PUUM
  # top is 18.55, 22.96, 32.42 for ABBY, GUAN, PUUM
  mutate(sensor_pos = ifelse(zOffset == min(zOffset), 'bot', 'top')) %>%
  ungroup() %>%
  # drop times between 11PM and 4AM
  filter(lubridate::hour(startDateTime) %in% c(4:22)) %>%
  # aggregate by monthly value
  group_by(siteID, m_y, sensor_pos) %>%
  summarize(
    PARmean = mean(PARMean, na.rm = T),
    PARmin = mean(PARMinimum, na.rm = T),
    PARmax = mean(PARMaximum, na.rm = T)
  ) %>%
  ungroup()

# Join ppt/temp by siteID/m_y
clim <- dplyr::left_join(ppt, temp, by = c('siteID', 'm_y'))

# Join wind/par by siteID/m_y
wpar <- dplyr::left_join(par, wind, by = c('siteID', 'm_y')) %>%
  mutate(wmax = ifelse(sensor_pos == 'bot', wmax, NA)) %>%
  mutate(wind_lab = rep('Wind Speed', nrow(.)))

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

# All days at GUAN with > 50 mm rain
# 2019-09-25, Tropical Storm Karen
# 2020-07-30, Hurricane Isaias
# 2020-11-11/12, pre-Hurricane Iota
# 2022-09-18/19, Hurricane Fiona
# 2022-10-26/27, pre-Hurricane Lisa
# 2022-11-05, pre-Hurricane Nicole
#
# All GUAN daily totals > 50mm rain were associated with tropical storms or weather
# systems that eventually became tropical storms.

yp <- c(0, 850)   # in this example, PAR
ys <- c(0, 5)    # in this example, wind speed
b <- diff(yp) / diff(ys)
a <- yp[1] - (b * ys[1])

# Plot of PAR/wind
wpar_plot <- ggplot(data = wpar, aes(x = m_y, y = PARmean, color = sensor_pos, shape = wind_lab)) +

  geom_line(size = 0.8) +
  geom_point(mapping = aes(y = a + wmax * b), color = 'black', fill = 'black') +

  facet_wrap(~ factor(siteID, levels = c('ABBY', 'PUUM', 'GUAN')), ncol = 1) +

  scale_y_continuous(
    bquote('PAR (μmol'~~s^-1~m^-2~')'),
    sec.axis = sec_axis(~ (. - a) / b, name = bquote('Wind Speed (m'~~s^-1~')'))
  ) +

  scale_color_manual(
    values = c('grey30', 'grey70'),
    labels = c('Lower Tower', 'Upper Tower'),
    guide = guide_legend(override.aes = list(
      linetype = c('solid', 'solid'),
      shape = c(NA, NA)
    ))
  ) +

  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black'),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey95'),
    strip.text = element_text(color = 'black')
  ) +
  labs(
    x = '',
    #y = bquote('PAR (μmol'~s^-1~m^-2~')'),
    color = 'PAR',
    shape = ''
  )
# maybe add hurricanes later???
HighResTiff(wpar_plot, wpar_pfile_1, 6, 6, 300)
