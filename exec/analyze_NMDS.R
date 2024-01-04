# BEM 12 May 2023
# updated 18 Dec 2023

# Setup
#install_github("admahood/neonPlantEcology")
library(neonPlantEcology)
#library(vegan)
library(NEON1)
library(gllvm)

set.seed(1)

# Filenames

# (1) How much of plant composition is explained by time or disturbance at Puʻu Makaʻala?

# Create matrix of cover data
data_div <- neonPlantEcology::npe_download(sites = c("PUUM"))
# can add WREF and GUAN if need be
data_df <- neonPlantEcology::npe_longform(data_div, scale = 'plot')
data_df$cover <- data_df$cover / 100
data_df <- data_df[-which(data_df$taxonID == '2PLANT'), ]
colnames(data_df)[which(colnames(data_df) == 'eventID')] <- 'year'
data_mat <- neonPlantEcology::npe_community_matrix(data_div)
data_mat <- data_mat[, -which(colnames(data_mat) == '2PLANT')]
#data_mat <- data_mat[grepl('2018', row.names(data_mat)), ]
#data_mat <- data_mat[rowSums(data_mat) > 0, ]
data_mat <- data_mat / 100

# PUUM
PUUM_mat <- data_mat[grepl('PUUM', row.names(data_mat)), ]
PUUM_mat <- PUUM_mat[, colSums(PUUM_mat) > 0]
#NEON1::NMDS_screeplot(PUUM_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions, can try to get away with 3
# this is UNCONSTRAINED ORDINATION
m_PUUM <- metaMDS(PUUM_mat, k = 3, try = 50, trymax = 200)
m_PUUM$stress # 0.1367
# for plotting/stats:
PUUM_g_y <- as.factor(sapply(strsplit(row.names(PUUM_mat), '_'), \(xx) xx[4]))
PUUM_g_p <- as.factor(sapply(strsplit(row.names(PUUM_mat), '_'), \(xx) xx[2]))
# significance:
PERMANOVA_PUUM <- vegan::adonis2(PUUM_mat ~ PUUM_g_p + PUUM_g_y, by = 'margin')
# plot and year are both significant
pairwise_PUUM_y <- pairwiseAdonis::pairwise.adonis(PUUM_mat, factors = PUUM_g_y)
pairwise_PUUM_y <- pairwise_PUUM_y[which(pairwise_PUUM_y$p.adjusted <= 0.05), ]
# no years are significantly different from one another
pairwise_PUUM_p <- pairwiseAdonis::pairwise.adonis(PUUM_mat, factors = PUUM_g_p)
pairwise_PUUM_p <- pairwise_PUUM_p[which(pairwise_PUUM_p$p.adjusted <= 0.05), ]
# no plots are significantly different from one another

# GUAN
GUAN_mat <- data_mat[grepl('GUAN', row.names(data_mat)), ]
GUAN_mat <- GUAN_mat[, colSums(GUAN_mat) > 0]
NEON1::NMDS_screeplot(GUAN_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions
m_GUAN <- metaMDS(GUAN_mat, k = 4, try = 50, trymax = 200)
m_GUAN$stress # 0.1647
# for plotting/stats:
GUAN_g_y <- as.factor(sapply(strsplit(row.names(GUAN_mat), '_'), \(xx) xx[4]))
GUAN_g_p <- as.factor(sapply(strsplit(row.names(GUAN_mat), '_'), \(xx) xx[2]))
# significance:
PERMANOVA_GUAN <- vegan::adonis2(GUAN_mat ~ GUAN_g_p + GUAN_g_y, by = 'margin')
# plot and year are both significant
pairwise_GUAN_y <- pairwiseAdonis::pairwise.adonis(GUAN_mat, factors = GUAN_g_y)
pairwise_GUAN_y <- pairwise_GUAN_y[which(pairwise_GUAN_y$p.adjusted <= 0.05), ]
# 28 year-pairs have significant differences, mostly 2015
pairwise_GUAN_p <- pairwiseAdonis::pairwise.adonis(GUAN_mat, factors = GUAN_g_p)
pairwise_GUAN_p <- pairwise_GUAN_p[which(pairwise_GUAN_p$p.adjusted <= 0.05), ]
# no plots are significantly different from one another

# WREF
WREF_mat <- data_mat[grepl('WREF', row.names(data_mat)), ]
WREF_mat <- WREF_mat[, colSums(WREF_mat) > 0]
NEON1::NMDS_screeplot(WREF_mat, k = 15)
# looks like stress drops are smaller after 4 dimensions, can try 3
m_WREF <- metaMDS(WREF_mat, k = 3, try = 50, trymax = 200)
m_WREF$stress # 0.1902
# for plotting/stats:
WREF_g_y <- as.factor(sapply(strsplit(row.names(WREF_mat), '_'), \(xx) xx[4]))
WREF_g_p <- as.factor(sapply(strsplit(row.names(WREF_mat), '_'), \(xx) xx[2]))
# significance:
PERMANOVA_WREF <- vegan::adonis2(WREF_mat ~ WREF_g_p + WREF_g_y, by = 'margin')
# plot and year are both significant
pairwise_WREF_y <- pairwiseAdonis::pairwise.adonis(WREF_mat, factors = WREF_g_y)
pairwise_WREF_y <- pairwise_WREF_y[which(pairwise_WREF_y$p.adjusted <= 0.05), ]
# 28 year-pairs have significant differences, mostly 2015
pairwise_WREF_p <- pairwiseAdonis::pairwise.adonis(WREF_mat, factors = WREF_g_p)
pairwise_WREF_p <- pairwise_WREF_p[which(pairwise_WREF_p$p.adjusted <= 0.05), ]
# no plots are significantly different from one another

# If both space and time are significantly related to species assemblages in all
# 3 plots, can a multilevel model predict which aspects of space and time are
# responsible for this? how many other variables are there?

# can we use an LVM with 2 latent variables? then the 2 latent variables can be
# assessed with a biplot, like a PCA/NMDS

# can we use the residual correlations between species to speak of associations?
# is residual species correlation a nuisance variable or variable of interest? this
# is a *model based approach to unconstrained ordination*, see:

# Walker, S.C. and Jackson, D.A. (2011) Random-effects ordina-
#  tion: describing and predicting multivariate correlations and co-
#  occurrences. Ecol. Monogr. 81, 635–663
# Hui, F.K.C. et al. (2015) Model-based approaches to uncon-
#  strained ordination. Methods Ecol. Evol. 6, 399–411

# can we exploit temporal and spatial autocorrelation for the purposes of seeing
# if there is temporal or spatial trends across the data? i.e. no temporal autocorrelation
# means that there is no inertia within the system, and it will respond rapidly
# to changes in precipitation or temperature. if there is no spatial autocorrelation,
# it means the plots within the site are operating independently of one another

# are site effects relevant? do we combine site and species effects into aggregate
# abundance, or do we seperate them into site abundance and ordinated composition?
# the latter method just needs an extra intercept (fixed or random) fit to site

# if we use an LVM, we assume that there is a 'latent' gradient that describes
# the sites. sites are normal along the gradient, and species are linearly related
# to gradient position


