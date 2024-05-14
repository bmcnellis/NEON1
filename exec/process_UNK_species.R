# BEM May 2024

# TODO: 'Rhynchospora chinensis' assumed native, b/c native ssp, unsure on other ssp presence

library(NEON1)

key_nativeStatusCode <- data.frame(taxonID = character(), nativeStatusCode = character())

data_div <- NEON1::PUUM_div
data_df <- neonPlantEcology::npe_longform(data_div, scale = 'plot')

unks <- unique(data_df$scientificName[which(data_df$nativeStatusCode == 'UNK')])
unkc <- unique(data_df$taxonID[which(data_df$scientificName %in% unks)])
#unk_code <- unique(data_df$taxonID[which(data_df$nativeStatusCode == 'UNK')])

res <- vector('list', 3)

for (i in seq_along(unks)) {

  ii <- unks[i]
  i0 <- sapply(strsplit(ii, ' '), \(xx) xx[1])

  idf <- data_df[grepl(i0, data_df$scientificName), ]

  res[[1]][[i]] <- as.data.frame(idf)
  res[[1]][[i]] <- res[[1]][[i]][order(res[[1]][[i]]$cover, decreasing = T), ]

  res[[2]][[i]] <- res[[1]][[i]][, c('taxonID', 'nativeStatusCode', 'scientificName')]
  res[[2]][[i]] <- res[[2]][[i]][!duplicated(res[[2]][[i]]), ]

  res[[3]][[i]] <- table(idf$nativeStatusCode)

}

names(res[[1]]) <- unks
names(res[[2]]) <- unks
names(res[[3]]) <- unks

# Asplenium
#   Asplenium has invasive and native members, hard to decide between them

# Carex
#   All of the Carex on the plot are native, assume the UNK are also native
key_nativeStatusCode <- rbind(key_nativeStatusCode, data.frame(taxonID = c('CAREX', 'CAREXSPP'), nativeStatusCode = c('N', 'N')))

# Fragaria sp/Fragaria spp
#   USDA Plants:
#     N: Fragaria chiloensis (L.) Mill. ssp. sandwicensis (Decne.) Staudt
#       "Less common on Hawaii than it was in past decades" per as of Wagner et al. (1999)
#       Wagner et al. (1999) says it occurs up to 3000ft, but others have its habit as near-ocean/beach
#       Only record on Hawaii island at NTBG is from 1855
#     N: Fragaria chiloensis (L.) Mill.
#     I: Fragaria vesca L.
#       NTBG has specimens from PUUM N.A.R. in 1988, and 1996 from Kohala
#     I: Fragaria vesca L. ssp. americana (Porter) Staudt
#
# Based on above, FRAGA/FRAGASPP is likely Fragaria vesca, and should be marked I

key_nativeStatusCode <- rbind(key_nativeStatusCode, data.frame(taxonID = c('FRAGA', 'FRAGASPP'), nativeStatusCode = c('I', 'I')))
