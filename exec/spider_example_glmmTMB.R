# from:
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html

library(mvabund)
library(glmmTMB)

data(spider)
## organize data into long format
sppTot <- sort(colSums(spider$abund), decreasing = TRUE)
tmp <- cbind(spider$abund, spider$x)
tmp$id <- 1:nrow(tmp)
spiderDat <- reshape(tmp,
                     idvar = "id",
                     timevar = "Species",
                     times =  colnames(spider$abund),
                     varying = list(colnames(spider$abund)),
                     v.names = "abund",
                     direction = "long")
## fit rank-reduced models with varying dimension
fit_list <- lapply(2:10,
                   function(d) {
                     fit.rr <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = d),
                                       data = spiderDat)
                   })
## compare fits via AIC
aic_vec <- sapply(fit_list, AIC)
aic_vec - min(aic_vec, na.rm = TRUE)

# just fit the one model
fit0 <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 2),
                data = spiderDat)
