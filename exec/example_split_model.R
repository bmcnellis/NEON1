

library(glmmTMB)

set.seed(1)

s0 <- glmmTMB::Salamanders
s1 <- s0
s1$count_logi <- s1$count < 0.5
s2 <- s0[which(s0$count > 0), ]

zinb <- glmmTMB(count ~ spp * mined + (1|site), zi = ~spp * mined, data = s0, family = nbinom2)

hnb <- glmmTMB(count ~ spp * mined + (1|site), zi = ~spp * mined, data = s0, family = truncated_nbinom2)

# The count < 0.5 converts the 0/1 to a TRUE/FALSE binary, with TRUE centered on 0
z0 <- glmmTMB(count < 0.5 ~ spp * mined, data = s0, family = binomial)
z1 <- glmmTMB(count_logi ~ spp * mined, data = s1, family = binomial)

pos <- glmmTMB(count ~ spp * mined + (1|site), data = s2, family = truncated_nbinom2)

logLik(pos)
#'log Lik.' -491.5107 (df=16)
logLik(z0)
#'log Lik.' -315.2394 (df=14)
logLik(z1)

logLik(pos) + logLik(z0)
logLik(pos) + logLik(z1)
#'log Lik.' -806.7501 (df=16)
logLik(hnb)
#'log Lik.' -806.7501 (df=30)
