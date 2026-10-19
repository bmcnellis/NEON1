# BEM October 2026

# Purpose: download the raw NEON data products and cache them on disk, so that
#          import_data.R can be re-run without touching the NEON API.

# Why this is a separate script:
#     As of June 2026 NEON requires an API token for every data download, so
#     loadByProduct() now fails outright without one. The token is a private
#     credential: it lives in exec/NEON_TOKEN.R, which is listed in .gitignore
#     and must never be committed. Keeping the downloads here means the token is
#     read by exactly one script, and that script is the only part of the
#     pipeline that needs network access. The downloads are also slow and return
#     provisional data that changes between NEON releases, so import_data.R now
#     works from a fixed local copy instead of re-downloading every time. Re-run
#     this script only when the data should be refreshed, and update the
#     "last run" date below when you do.

# Last run:  manually by BEM 19/10/26

# NOTE: this script is run from NEON1/exec/, not from the package root like the
#       rest of the pipeline, because it writes outside the package.

library(neonUtilities) # 4.0.2

### NEON API token
token_fl <- 'NEON_TOKEN.R'
if (!file.exists(token_fl)) {
  stop('NEON API token not found at ', token_fl,
       ' - set the working directory to NEON1/exec')
}
source(token_fl)
stopifnot(exists('NEON_TOKEN'), nzchar(NEON_TOKEN))

### Output directory
neon_dir <- '../../data_NEON'
if (!dir.exists(neon_dir)) dir.create(neon_dir, recursive = T)

### Download data products
# d0 DP1.10058.001 plant presence and percent cover
# d1 DP1.10098.001 vegetation structure
# d2 DP1.10047.001 soil chemical properties, distributed initial characterization
# d3 DP1.10086.001 soil physical and chemical properties, periodic
# d4 DP1.10026.001 plant foliar traits
d0 <- neonUtilities::loadByProduct('DP1.10058.001', 'PUUM', include.provisional = T, check.size = F, token = NEON_TOKEN)
d1 <- neonUtilities::loadByProduct('DP1.10098.001', 'PUUM', include.provisional = T, check.size = F, token = NEON_TOKEN)
d2 <- neonUtilities::loadByProduct('DP1.10047.001', 'PUUM', include.provisional = T, check.size = F, token = NEON_TOKEN)
d3 <- neonUtilities::loadByProduct('DP1.10086.001', 'PUUM', include.provisional = T, check.size = F, token = NEON_TOKEN)
d4 <- neonUtilities::loadByProduct('DP1.10026.001', 'PUUM', include.provisional = T, check.size = F, token = NEON_TOKEN)

### Save, one file per product, named for the object import_data.R expects
for (ii in paste0('d', 0:4)) {
  save(list = ii, file = file.path(neon_dir, paste0(ii, '.rda')))
}

# record what was downloaded and when, since the products are provisional
writeLines(
  c(
    paste('downloaded', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')),
    paste('neonUtilities', as.character(utils::packageVersion('neonUtilities'))),
    paste0('d', 0:4, ': ', c('DP1.10058.001', 'DP1.10098.001', 'DP1.10047.001', 'DP1.10086.001', 'DP1.10026.001'),
           ' - ', sapply(list(d0, d1, d2, d3, d4), \(xx) paste(length(xx), 'tables')))
  ),
  file.path(neon_dir, 'download_log.txt')
)

