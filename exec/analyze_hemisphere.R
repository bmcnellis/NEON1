# BEM 1/8/2024
library(NEON1)
library(neonUtilities)
library(reticulate)

data_dir <- 'C:/Users/BrandonMcNellis/Documents/NEON_data'
img_dir <- 'C:/Users/BrandonMcNellis/Documents/scratch/R_temp'
proc_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/data/processed_hemisphere_photos'

# Use the following citation for this method:
# Brown, L.A., Morris, H., Leblanc, S., Bai, G., Lanconelli, C., Gobron, N.,
#     Meier, C., Dash, J. HemiPy: A Python module for automated estimation of
#     forest biophysical variables and uncertainties from digital hemispherical
#     photographs, Methods Ecol. Evol.

# Uses python package at: https://github.com/luke-a-brown/hemipy

# Reticulate setup - only needs to be run once per script run, but should
# be deleted afterwards to reduce version control size.
ret_dir <- file.path(getwd(), 'r-python')
reticulate::virtualenv_create(ret_dir)
lib_dir <- file.path(ret_dir, 'Lib/site-packages')
system2(paste0('pip install https://github.com/luke-a-brown/hemipy/archive/refs/tags/v0.1.2.zip --upgrade --target="', lib_dir, '"'))
# If `system2` fails, you can run the pip command in git bash, e.g.:
# pip install https://github.com/luke-a-brown/hemipy/archive/refs/tags/v0.1.2.zip --upgrade --target="C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/NEON1/r-python/Lib/site-packages"

hemi0 <- neonUtilities::stackByTable(filepath = file.path(data_dir, 'NEON_hemispheric-photos-veg.zip'), savepath = 'envt')
hemi1 <- hemi0$dhp_perimagefile
# image files are in hemi1$imageFileUrl

for (i in seq(nrow(hemi1))) {
  if (i > 1) break

  f_URL <- hemi1[i, 'imageFileUrl']
  f_loc <- file.path(img_dir, hemi1[i, 'imageFileName'])
  f_out <- file.path(proc_dir, paste0(basename(hemi1[i, 'imageFileName']), 'processed'))

  # Only run if the output file is not present
  if (!file.exists(f_out)) {

    # If local file is missing, download it
    if (!file.exists(f_loc)) {
      utils::download.file(f_URL, f_loc)
    }

    # Process the local file

    # Clean-up
    file.remove(f_out)
  }

}

