# BEM 1/8/2024
library(NEON1)
library(neonUtilities)
library(reticulate)
library(exifr)


# this script is shit and reticulate is always a waste of time

#data_dir <- 'C:/Users/BrandonMcNellis/Documents/NEON_data'
data_dir <- '/media/bem/data/NEON'
#img_dir <- 'C:/Users/BrandonMcNellis/Documents/scratch/R_temp'
img_dir <- '/media/bem/data/NEON/temp_img_dir'
#proc_dir <- 'C:/Users/BrandonMcNellis/OneDrive - USDA/NEON1/data/processed_hemisphere_photos'
proc_dir <- '/media/bem/data/NEON/processed_hemisphere_photos'

# Uses python package at: https://github.com/luke-a-brown/hemipy

# Reticulate setup - only needs to be run once per script run, but should
# be deleted afterwards to reduce version control size.
ret_dir <- file.path(getwd(), 'r-python')
#reticulate::virutalenv_create(ret_dir) # for Windows
reticulate::virtualenv_create(ret_dir, python = install_python()) # for Ubuntu
#lib_dir <- file.path(ret_dir, 'Lib/site-packages') # for Windows
lib_dir <- file.path(ret_dir, 'lib/python3.9/site_packages') # for Ubuntu
system2('./r-python/bin/pip install https://github.com/luke-a-brown/hemipy/archive/refs/tags/v0.1.2.zip --upgrade')

# Python env setup - used each time the R script is used
reticulate::use_virtualenv('./r-python')
import('hemipy')

hemi0 <- neonUtilities::stackByTable(filepath = file.path(data_dir, 'NEON_hemispheric-photos-veg.zip'), savepath = 'envt')
hemi1 <- hemi0$dhp_perimagefile
# image files are in hemi1$imageFileUrl

for (i in seq(nrow(hemi1))) {
  if (i > 1) break

  f_URL <- hemi1[i, 'imageFileUrl']
  f_loc <- file.path(img_dir, hemi1[i, 'imageFileName'])
  f_out <- file.path(proc_dir, paste0(tools::file_path_sans_ext(hemi1[i, 'imageFileName']), '_processed'))

  # Only run if the output file is not present
  if (!file.exists(f_out)) {

    # If local file is missing, download it
    if (!file.exists(f_loc)) {
      utils::download.file(f_URL, f_loc)
    }

    # Process the local file
    ex0 <- exifr::read_exif(f_loc)
    h <- ex0$ImageHeight
    w <- ex0$ImageWidth
    s <- prod(h, w)

    # just need to write the process script in python here
    #reticulate::py_run_string('hemipy.zenith()')

    # Clean-up
    if (file.exists(f_out)) {
      file.remove(f_loc)
    }
  }

}

