# BEM 1/8/24

# Uses python package at: https://github.com/luke-a-brown/hemipy
import hemipy

data_dir = '/media/bem/data/NEON'
img_dir = '/media/bem/data/NEON/temp_img_dir'
proc_dir = '/media/bem/data/NEON/processed_hemisphere_photos'

# Use the following citation for this method:
# Brown, L.A., Morris, H., Leblanc, S., Bai, G., Lanconelli, C., Gobron, N.,
#     Meier, C., Dash, J. HemiPy: A Python module for automated estimation of
#     forest biophysical variables and uncertainties from digital hemispherical
#     photographs, Methods Ecol. Evol.

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

