# BEM 1/8/24

# Uses python package at: https://github.com/luke-a-brown/hemipy
# Deconstructed for simplicity and portability

# Imports for parent  script:
import pandas, os
from urllib.request import urlretrieve

# Imports for hemipy_mod:
import glob, datetime, math, rawpy
import numpy as np
from skimage import filters, measure
import imageio as iio
from uncertainties import unumpy, umath

# Local import for hemipy functions:
import hemipy_mod

# Input data and directories:
df = pandas.read_csv('dhp_perimagefile.csv')
img_dir = '/media/bem/data/NEON/temp_img_dir'
proc_dir = '/media/bem/data/NEON/processed_hemisphere_photos'

# Work loop:
for i in range(1, len(df.index)):
  
  # Debug break:
  if i == 2:
    break
  
  # Files:
  f_URL = df.loc[i, 'imageFileUrl']
  f_loc = os.path.join(img_dir, df.loc[i, 'imageFileName'])
  f_out = os.path.join(proc_dir, os.path.splitext(df.loc[i, 'imageFileName'])[0] + '_processed')

  # Only run if the output file is not present
  if not os.path.isfile(f_out):
    urlretrieve(f_URL, f_loc)
    
    # Calculate zenith
    # Requires: size in pixels of image, opitcal center of image, and polynomial projection
    
    # Calculate azimuth

    # OLD R CODE BELOW
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

