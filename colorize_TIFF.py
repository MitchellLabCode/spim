"""
Simple script to colorize a TIFF file based on custom colors for endoderm/circular muscle/longitudinal muscle
"""

import os
import h5py
import numpy as np
import tifffile as tiff
from pathlib import Path
from scipy.ndimage import zoom as ndi_zoom

# ==== USER CONFIG ============================================================
datdir = Path("/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/figures/byn_domain")  # Base directory for data
H5_PATH = os.path.join(datdir, "APDV_Time_000005_8x_Probabilities.h5")        # Path to HDF5 file with probabilities
H5_DATASET = "/exported_data"                # Dataset path within H5 file
PROB_CHANNEL = 0                             # Which channel/class to use (0 for first)
DOWNSAMPLE_FACTOR = 4                        # Downsampling factor of prob data (e.g., 4 = 4x downsampled)
THRESHOLD = 0.5                              # Probability threshold for mask (0.0-1.0)
TIFF_INPUT = os.path.join(datdir, "3DProjections_of_ResliceofAPDV_Time_000005_2x_masked.tif")             # Input TIFF file
TIFF_OUTPUT = os.path.join(datdir, "colored_3DProjections_of_ResliceofAPDV_Time_000005_2x_masked.tif")    # Colorized TIFF file
TIFF_CHANNEL_IDX = 0                         # Which channel to extract if 4D (0 for first)
TIFF_CHANNEL_AXIS = 1                       # Which axis holds channels if 4D (-1=last, 0=first, 1=second, etc)
PRESERVE_DTYPE = True                        # Cast output back to input dtype
# ============================================================================

