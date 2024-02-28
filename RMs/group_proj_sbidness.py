#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
SBID lists to be defined and then read in by groups_basic_pipeline_output_fdf_analyser_ver10_collab_restructure.py
or later versions.
"""

##########################################################################################
###################### SBID ACCEPTANCE / REJECTION HANDLING ##############################
##########################################################################################

import numpy as np
import os

basedir = os.path.abspath(".").replace("\\","/")

# List of SBIDs that we want to ignore -- Galactic, band 2, rejected, etc.
ignore_SBIDs_list_initial = np.load(basedir+'ignore_sbids/ignore_SBIDs_list.npy')
new_bad_SBIDs = np.load(basedir+'ignore_sbids/new_bad_SBIDs.npy')
wallaby_SBIDs = np.load(basedir+'ignore_sbids/wallaby_SBIDs.npy')
POSSUM_pilot_band2_SBIDs = np.load(basedir+'ignore_sbids/POSSUM_pilot_band2_SBIDs.npy')
ignore_SBIDs_list = list(np.concatenate((ignore_SBIDs_list_initial, new_bad_SBIDs, wallaby_SBIDs, POSSUM_pilot_band2_SBIDs)))

# Ensures we use the same pool of SBIDs for the real and control experiments, so we don't for example, add in a new field with higher / lower RM variance
SBIDs_in_non_shifted_sample = np.load(basedir+'ignore_sbids/SBIDs_in_non_shifted_sample.npy')

# ALT --- just specify main survey sbids that are at |b|>20 that are in the POSSUM Pipeline Status Jun 27 spreadsheet
include_sbids_band1_POSSUM_full = [
	9287, # EMU pilot 1 band 1
	9325, # EMU pilot 1 band 1
	9351, # EMU pilot 1 band 1
	9410, # EMU pilot 1 band 1
	9434, # EMU pilot 1 band 1
	9437, # EMU pilot 1 band 1
	9442, # EMU pilot 1 band 1
	9501, # EMU pilot 1 band 1
	10083, # EMU pilot 1 band 1
	10635, # EMU pilot 1 band 1
	31375, # pilot 2 band 1
	#33287, # pilot 2 band 1        ---> Faraday Spectra not available on CASDA
	33370, # pilot 2 band 1 I
	#33400, # pilot 2 band 1        ---> Faraday Spectra not available on CASDA
	#33459, # pilot 2 band 1        ---> Faraday Spectra not available on CASDA
	33460, # pilot 2 band 1 I
	33482, # pilot 2 band 1 I
	33509, # pilot 2 band 1 I
	33553, # pilot 2 band 1 I
	43137, # pilot 2 band 1 I
	43237, # pilot 2 band 1 I
	43773, # pilot 2 band 1 I
	45761, # >>>>> Start full survey SBIDs
	45781,
	45811,
	46925,
	46943,
	46946,
	46951,
	#46955, Very strange source dist. Surely RFI/Solar affected. See /Users/canders/Dropbox/Scripts/python/RM_interpol/Runs/plots/InferenceParams_46955/data/angular_position_data_projection.png
	#46957, polar cap obs, whose leakage corrections are BAD
	46962,
	46966,
	46971,
	46976,
	#46978, polar cap obs, whose leakage corrections are BAD
	46982,
	46984,
	46986,
	47034,
	47130,
	47136,
	49990,
	49992,
	#50008, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	50011,
	#50047, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	#50048, #Rejected by POSSUM. Reason: No beamwise validation plots available. 
	#50049, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	#50050, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	#50181, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	#50182, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	#50220, #Northern declination field, two interleaves, greater uncertainty, exclude for time being
	50230,
	50413,
	50415,
	50423, #Not yet released for POSSUM
	50538, 
	50787, #Not yet released for POSSUM
	51431, #Not yet released for POSSUM
	51574, #Not yet released for POSSUM
	51797, #Not yet released for POSSUM
	51818, #Not yet released for POSSUM
	51819] #Not yet released for POSSUM
51927 #Abell one.	
	


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import zoom

def downsample_fits_to_jpeg(fits_filename, jpeg_filename, downscale_factor=10):
    """
    Reads a FITS file and produces a downsampled JPEG.
    
    Parameters:
    - fits_filename: Path to the FITS file.
    - jpeg_filename: Path where the output JPEG will be saved.
    - downscale_factor: Factor by which the image/cube should be downscaled.
    """
    # Load the FITS file
    with fits.open(fits_filename) as hdulist:
        data = hdulist[0].data

    # Check if there are any singleton dimensions and squeeze them out
    data = np.squeeze(data)

    # If the data is a 3D cube, average along the third axis
    if len(data.shape) == 3:
        data = np.mean(data, axis=0)
    
    # Downsample the 2D image
    data_downsampled = zoom(data, 1/downscale_factor)

    # Save the downsampled image as a JPEG
    plt.imsave(jpeg_filename, data_downsampled, cmap='gray')

