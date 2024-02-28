#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
Functions to be read in by groups_basic_pipeline_output_fdf_analyser_ver10_collab_restructure.py
or later versions.
"""

import os
import glob
import warnings
from copy import deepcopy
from subprocess import call

import numpy as np
import pandas as pd
import pylab as pl

from astropy.io import fits, ascii
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import WCS
from astropy.table import Table, unique
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.stats import mad_std, sigma_clip

from astropy.cosmology import Planck18 as assumed_comology
from astroquery.simbad import Simbad
from astropy.utils.exceptions import AstropyWarning

from RMutils.util_RM import measure_FDF_parms

from matplotlib import pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.pyplot import cm
import matplotlib.patches as mpatches
from adjustText import adjust_text

from scipy.constants import c as c
from scipy.stats import median_abs_deviation as mad, iqr, stats, bootstrap, binned_statistic_2d
from statsmodels.robust.scale import huber

import angles as ang

import random

import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

#import healpy as hp
import h5py
import dask
from dask.distributed import Client


##########################################################################################
###################################### FUNCITONS #########################################
##########################################################################################

def assert_all_sublists_same_size(master_list):
	"""
	Check if all lists in a master list have the same size.

	Parameters:
	master_list (list): A list of lists to be checked for equal size.

	Returns:
	bool: True if all lists in the master list have the same length, False otherwise.
	"""
	# The all() function checks if all values in an iterable are True and returns True, else returns False.
	# This function checks if the length of each list in the master list is equal to the length of the first list.
	assert all(len(lst) == len(master_list[0]) for lst in master_list)


def round_to_sig_figs(x, n_sig_figs=2):
	"""
	Round the number to a specified number of significant figures.

	Parameters
	----------
	x : float
		The number to round
	n_sig_figs : int, optional
		The number of significant figures to round to, default is 2

	Returns
	-------
	float
		The number rounded to n_sig_figs significant figures
	"""
	return np.round(x, n_sig_figs-int(np.floor(np.log10(np.abs(x))))-1)


def get_wcs_object_from_3D_fits_header(header_object):
	"""
	Get WCS object from a 3D FITS header.

	Parameters
	----------
	header_object : astropy.io.fits.Header
		The FITS header

	Returns
	-------
	astropy.wcs.WCS
		The World Coordinate System (WCS) object
	"""
	from astropy.wcs import WCS
	wcs = WCS(header_object)
	wcs = wcs.dropaxis(2)
	return wcs


# def smooth(y, box_pts):
# 	"""
# 	Smooth the data using a boxcar window.
# 
# 	Parameters
# 	----------
# 	y : numpy.ndarray
# 		The input array
# 	box_pts : int
# 		The number of points in the boxcar window
# 
# 	Returns
# 	-------
# 	numpy.ndarray
# 		The smoothed array
# 	"""
# 	box = np.ones(box_pts)/box_pts
# 	y_smooth = np.convolve(y, box, mode='same')
# 	return y_smooth

def smooth(y, box_pts):
	"""
	Smooth the data using a boxcar window. Includes dynamic window sizing at edges to mitigate
	edge effects.

	Parameters
	----------
	y : numpy.ndarray
		The input array
	box_pts : int
		The number of points in the boxcar window

	Returns
	-------
	numpy.ndarray
		The smoothed array
	"""
	y_smooth = np.copy(y)
	half_window = box_pts // 2

	# Apply smoothing with full window size for the central part of the array
	box = np.ones(box_pts) / box_pts
	y_smooth[half_window:-half_window] = np.convolve(y, box, mode='valid')

	# Adjust window size for the edges
	for i in range(half_window):
		box_size = half_window + i + 1
		box = np.ones(box_size) / box_size
		y_smooth[i] = np.convolve(y[:box_size], box, mode='valid')[-1]
		y_smooth[-(i+1)] = np.convolve(y[-box_size:], box, mode='valid')[0]

	return y_smooth


# def outlier_killer(points,thresh=5):
# 	"""
# 	Identify outliers using the Modified Z-score method.
# 
# 	Parameters
# 	----------
# 	points : numpy.ndarray
# 		An array of observations
# 	thresh : float, optional
# 		The modified z-score to use as a threshold, default is 5
# 
# 	Returns
# 	-------
# 	numpy.ndarray
# 		A boolean array, True if points are outliers and False otherwise
# 	"""
# 	if len(points.shape) == 1:
# 		points = points[:,None]
# 	median = np.median(points, axis=0)
# 	diff = np.sum((points - median)**2, axis=-1)
# 	diff = np.sqrt(diff)
# 	med_abs_deviation = np.median(diff)
# 
# 	modified_z_score = 0.6745 * diff / med_abs_deviation
# 
# 	return modified_z_score > thresh


def is_RM_outlier(rms, specific_rm, thresh=3, tolerance=1e-5):
	"""
	Check if a specific RM is an outlier using the sigma clipping method.

	Parameters
	----------
	rms : numpy.ndarray
		An array of RM values
	specific_rm : float
		The value of the specific RM to check
	thresh : float, optional
		The number of standard deviations to use for the clipping threshold, default is 3
	tolerance : float, optional
		The tolerance level for comparing floating-point RM values, default is 1e-5 rad/m/m

	Returns
	-------
	bool
		True if the specific RM is an outlier, False otherwise
	"""
	# Apply sigma clipping to the RMs
	clipped_data = sigma_clip(rms, sigma=thresh, maxiters=None, cenfunc=np.median)

	# Check if the specific RM is among the masked (outlier) values
	is_outlier = np.any(clipped_data.mask & (np.abs(rms - specific_rm) < tolerance))

	return is_outlier


def is_RM_outlier_zscore(rms, rm_errors, specific_rm, thresh=20, tolerance=1e-5):
	"""
	Check if a specific RM is an outlier using the sigma clipping method, considering uncertainties.

	Parameters
	----------
	rms : numpy.ndarray
		An array of RM values
	rm_errors : numpy.ndarray
		An array of uncertainties associated with the RM values
	specific_rm : float
		The value of the specific RM to check
	thresh : float, optional
		The number of standard deviations to use for the clipping threshold, default is 3
	tolerance : float, optional
		The tolerance level for comparing floating-point RM values, default is 1e-5 rad/m/m

	Returns
	-------
	bool
		True if the specific RM is an outlier, False otherwise
	"""
	# Normalize the RM values by their uncertainties
	normalized_rms = rms / rm_errors

	# Apply sigma clipping to the normalized RMs
	clipped_data = sigma_clip(normalized_rms, sigma=thresh, maxiters=None, cenfunc=np.median)

	# Determine the index of the specific RM
	idx = np.where(np.abs(rms - specific_rm) < tolerance)[0]

	# Check if the specific RM is among the masked (outlier) values
	is_outlier = clipped_data.mask[idx]

	return is_outlier


def read_casda_fits(file):
	"""
	Opens a FITS file and returns the header and data.

	Parameters:
	----------
	file (str): The filename of the FITS file.

	Returns:
	-------
	tuple: A tuple containing the header, the HDU (Header/Data Unit) list, and the parent sbid.
	"""
	from astropy.io import fits

	#print('Opening %s'%print(file))
	hdu = fits.open(file)
	header = hdu[0].header
	sbid = header['SBID'].strip()
	return header, hdu, sbid

def extract_coordinates(header):
	"""
	Extracts coordinates from a FITS header.

	Parameters:
	header (Header): The FITS header.

	Returns:
	tuple: A tuple containing two SkyCoord objects, the source and pointing centre coordinates.
	"""
	wcs = get_wcs_object_from_3D_fits_header(header)
	wx, wy = wcs.wcs_pix2world(1, 1, 1) 
	coord_obj_source = SkyCoord(wx,wy,unit=(u.deg, u.deg),frame='icrs')

	ra_pointing_centre = header['CRVAL1']
	dec_pointing_centre = header['CRVAL2']
	coord_obj_pointing_centre = SkyCoord(ra_pointing_centre,dec_pointing_centre,unit=(u.deg, u.deg),frame='icrs')

	return coord_obj_source, coord_obj_pointing_centre


def load_and_parse_fits(file, include_sbids_band1_POSSUM_full):

	"""
	Loads and parses FITS data from a file.

	Parameters
	----------
	file : str
		The file path of the FITS file.
	include_sbids_band1_POSSUM_full : list
		A list of SBIDs to include in the analysis.

	Returns
	-------
	tuple
		A tuple containing the SkyCoord objects for the source and the pointing centre, the pixel offsets in the RA and Dec directions, and the pixel scale in the mosaic. Returns None if the SBID is in the ignore list.
	"""

	# Read in the source FDF fits
	header, hdu, sbid = read_casda_fits(file)

	# Ignore this FDF if the parent SBID is in the ignore list
	if int(sbid) not in include_sbids_band1_POSSUM_full:
		#print('SBID is not in inclusion list. Returning None.')
		return None

	# Get ra, decl coordinates of source associated with this FDF
	wcs = get_wcs_object_from_3D_fits_header(header)
	wx, wy = wcs.wcs_pix2world(1, 1, 1) # The fits cubes are 1 x 1 pix in RA and dec, so calculate the position of that pixel
	coord_obj_source = SkyCoord(wx,wy,unit=(u.deg, u.deg),frame='icrs')

	# Get ra, decl mosaic reference value
	ra_pointing_centre = header['CRVAL1']
	dec_pointing_centre = header['CRVAL2']
	ra_direction_pix_offset_in_mosaic = header['CRPIX1'] 
	dec_direction_pix_offset_in_mosaic = header['CRPIX2'] 
	ra_pixel_scale_in_mosaic = np.abs(header['CDELT1']) 
	dec_pixel_scale_in_mosaic = np.abs(header['CDELT2'])
	coord_obj_pointing_centre = SkyCoord(ra_pointing_centre,dec_pointing_centre,unit=(u.deg, u.deg),frame='icrs')

	# Check if pixel scales are equal
	assert ra_pixel_scale_in_mosaic==dec_pixel_scale_in_mosaic,'The pixel deltas in RA and Decl are different!'

	return coord_obj_source, coord_obj_pointing_centre, ra_direction_pix_offset_in_mosaic, dec_direction_pix_offset_in_mosaic, ra_pixel_scale_in_mosaic, sbid, header, hdu

def load_and_parse_noise_spectra_fits(file):

	"""
	Loads and parses FITS data from a file that contains QU noise.

	Parameters
	----------
	file : str
		The file path of the FITS file.

	Returns
	-------
	tuple
		A tuple containing the SkyCoord objects for the source and the pointing centre, the pixel offsets in the RA and Dec directions, and the pixel scale in the mosaic. Returns None if the SBID is in the ignore list.
	"""

	# Read in the source FDF fits
	header, hdu, sbid = read_casda_fits(file)

	# Get noise spectrum
	noise_spec = hdu[0].data.flatten()

	# Return it
	return noise_spec

def calculate_nearest_beam(beam_offsets_in_mosaics_degs, ra_direction_pix_offset_in_mosaic, dec_direction_pix_offset_in_mosaic, pixel_scale_in_mosaic):
	"""
	Calculates the nearest ASKAP beam number.

	Parameters
	----------
	beam_offsets_in_mosaics_degs : numpy.ndarray
		The beam offsets in degrees in the mosaics.
	ra_direction_pix_offset_in_mosaic : float
		The pixel offset in the RA direction in the mosaic.
	dec_direction_pix_offset_in_mosaic : float
		The pixel offset in the Dec direction in the mosaic.
	pixel_scale_in_mosaic : float
		The pixel scale in the mosaic.

	Returns
	-------
	list
		A list containing the beam number, the distance to the nearest beam, and the offsets in the X and Y directions.
	"""
	# Init variables
	nearest_beam_recorder = [] #init
	nearest_beam_dist = 9e9 # Init

	for row_idx,row in enumerate(beam_offsets_in_mosaics_degs):
		offsets_pix = (row[1:]/pixel_scale_in_mosaic).astype(int) #Nearest pixel is fine for our purposes
		xoff = np.abs(ra_direction_pix_offset_in_mosaic-offsets_pix[0])
		yoff = np.abs(dec_direction_pix_offset_in_mosaic-offsets_pix[1])
		cartesian_distance_from_current_source_pixels = np.sqrt( xoff**2 + yoff**2)
		if cartesian_distance_from_current_source_pixels < nearest_beam_dist:
			nearest_beam_recorder = [row[0],cartesian_distance_from_current_source_pixels,xoff,yoff]
			nearest_beam_dist = cartesian_distance_from_current_source_pixels #update nearest beam

	return nearest_beam_recorder

def do_skip_RMs_in_bad_positions(sbid,coord_obj_source):

	"""
	Checks whether the current source should be skipped if the RM falls near certain bright
	/resolved objects that could cause corruption of QU spectra. These were determined by
	manual inspection of EMU MFS images.

	Parameters
	----------
	sbid : int
		The scheduling block ID (SBID) for the source.
	coord_obj_source : astropy.coordinates.sky_coordinate.SkyCoord
		The SkyCoord object for the source.

	Returns
	-------
	bool
		True if the source should be skipped, False otherwise.
	"""
	
	### Any of these individually is cause to reject if True, so can structure as a series of if...elif.. etc.

	#LMC
	LMC = SkyCoord(79.0*u.deg, -69.3*u.deg)
	separation_LMC = coord_obj_source.separation(LMC).deg
	if separation_LMC < 10:
		return True
	#M83
	elif int(sbid)==34120 and ( coord_obj_source.separation(SkyCoord('13:37:00.50145 -29:51:49.4944',unit=(u.hourangle, u.deg),frame='icrs')).deg < 0.1): 
		return True
	#LMC
	elif int(sbid)==46978 and ( coord_obj_source.dec.deg > -72): 
		return True
	#LMC again, as well as a very bright point source in SB46957
	elif int(sbid)==46957 and ( (coord_obj_source.dec.deg > -72 and coord_obj_source.ra.deg > 71) or (coord_obj_source.separation(SkyCoord('04:08:48.74843 -75:07:25.8386',unit=(u.hourangle, u.deg),frame='icrs')).deg < 1) ): 
		return True
	#V. bright uncleaned sidelobes in north of SB 46962
	elif int(sbid)==46962 and ( coord_obj_source.dec.deg > -68): 
		return True	
	#Bright side-lobey source in SW corner of SB 46982
	elif int(sbid)==46982 and ( coord_obj_source.separation(SkyCoord('02:52:44.58844 -71:04:38.8836',unit=(u.hourangle, u.deg),frame='icrs')).deg < 1): 
		return True	
	#In SB 47034, two brightish sources in the east, and one in the west
	elif int(sbid)==47034 and ( (coord_obj_source.separation(SkyCoord('02:43:28.99679 -51:08:44.4127',unit=(u.hourangle, u.deg),frame='icrs')).deg < 0.5) or (coord_obj_source.separation(SkyCoord('02:10:46.48579 -51:01:04.4853',unit=(u.hourangle, u.deg),frame='icrs')).deg < 0.3) ): 
		return True			
	#Bright, v. resolved radio galaxy in SB 46986
	elif int(sbid)==46986 and ( coord_obj_source.separation(SkyCoord('01:33:58.10282 -36:29:50.3988',unit=(u.hourangle, u.deg),frame='icrs')).deg < 1): 
		return True	
	else:
		return False

def parse_polarisation(header, FDF, band_average_stokes_noise, c=None, freq_lo_hz=799990740.7407, freq_hi_hz=1086990740.7407):
	"""
	Given the FDF FITS header and frequency range (the latter is not in the former, preposterously), 
	calculate the Faraday depth sampling and full width half maximum (FWHM) of the RMSF. 
	Additionally, it extracts RM, PSNR, uncertainties, and other parameters from the FITS 
	HDU.

	Parameters:
	header : dict
		A header from an FDF FITS file containing the Faraday depth information.
	FDF : Complex Array
		A complex-valued FDF
	c : float, optional
		The speed of light in m/s (default is 3e8 m/s).
	freq_lo_hz : float, optional
		The lower bound of the frequency range in Hz (default is 799990740.7407 Hz).
	freq_hi_hz : float, optional
		The upper bound of the frequency range in Hz (default is 1086990740.7407 Hz).

	Returns:
	tuple
		The Faraday depth sampling, the FWHM of the RMSF, and a dictionary containing the 
		fitted peak Faraday depth, its uncertainty, and the fitted signal-to-noise ratio.
	"""
	if c is None:
		from scipy.constants import c
	crval3 = header['CRVAL3'] #cube Faraday depth reference
	cdelt3 = header['CDELT3'] #cube Faraday depth reference
	naxis3 = header['NAXIS3'] #number of faraday depth samples along FD axis
	phiArr = np.array([crval3+(i*cdelt3) for i in range(naxis3)]) #FD sampling

	#band_average_stokes_noise = 15e-6

	lsq_freq_lo = (c/freq_lo_hz)**2
	lsq_freq_hi = (c/freq_hi_hz)**2
	fwhmRMSF = (2*np.sqrt(3))/(lsq_freq_lo-lsq_freq_hi) 

	lamSqArr_m2 = (c/np.array([7.999907407407E+08+(1e6*i) for i in range(288)]))**2

	#FDF = np.squeeze(hdu[0].data)
	FDF_params_dict = measure_FDF_parms(FDF, phiArr, fwhmRMSF, dFDF=band_average_stokes_noise, lamSqArr_m2=lamSqArr_m2,lam0Sq=np.nanmean(lamSqArr_m2), snrDoBiasCorrect=5.0)

	return phiArr, fwhmRMSF, FDF_params_dict


def from_J_format(coord_string):
	"""
	Converts coordinates from J format (e.g., "J000022.6+392944")
	to a more standard format (e.g., "00:00:22.6 +39:29:44").

	Parameters:
	- coord_string (str): Input string in J format.

	Returns:
	- str: Formatted string.
	"""

	# Ensure the string starts with 'J'
	if not coord_string.startswith("J"):
		raise ValueError("Input string does not have the correct format.")

	# Parse the string for RA and Dec
	ra_string = coord_string[1:9]  # Extracts "000022.6"
	dec_string = coord_string[9:]  # Extracts "+392944"

	# Convert to the desired format
	ra_formatted = f"{ra_string[0:2]}:{ra_string[2:4]}:{ra_string[4:]}"
	dec_formatted = f"{dec_string[0:3]}:{dec_string[3:5]}:{dec_string[5:]}"

	return f"{ra_formatted} {dec_formatted}"



def ra_dec_to_theta_phi(skycoord_column):
	"""
	For the Hutschenreuter Bayesian foreground fitting. 
	Converts a column of SkyCoord objects representing Right Ascension (RA) and Declination (Dec)
	to arrays of colatitude (theta) and longitude (phi) in radians, as per HealPy coordinate conventions.
	See notes in python/RM_interpol/src/Functions/data.py for more.

	Parameters:
	skycoord_column (Column): Astropy table column containing SkyCoord objects.

	Returns:
	array: An array of tuples containing corresponding colatitude (theta) and longitude (phi) in radians.
	"""

	theta_phi_coords = []
	for position in skycoord_column:
		ra_radians = position.ra.to_value(u.rad)  # Convert RA to radians
		dec_radians = position.dec.to_value(u.rad) # Convert Dec to radians
		theta = np.pi / 2 - dec_radians
		phi = ra_radians
		theta_phi_coords.append((theta, phi))

	return np.array(theta_phi_coords)


def ra_dec_to_coords_packager(skycoord_column):
	"""
	For the Hutschenreuter Bayesian foreground fitting. 
	Converts a column of SkyCoord objects representing Right Ascension (RA) and Declination (Dec)
	to arrays of those coords, for input into the inference engine.

	Parameters:
	skycoord_column (Column): Astropy table column containing SkyCoord objects.

	Returns:
	array: An array of tuples containing corresponding ra and dec in degrees.
	"""

	ra_dec_deg_coords = []
	for position in skycoord_column:
		ra_deg = position.ra.to_value(u.deg)  # Convert RA to radians
		dec_deg = position.dec.to_value(u.deg) # Convert Dec to radians
		ra_dec_deg_coords.append((ra_deg, dec_deg))

	return np.array(ra_dec_deg_coords)


def generate_hutsch_param_file(file_name, base_name_out, centre_list, random_seed):

	"""
	Writes text to a file, including specified base name and center list, in the format
	demanded by the Hutschenreuter inference engine.

	:param file_name: str, the name of the file where the text will be written
	:param base_name_out: str, the base name to be included in the text
	:param centre_list: list, the center list to be included in the text
	"""

	text = f"""
base_name_out = '{base_name_out}'
centre_list = {centre_list}

run_params = {{
	'seed': {random_seed},
	'run_name': base_name_out,
	'use_mock_data': False,
	'data_name': base_name_out,
	'do_plot': True,
	'n_iterations': 15
}}

domain_params = {{
	'full_sphere': False,
	'nside': 128,
	'nx': 175,
	'ny': 175,
	'dx': 0.04,
	'dy': 0.04,
	'center': centre_list,
}}

# Model parameters

# Parameters for the log-amplitude sky map
amplitude_params = {{# prior parameters on the shape of the power spectrum ([mean, std])
					'asperity': [.001, 0.001], 'flexibility': [.001, 0.001], 'fluctuations': [2., .5], 'loglogavgslope': [-4., .5],
					# Parameters on the mean of the log amplitude field
					'offset_mean': 5., 'offset_std': [3., 2.]
					}}

# Parameters for the sign sky map
sign_params = {{# Prior parameters on the shape of the power spectrum ([mean, std])
			   'asperity': [.001, 0.001], 'flexibility': [1., .5], 'fluctuations': [1., 1.], 'loglogavgslope': [-4, .5],
				#  Parameters on the mean of the sign field
				'offset_mean': 0, 'offset_std': [1., .5],
			  }}

extragal_params = {{'type': 'explicit', # Model for the generation of the extragalactic component, choices are explicit, implicit or <None>
				   'alpha': 1.0, 'q': 2.0 # Hyperparameters for the explicit and implicit cases
				   }}
"""

	with open(file_name, 'w') as file:
		file.write(text)


def read_values_from_hutsch_param_file(file_name):
	"""
	Reads a file and extracts the values of nx, ny, dx, dy, and center.

	:param file_name: str, the name of the file to read
	:return: tuple, the extracted values of nx, ny, dx, dy, and center
	"""
	import re

	# Define the patterns to search for
	nx_pattern = r"'nx':\s*(\d+)"
	ny_pattern = r"'ny':\s*(\d+)"
	dx_pattern = r"'dx':\s*([\d.]+)"
	dy_pattern = r"'dy':\s*([\d.]+)"
	center_pattern = r"centre_list\s*=\s*\[([\d.-]+),\s*([\d.-]+)\]"

	# Read the file content
	with open(file_name, 'r') as file:
		file_content = file.read()

	# Extract values using regular expressions
	nx = int(re.search(nx_pattern, file_content).group(1))
	ny = int(re.search(ny_pattern, file_content).group(1))
	dx = float(re.search(dx_pattern, file_content).group(1))
	dy = float(re.search(dy_pattern, file_content).group(1))
	center = [float(x) for x in re.search(center_pattern, file_content).groups()]

	return nx, ny, dx, dy, center


def hutschen_calc_pixels(colat_rad, lon_rad, centre, dx, dy, nx, ny):
	"""
	This is here to expose the pixel mapper from the 
	Hutschenreuter PlaneProjector operator. Calculates the 
	pixel numbers corresponding to the data points within the plane cutout.
	"""
	assert nx==ny, "nx must equal ny!"
	assert dx==dy, "dx must equal dy!"
	distances = np.array([dx,dy])
	shape = np.array([nx,ny])
	lower_corner = centre - np.asarray(distances)*np.asarray(shape)*0.5
	#print(lower_corner)
	indices = []
	pixels = []
	for i in range(len(colat_rad)):
		pix_colat = np.floor((colat_rad[i] - lower_corner[0])/distances[0])
		pix_lon = np.floor((lon_rad[i] - lower_corner[1])/distances[1])
		#print([pix_colat, pix_lon])
		if (0 <= pix_colat < nx) & (0 <= pix_lon < ny):
			pixels.append([pix_colat, pix_lon])
			indices.append(i)
		else:
			pixels.append([np.nan, np.nan])
			indices.append(i)
	return np.asarray(pixels), indices
