#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
Functions to be read in by groups_basic_pipeline_output_fdf_analyser_ver10_collab_restructure.py
or later versions.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.constants import c
import astropy.units as u

# Only include this if measure_FDF_parms is indeed utilized in the script.
from RMutils.util_RM import measure_FDF_parms


##########################################################################################
###################################### FUNCITONS #########################################
##########################################################################################

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

def load_and_parse_fits(file):

	"""
	Loads and parses FITS data from a file.

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





