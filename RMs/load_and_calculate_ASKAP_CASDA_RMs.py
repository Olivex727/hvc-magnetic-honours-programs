#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
~MWE script for Olivia demonstrating read-in of ASKAP FDFs, calculate RMs.
"""

##########################################################################################
###################################### IMPORTS ###########################################
##########################################################################################

import os
import glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import warnings

#warnings.simplefilter('ignore', category=VerifyWarning)
warnings.filterwarnings('ignore', category=UserWarning, append=True)

##########################################################################################
######################### SETTINGS, PARAMS, FUNCTIONS ####################################
##########################################################################################

from group_proj_settings import *
from group_proj_parameters import *
from group_proj_functions import *
from group_proj_sbidness import *
#warnings.simplefilter('ignore', category=VerifyWarning)
warnings.filterwarnings('ignore', category=UserWarning, append=True)
#warnings.simplefilter('ignore', category=AstropyWarning)

##########################################################################################
###################################### INITS #############################################
##########################################################################################

#pl.ion()
#Init results arrays
radec_obj_accumulator=[]
radec_accumulator=[]

##########################################################################################
###################################### MAIN ##############################################
##########################################################################################


"""
Read in FDFs + Stokes spectra noise and calculate RMs for data outputted from the raw CASDA validation spectra.
"""

if True:
	# Get into data directory
	os.chdir(basedir+'POSSUM_pilot_askapsoft-pipeline_fdfs/fdf_amp_spectra')

	# List (FDF FITS files) in data directory
	files = glob.glob('./pol_FDF_amp_*.fits')
	files.sort()

	# For each FDF file...
	for fidx,file in enumerate(files):

		print('Loading in source #%d (of %d)'%(fidx,len(files)))

		# Load FITS file, and return results. Continue if current FDF comes from one of the SBIDs to ignore (handled by load_and_parse_fits returning 'None').
		result_FDF_amp = load_and_parse_fits(file, include_sbids_band1_POSSUM_full)
		result_FDF_phase = load_and_parse_fits('../fdf_phase_spectra/'+file.replace('amp','phase'), include_sbids_band1_POSSUM_full)

		# Parse out 'results' --- quantities that can be derived from the FITS headers.
		coord_obj_source, coord_obj_pointing_centre, ra_direction_pix_offset_in_mosaic, dec_direction_pix_offset_in_mosaic, pixel_scale_in_mosaic, sbid, FDF_header, hdu_amp = result_FDF_amp
		_, _, _, _, _, _, _, hdu_phase = result_FDF_phase
	
		# Create complex-valued FDF from amp, phase FDF data
		amplitudes = np.squeeze(hdu_amp[0].data)
		phases = np.squeeze(hdu_phase[0].data)
		complex_FDF = amplitudes * (np.cos(phases) + 1j * np.sin(phases))

		#Get QU noise
		BAN_recorder = None	
		for stokes in ['Q','U']:
			stokes_noise_spectrum_file = basedir+'POSSUM_pilot_askapsoft-pipeline_fdfs/iquv_noise_spectra/'+file.replace('pol_FDF_amp','pol_noise_%s'%stokes)
			noise_spectrum = load_and_parse_noise_spectra_fits(stokes_noise_spectrum_file)
			band_av_noise = np.nanmedian(noise_spectrum) / np.sqrt(len(noise_spectrum))
			if stokes == 'Q':
				BAN_recorder = band_av_noise
		band_av_noise_stokes_av = np.mean([band_av_noise,BAN_recorder])

		## Get polarisation goodies from source
		phiArr, fwhmRMSF, FDF_params_dict = parse_polarisation(FDF_header, complex_FDF, band_av_noise_stokes_av)
		phi_peak_fitted_radmm = FDF_params_dict['phiPeakPIfit_rm2']
		dphi_peak_fitted_radmm = FDF_params_dict['dPhiPeakPIfit_rm2']
		SN_PI_fitted = FDF_params_dict['snrPIfit']
		polint_PI_fitted = FDF_params_dict['ampPeakPIfit']
		polint_err_PI_fitted = FDF_params_dict['dAmpPeakPIfit']


		#Append results to accumulator lists if the source/tile meets certain criteria
		radec_obj_accumulator.append(coord_obj_source)
		radec_accumulator.append([coord_obj_source.ra.deg,coord_obj_source.dec.deg])

## Post-loop --- merge results into single array/table. Save to external file if we want.
radec_obj_accumulator_arr = np.array(radec_obj_accumulator)
radec_accumulator_arr = np.array(radec_accumulator)


## Create astropy table of results for individual spectra
#t = Table([radec_accumulator_arr, radec_pointing_centre_arr, gal_lat_lon_accumulator_arr, phi_peak_accumulator_arr, phi_peak_uncertainty_accumulator_arr, sbid_accumulator_arr, crpix_accumulator_arr, cdelt_accumulator_arr, nearest_beam_accumulator_arr, nearest_beam_accumulator_arr[:,0].astype(int), stokes_I_peak_flux_density_accumulator_arr, stokes_I_peak_flux_density_err_accumulator_arr,stokes_I_flag_accumulator_arr,stokes_I_has_siblings_accumulator_arr,polint_accumulator_arr,polint_uncertainty_accumulator_arr], names=('ra_dec_deg', 'ra_dec_pointing_centre_for_SBID_deg', 'gal_lat_lon_deg', 'faraday_depth_radmm', 'faraday_depth_err_radmm', 'sbid', 'crpix_mosaic_offsets', 'cdelts_mosaic', 'nearest_beam_info', 'nearest_beam_no','stokes_I_peak_flux_density','stokes_I_peak_flux_density_err','stokes_I_flag','has_stokes_I_siblings','pol_int','pol_int_err'))
table_data = {
	'ra_dec_obj': radec_obj_accumulator_arr,
	'ra_dec_deg': radec_accumulator_arr,
}

t = Table(table_data)

"""
Add descriptions to individual columns metadata
"""

t['ra_dec_obj'].description = 'RA/DEC object coordinates'
t['ra_dec_deg'].description = 'RA/DEC coordinates in degrees'

