#!/usr/bin/env python3

from astropy.coordinates import SkyCoord

from astropy.wcs import WCS
import numpy as np

import h5py
from astropy.io import fits
import reproject as repo

class image_transform:
    def get_pixel(wcs, coords):
        car = 0
        if "galactic" in wcs.world_axis_physical_types[0]:
            coords = coords.galactic
            car = np.array([[coords.l.deg],[coords.b.deg]]).T
        else:
            coords = coords.icrs
            car = np.array([[coords.ra.deg],[coords.dec.deg]]).T

        return WCS.wcs_world2pix(wcs,car,1)[0]
    
    # Image must be in FITS format
    def get_flux_at_point(image, coords, healpix=False):
        if healpix:
            return 0
        else:
            wcs = WCS(image.header)
            pix = image_transform.get_pixel(wcs, coords)
            return image.data[int(pix[1])-1][int(pix[0])-1]
    
    # Read a h5py healpix file and write to a fits image
    # Outputs footprint of reprojection
    # IMPORTANT: function remains untested
    def healpy_to_fits(load_file, save_file):
        f = h5py.File('../data_catalog/'+load_file+'.hdf5', 'r')

        healpix_ringorder_N128_mean = f['faraday sky']['mean']
        healpix_ringorder_N128_std = f['faraday sky']['std']

        mean, footprint_mean = repo.reproject_from_healpix((np.array(healpix_ringorder_N128_mean), 'galactic'), target_header, nested=False)
        std, footprint_std = repo.reproject_from_healpix((np.array(healpix_ringorder_N128_std), 'galactic'), target_header, nested=False)

        hdu_mean = fits.PrimaryHDU(mean, target_header)
        hdu_mean.writeto("../data_catalog/"+save_file+"_mean.fits")

        hdu_std = fits.PrimaryHDU(std, target_header)
        hdu_std.writeto("../data_catalog/"+save_file+"_std.fits")

        return footprint_mean, footprint_std
    
    def wcs(header):
        return WCS(header)
    
target_header = fits.Header.fromstring("""
SIMPLE  =                    T /     
BITPIX  =                  -32 / Number of bits per data pixel                  
NAXIS   =                    2 / Number of data axes                            
NAXIS1  =                 8640 /                                                
NAXIS2  =                 4320 /                                                
DATE    = '2024-04-12'         / Creation UTC (CCCC-MM-DD) date of FITS header  
COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy  
COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H    
CTYPE1  = 'GLON-CAR'           /                                                
CRPIX1  =              4320.50 /                                                
CRVAL1  =              0.00000 /                                                
CTYPE2  = 'GLAT-CAR'           /                                                
CRPIX2  =              2160.50 /                                                
CRVAL2  =              0.00000 /                                                
CD1_1   =     -0.0416666666667 /                                                
CD1_2   =        0.00000000000 /                                                
CD2_1   =        0.00000000000 /                                                
CD2_2   =      0.0416666666667 /                                                
LONPOLE =                  180 /                                                
PROJP1  =                    0 /                                                
PROJP2  =                    0 /                                                
EQUINOX =              2000.00 /                                                
BUNIT   = 'R       '           /                                                
VERSION = '1.1     '           /   
""", sep='\n')