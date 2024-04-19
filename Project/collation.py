#!/usr/bin/env python3

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')

import numpy as np
import math

from astropy.io import ascii

from astropy.wcs import WCS
from astropy.io import fits

from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table

from astropy.units import astrophys as astru

from foreground import foreground_remover as fr, interpolate as intp
from wcs import image_transform as it

class file_find:

    def get_RMs(file="../data_catalog/CASDA_RMs.ecsv"):
        return ascii.read(file)

    def get_HI_emission(h1_img):
        hdu = fits.open(h1_img)[0]
        return hdu
    
    def get_HVC_locations(hvc_area_range=(1, np.pi), full_hvc_range=False):
        def format_sexagesimal(coord):
            ra, dec = coord.split()[:3], coord.split()[3:]
            ra_formatted = ':'.join(ra)
            dec_formatted = f"{dec[0]}:{dec[1]}:{dec[2]}"
            return f"{ra_formatted} {dec_formatted}"

        # Read the VOTable file
        table = parse_single_table('../data_catalog/vizier_Moss2013_HVCs.vot').to_table()

        # Mask to decent sized area
        if (not full_hvc_range):
            mask = (table['Area']>hvc_area_range[0]) & (table['Area']<hvc_area_range[1])
            table = table[mask]

        # Assuming the coordinates are in columns 'ra' and 'dec', adjust as needed
        ra_dec = [f"{row['RAJ2000']} {row['DEJ2000']}" for row in table]
        ra_dec = [format_sexagesimal(f"{row['RAJ2000']} {row['DEJ2000']}") for row in table]

        ra_dec = SkyCoord(ra_dec,frame='fk5',unit=('hour','deg')).icrs

        # Data is now in desired structure
        table.add_column(ra_dec, index=1, name="SkyCoord")

        return table

    def get_H_alpha():
        img = fits.open("../data_catalog/Halpha_map.fits")[0]
        err = fits.open("../data_catalog/Halpha_error.fits")[0]
        return img, err
    
    def get_interpolation(RMs=0, use_H_alpha=True, H_alpha=0, calculate_interpolation=True):
        interpolation, error = intp.interpolate(RMs, use_H_alpha, H_alpha, calculate_interpolation)
        return [interpolation, error, fr.get_k_space(interpolation.data)]

class collator:

    def data_whole_sky(calculate_interpolation, hvc_area_range=(1, np.pi), full_hvc_range=False, save_data="", load_data="", h1_img="../data_catalog/hi4pi-hvc-nhi-car.fits"):
        print("Gathering data ...")
        print("Getting H-alpha emission")
        H_alpha = file_find.get_H_alpha()
        print("Collating RMs")
        RMs = 0
        if load_data:
            #"../data_processed/proc_rms.ecsv"
            RMs = file_find.get_RMs(load_data+".ecsv")
        else:
            RMs_raw = file_find.get_RMs()
            RMs = collator.collate(RMs_raw, H_alpha[0], H_alpha[1])
        print("Getting HVC location data")
        HVCs = file_find.get_HVC_locations(hvc_area_range, full_hvc_range)
        print("Getting HI emission")
        HIem = file_find.get_HI_emission(h1_img)
        print("Interpolating")
        interp = file_find.get_interpolation(calculate_interpolation=calculate_interpolation)
        if save_data:
            print("Saving processed RM table")
            collator.write_processed(RMs, save_data)
        print("Collation complete")
        return RMs, HVCs, HIem, H_alpha[0], interp
    
    def collate(RMs, Ha_img, Ha_err):
        Ha_eco = []
        Ha_col = []
        l = len(RMs)
        for row in range(l):
            ha = it.get_flux_at_point(Ha_img, RMs[row][0])
            he = it.get_flux_at_point(Ha_err, RMs[row][0])

            # H-alpha data is in rayleighs. 1 R = 1e10 m-2 sr-1 sec-1
            Ha_col.append(ha * astru.R)
            Ha_eco.append(he * astru.R)

            print(str(int(row/l*100))+"% \r", sep="", end="", flush=True)

        RMs.add_column(Ha_col, index=1, name="H-alpha flux")
        RMs.add_column(Ha_eco, index=1, name="H-alpha flux [Error]")
        return RMs
    
    def write_processed(RMs, file):
        RMs.write(file+".ecsv", format='ascii.ecsv', overwrite=True)  

class hvc_snapshot:

    def take_snapshot(hvc_index, RMs, HVCs, HIem, H_alpha, interp):
        return 0
    
    def get_corners(hvc_index, HVCs):
        return 0
    
    def crop_wcs(corners, image):
        return 0
    
    def rm_filter(corners, RMs):
        return 0
    
    def foreground_correction(hvc_index, HVCs, interpolation, RMs):
        return 0
    