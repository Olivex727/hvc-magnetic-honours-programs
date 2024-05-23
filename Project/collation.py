#!/usr/bin/env python3

# TODO:
# - Find scaling constant to convert from physical sizes to k-space
# - Find HVCs in RM range
# - Test RM function
# - Test the corner functions

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')

import numpy as np
import copy

from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table

from astropy.coordinates import SkyCoord, Angle
from astropy.io.votable import parse_single_table

from astropy.units import astrophys as astru
from astropy import units as u

from foreground import foreground_remover as fgrm, interpolate as intp
from wcs import image_transform as it

from plotting import honours_plot as hplt

import warnings
from astropy.wcs import FITSFixedWarning as fitswarn
verwarn = fits.verify.VerifyWarning

class file_find:

    def get_RMs(file='../data_catalog/main_table_17May2024', ext=".fits"):
        if ext == ".ecsv": return ascii.read(file+ext) #"../data_catalog/CASDA_RMs.ecsv"
        if ext == ".fits":
            # Must manually add the ra_dec_obj column
            t = Table(fits.open(file+ext)[1].data)
            t.remove_column('ra_dec_obj.ra')
            t.remove_column('ra_dec_obj.dec')
            t['RM'].unit = u.rad / (u.m ** 2)
            t['RM_uncert'].unit = u.rad / (u.m ** 2)
            tobj = []
            l = len(t)
            for rm_index in range(len(t)):
                rm_point = t[rm_index]['ra_dec_deg']
                tobj.append(SkyCoord(ra=rm_point[0]*u.degree, dec=rm_point[1]*u.degree, frame='icrs'))
                print(str(int(rm_index/l*100))+"% \r", sep="", end="", flush=True)
            print("Converting coordinates")
            t.add_column(tobj, 0,'ra_dec_obj')
            return t

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
        return {"interpolation":interpolation, "error":error, "k-space":fgrm.get_k_space(interpolation.data)}

class collator:

    # WARNING: Collating the RMs from scratch may take some time, make sure to specify a file to save the list in, and load on all future usages of this function, so that it only needs to run once.
    def data_whole_sky(calculate_interpolation, hvc_area_range=(1, np.pi), full_hvc_range=False, save_data="", load_data="", h1_img="../data_catalog/hi4pi-hvc-nhi-car.fits"):
        with warnings.catch_warnings(action="ignore", category=verwarn) and warnings.catch_warnings(action="ignore", category=fitswarn):
            print("=== WHOLE-SKY DATA COLLATION ===")
            print("Gathering data ...")
            print("Getting H-alpha emission")
            H_alpha = file_find.get_H_alpha()
            print("Extracting RMs")
            RMs = 0
            if load_data:
                #"../data_processed/proc_rms.ecsv"
                print("Collating RMs")
                RMs = file_find.get_RMs(load_data, ".ecsv")
            else:
                RMs_raw = file_find.get_RMs()
                print("Collating RMs")
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
            return {"RMs":RMs, "HVCs":HVCs, "HI":HIem, "H-alpha":H_alpha[0], "interpolation":interp}
    
    def collate(RMs, Ha_img, Ha_err):
        Ha_eco = []
        Ha_col = []
        l = len(RMs)
        for row in range(l):
            ha = it.get_flux_at_point(Ha_img, RMs[row]['ra_dec_obj'])
            he = it.get_flux_at_point(Ha_err, RMs[row]['ra_dec_obj'])

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

    # IMPORTANT: fourier transforms are linear, thus the uncertainties in the interpolation remain unchanged after filtering
    def take_snapshot(hvc_index, RMs, HVCs, HIem, H_alpha, interp, custom_selection=False, hvc_area_range=(1, np.pi), plot=False):
        with warnings.catch_warnings(action="ignore", category=fitswarn):
            print("=== HVC SNAPSHOT ===")
            print("Gathering data ...")
            if not custom_selection:
                selected_HVC = HVCs[hvc_index]
            else:
                selected_HVC = custom_selection
                print("Determining corners")
            corners = hvc_snapshot.get_corners(selected_HVC)
            if not corners:
                print("Could not resolve - HVC is on edge of sky")
                return 0
            print("Cropping H-alpha")
            Ha = hvc_snapshot.crop_wcs(corners, H_alpha, index=hvc_index, plot=plot)
            print("Cropping HI")
            H1 = hvc_snapshot.crop_wcs(corners, HIem, index=hvc_index, plot=plot)
            print("Cropping interpolation")
            intp_map = hvc_snapshot.crop_wcs(corners, interp["interpolation"], index=hvc_index, plot=plot)
            print("Filtering RMs")
            RMs_filtered = hvc_snapshot.rm_filter(corners, RMs)
            print("Correcting foreground")
            intp_cor, RMs_filtered = hvc_snapshot.foreground_correction(corners, interp["interpolation"], interp["k-space"], interp["error"], RMs_filtered, hvc_area_range, plot=plot, index=hvc_index)
            print("Snipping complete")
            return {"corners":corners, "H-alpha":Ha, "HI":H1, "interpolation_raw":intp_map, "interpolation_corrected":intp_cor, "RMs":RMs_filtered}
    
    # FIXME: Test if correct
    def get_corners(selected_HVC):
        dx = selected_HVC['dx']
        dy = selected_HVC['dy']

        # Twice the area of the HVC to ensure it's including all of the HVCs
        di = max(dx, dy) * u.deg

        centre_coord = selected_HVC['SkyCoord'].galactic

        # Calculate upper corner coordinate
        new_coord_up = SkyCoord(centre_coord.l+di, centre_coord.b+di, frame='galactic')
        new_coord_down = SkyCoord(centre_coord.l-di, centre_coord.b-di, frame='galactic')

        if new_coord_down.l > new_coord_up.l or new_coord_down.b > new_coord_up.b:
            return False
        
        return [new_coord_down, new_coord_up, centre_coord]
    
    def crop_wcs(corners, image, plot=False, index=-1):

        # Get corners in terms of images
        wcs = it.wcs(image.header)
        pix_down = list(np.array(list(map(int, it.get_pixel(wcs, corners[0]))))-1)
        pix_up = list(np.array(list(map(int, it.get_pixel(wcs, corners[1]))))-1)
        if len(corners) > 2: pix_c = list(np.array(list(map(int, it.get_pixel(wcs, corners[2]))))-1)
        else: pix_c = [0,0]

        def crop_img_1d(img, pixel_up, pixel_down):
            img = img[pixel_down[1]:pixel_up[1]]
            return img

        def crop_img(img, pixel_up, pixel_down):
            img = img[pixel_down[1]:pixel_up[1]] if pixel_down[1]<pixel_up[1] else img[pixel_up[1]:pixel_down[1]]
            img = np.transpose(img)
            img = img[pixel_down[0]:pixel_up[0]] if pixel_down[0]<pixel_up[0] else img[pixel_up[0]:pixel_down[0]]
            img = np.transpose(img)
            return img
        
        if plot:
            print(pix_up)
            print(pix_down)
            hplt.plot_image_crop(image, crop_img(image.data, pix_up, pix_down), pix_up, pix_down, pix_c, index)
        
        return crop_img(image.data, pix_up, pix_down)

    # FIXME: Test for validity
    def rm_filter(corners, RMs):

        # Create mask
        gal_RM_locations = RMs["ra_dec_obj"].galactic
        mask = list(map(lambda rm_loc: corners[0].l < rm_loc.l < corners[1].l and corners[0].b < rm_loc.b < corners[1].b, gal_RM_locations))

        print(str(sum(mask))+" RM grid points found")

        # Filter the RMs
        filtered_RMs = copy.deepcopy(RMs[mask])

        return filtered_RMs
    
    # Assumes RMs are already corrected
    def foreground_correction(corners, interpolation, k_space, interpolation_std, RMs, hvc_area_range=(1, np.pi), index=-1, plot=False):

        # Filter the k-space and get corrected foreground
        new_k_space = fgrm.filter_k_space(k_space, hvc_area_range)
        cor_fg = fgrm.restore_foreground(new_k_space, interpolation)

        # Correct RMs
        cor_RMs = fgrm.get_corrected_RMs(interpolation, cor_fg, interpolation_std, RMs)

        # Take snapshots of corrected interpolation
        snap_cor = hvc_snapshot.crop_wcs(corners, cor_fg, index=index, plot=plot)

        return snap_cor, cor_RMs

class collation_tools:

    # WARNING: This function takes a long time to execute, make sure to specify a file to save the list in, so that it only needs to run once.
    def calculate_HVC_seperation(collated_data, save_file = "../output.txt"):
        l = len(collated_data['HVCs']) - 1
        lr = len(collated_data['RMs']) - 1
        row = 0
        row2 = 0
        val = 0

        sep_arr = []

        for hvc in collated_data['HVCs']:
            val = int(row/l*100)
            row2 = 0

            minimum = -1
            for rm in collated_data['RMs']:
                new = hvc['SkyCoord'].separation(rm['ra_dec_obj'])
                if minimum == -1 or minimum > new:
                    minimum = new
        
                row2 = row2 + 1
                print("RM loop: "+str(int(row2/lr*100))+"%; HVC loop: "+str(val)+"% \r", sep="", end="", flush=True)

            sep_arr.append(minimum)

            print("RM loop: 100%; HVC loop: "+str(val)+"% \r", sep="", end="", flush=True)
            row = row + 1

        with open("output.txt", "w") as txt_file:
            for line in sep_arr:
                txt_file.write(str(line)+"\n")
        
    def add_HVC_seperation(hvcs, max_sep=2*np.pi, load_file="../output.txt"):
        sep_arr = []

        with open("output.txt", "r") as txt_file:
            sep_arr = txt_file.read().split("\n")[:-1]

        sep_arr = list(map(Angle, sep_arr))

        hvcs = copy.deepcopy(hvcs)

        hvcs.add_column(sep_arr, name="Nearest RM")

        hvcs = hvcs[hvcs["Nearest RM"].value < max_sep]

        return hvcs