#!/usr/bin/env python3

from astropy.io import fits
from scipy import fft
import numpy as np
from dict2obj import Dict2Obj

from wcs import image_transform as it

class interpolate:
    def interpolate(RMs=0, use_H_alpha=True, H_alpha=0, calculate_interpolation=True):
        if calculate_interpolation:
            return 0
        else:
            if use_H_alpha:
                hdu_main = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_wff_mean.fits")[0]
                hdu_err = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_wff_std.fits")[0]
                return hdu_main, hdu_err
            else:
                hdu_main = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_woff_mean.fits")[0]
                hdu_err = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_woff_std.fits")[0]
                return hdu_main, hdu_err

# Return set of corrected RM points
class foreground_remover:
    def get_k_space(interpolation):
        return fft.fft2(interpolation)
    
    def punch_annulus(base, inner_radius, outer_radius, centre=False):
        if not centre:
            centre = (len(base[0])/2, len(base)/2)
        
        for y in range(len(base)):
            for x in range(len(base[y])):
                if inner_radius ** 2 < ((x-centre[0]) ** 2 + (y-centre[1]) ** 2) < outer_radius ** 2:
                    base[y][x] = 0
    
        return base
    
    def filter_k_space(k_space, size_params=(1,np.pi)):
        annulus_space = (k_space.real * 0) + 1
        annulus_space = foreground_remover.punch_annulus(annulus_space, scale(size_params[0]), scale(size_params[1]))
        new_k_space = annulus_space*k_space
        return new_k_space
    
    def get_corrected_image(k_space):
        return fft.ifft2(k_space)
    
    # Assumes a pre-filtered k-space
    def restore_foreground(k_space, interpolation):
        corr_fg = Dict2Obj({"header":interpolation.header, "data":foreground_remover.get_corrected_image(k_space).real})
        return corr_fg
    
    def get_corrected_RMs(interpolation_pre, interpolation_post, interpolation_std, RMs):
        pre_fg_list = []
        pst_fg_list = []
        std_fg_list = []

        for entry in RMs:
            coords = entry["ra_dec_obj"].galactic
            strength = it.get_flux_at_point(interpolation_pre, coords)
            pre_fg_list.append(strength)
            strength = it.get_flux_at_point(interpolation_post, coords)
            pst_fg_list.append(strength)
            strength = it.get_flux_at_point(interpolation_std, coords)
            std_fg_list.append(strength)

        RMs.add_column(pre_fg_list, name="interpolation_raw")
        RMs.add_column(pst_fg_list, name="interpolation_cor")
        RMs.add_column(std_fg_list, name="interpolation_unc")
        return RMs
    
# Define a scaling function to convert real size -> k-space size
def scale(value):
    scl = 3000
    return value * scl