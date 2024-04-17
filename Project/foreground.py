#!/usr/bin/env python3

from astropy.io import fits
from scipy import fft

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
    
    def punch_annulus(base, inner_radius, outer_radius, centre=-1):
        if centre == -1:
            centre = (len(base)/2, len(base)/2)

        for x in range(len(base)):
            for y in range(len(base[x])):
                if inner_radius ** 2 < ((x-centre[0]) ** 2 + (y-centre[1]) ** 2) < outer_radius ** 2:
                    base[x][y] = 0
    
        return base
    
    def filter_k_space(interpolation):
        return 0
    
    def get_corrected_image(k_space):
        return 0
    
    def get_corrected_RMs(k_space, RMs):
        return 0