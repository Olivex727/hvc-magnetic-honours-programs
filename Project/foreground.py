#!/usr/bin/env python3

from astropy.io import fits
from scipy import fft
import numpy as np
from dict2obj import Dict2Obj

class interpolate:
    def interpolate(RMs=0, use_H_alpha=True, H_alpha=0, calculate_interpolation=True):
        if calculate_interpolation:
            return 0
        else:
            if use_H_alpha:
                hdu_main = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_wff_mean.fits")[0]
                hdu_err = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_wff_std.fits")[0]
                return hdu_main, hdu_err, foreground_remover.get_k_space(hdu_main.data), foreground_remover.get_k_space(hdu_err.data)
            else:
                hdu_main = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_woff_mean.fits")[0]
                hdu_err = fits.open("../data_preprocessed/Hutschenreuter_2020_faraday_sky_woff_std.fits")[0]
                return hdu_main, hdu_err, foreground_remover.get_k_space(hdu_main.data), foreground_remover.get_k_space(hdu_err.data)
    
    def fourier_interpolate(interpolation, k_space, hvc_area_range):
        filtered_k_space = foreground_remover.filter_k_space(k_space, hvc_area_range)
        return foreground_remover.restore_foreground(filtered_k_space, interpolation)

# Return set of corrected RM points
class foreground_remover:
    def get_k_space(interpolation):
        return fft.fft2(interpolation)
    
    def get_freq_domains(k_space):
        fy = fft.fftfreq(k_space.shape[0])
        fx = fft.fftfreq(k_space.shape[1])
        return fx, fy
    
    def punch_annulus(base, inner_radius, outer_radius, centre=False):
        if not centre:
            centre = (len(base[0])/2, len(base)/2)
        
        for y in range(len(base)):
            for x in range(len(base[y])):
                if inner_radius ** 2 < ((x-centre[0]) ** 2 + (y-centre[1]) ** 2) < outer_radius ** 2:
                    base[y][x] = 0
    
        return base
    
    def convert_size_to_frequency(theta_hvc, l_shape=8640):
        R = l_shape / 360
        f_hvc = 1/(2*theta_hvc*R)
        return f_hvc

    def convert_frequency_to_pixel(f_range, fftfreq):
        fftfreq_range = list(map(lambda x: f_range[0] < x < f_range[1], fftfreq))
        return fftfreq_range

    def punch_crosshatch(base, x_mask, y_mask):
        for y in range(len(base)):
            for x in range(len(base[y])):
                if x_mask[x] or y_mask[y]:
                    base[y][x] = 0
    
        return base
    
    def filter_k_space(k_space, size_params=(1,np.pi)):
        fx, fy = foreground_remover.get_freq_domains(k_space)

        crosshatch_space = ((k_space.real * 0) + 1)
        hvc_f_range_pos = tuple(map(foreground_remover.convert_size_to_frequency, size_params))[::-1]
        hvc_f_range_neg = tuple(map(lambda x: - x, map(foreground_remover.convert_size_to_frequency, size_params)))

        x_range = foreground_remover.convert_frequency_to_pixel(hvc_f_range_pos, fx)
        y_range = foreground_remover.convert_frequency_to_pixel(hvc_f_range_pos, fy)

        crosshatch_space = foreground_remover.punch_crosshatch(crosshatch_space, x_range, y_range)

        x_range = foreground_remover.convert_frequency_to_pixel(hvc_f_range_neg, fx)
        y_range = foreground_remover.convert_frequency_to_pixel(hvc_f_range_neg, fy)

        crosshatch_space = foreground_remover.punch_crosshatch(crosshatch_space, x_range, y_range)
        new_k_space = crosshatch_space*k_space
        return new_k_space
    
    def get_corrected_image(k_space):
        return fft.ifft2(k_space)
    
    # Assumes a pre-filtered k-space
    def restore_foreground(k_space, interpolation):
        corr_fg = Dict2Obj({"header":interpolation.header, "data":foreground_remover.get_corrected_image(k_space).real})
        return corr_fg
    
class interpolation_comparison:
    def compare_interpolations():
        return 0