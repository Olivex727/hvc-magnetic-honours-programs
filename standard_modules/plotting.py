#!/usr/bin/env python3

from PIL import Image
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class honours_plot:

    def plot_RMs(rms, scale=1, show=False):
        rms_pos = np.array(list(filter(lambda val: val[2] > 0, rms.transpose()))).transpose()
        rms_neg = np.array(list(filter(lambda val: val[2] < 0, rms.transpose()))).transpose()

        plt.scatter(rms_pos[0] * scale, rms_pos[1] * scale, s=rms_pos[2], color=(1, 1, 1, 0), edgecolors='red', linewidth=2)
        plt.scatter(rms_neg[0] * scale, rms_neg[1] * scale, s=-1 * rms_neg[2], color=(1, 1, 1, 0), edgecolors='blue', linewidth=2)

        if show:
            plt.show()
    
    def plot_RMs_overlay(x_pixels, y_pixels, rms, h1, corners=[], index=-1):
        rm_overlay = np.array([
            x_pixels,
            y_pixels,
            rms
        ])
        honours_plot.plot_fits_RM_overlay(rm_overlay, h1, show=True, index=index, pixel_corners=corners)

    def plot_image_crop(image, image_cropped, pix_up, pix_down, pix_c=[0, 0], index=-1, show=True):
        if index >= 0: plt.title("Location of HVC index "+str(index))
        plt.xticks([])
        plt.yticks([])
        plt.imshow(image.data)
        if not pix_c == [0, 0]: plt.plot(*pix_c, 'bx', ms=1)
        plt.gca().add_patch(Rectangle(
            (pix_up[0], pix_down[1]),
            np.abs(pix_up[0]-pix_down[0]),
            np.abs(pix_up[1]-pix_down[1]),
            linewidth=1,edgecolor='r',facecolor='none'
            ))
        if show:
            plt.show()

        plt.xticks([])
        plt.yticks([])
        plt.imshow(image_cropped)
        if not pix_c == [0, 0]: plt.plot((pix_c[0]-pix_up[0]),(pix_c[1]-pix_down[1]), 'rx')
        if index >= 0: plt.title("HVC index "+str(index))

        if show:
            plt.show()

    def plot_fits(image, show=True):
        plt.xticks([])
        plt.yticks([])
        plt.imshow(image.data)

        if show:
            plt.show()
    
    def plot_fits_RM_overlay(rms, image, scale=1, show=True, index=-1, pixel_corners=[]):
        if index >= 0: plt.title("RM field for HVC "+str(index))
        
        honours_plot.plot_RMs(rms, scale=scale)
        honours_plot.plot_fits(image, False)

        if pixel_corners: plt.scatter((pixel_corners[2][0]-pixel_corners[1][0]),(pixel_corners[2][1]-pixel_corners[0][1]), marker='x', color=(0, 0, 0, 1))
        
        if show:
            plt.show()
    
    def plot_colourspace_glat(rms, show=True, scale=0.1, color_map="gist_rainbow", x_col='RM', y_col="interpolation_raw", show_colorbar=True, xlabel=r"Faraday depth (Actual) [$rad m^{-2}$)]", ylabel=r"Faraday depth (Interpolation) [$rad m^{-2}$]", title="Comparison of RMs"):
        b_list = rms['ra_dec_obj'].galactic.b

        colormap = plt.colormaps[color_map]

        scatter = plt.scatter(rms[x_col], rms[y_col], marker='s', s=scale, c=b_list, cmap=colormap)#colors)#, xerr=filtered['faraday_depth_err_radmm'], yerr=filtered['interpolation_err'], ecolor = "black")
        if show_colorbar: plt.colorbar(scatter, label=r"Galactic Latitude [$deg$]")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if show:
            plt.show()

    def plot_heatmap_single(x_col, y_col, xlabel="", ylabel="", title="", bins=20, rng=[[-80, 80], [-180, 180]], show=True):
        plt.hist2d(x_col, y_col, bins=bins, range=rng)

        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)
        if title: plt.title(title)

        if show:
            plt.show()
    
    def plot_H_alpha(Ha, bins=100, rng=(0,5), show=True):
        plt.hist(Ha, bins, rng, color=[0.8, 0.1, 0.8, 0.6])
        plt.title("RM gird point H-alpha fluxes")
        plt.xlabel(r"Flux [$R$]")
        plt.ylabel("Counts")

        if show:
            plt.show()

    def plot_RM_histogram_single(set_1, set_2, set_1_name="", set_2_name="", title="", bounds=(-100,100), bins=100, show=True, ylabel="", xlabel=r"Faraday depth [$rad m^{-2}$]"):
        if set_1_name: plt.hist(set_1, bins, bounds, label=set_1_name, color=[0.8, 0.1, 0.1, 0.4])
        if set_2_name: plt.hist(set_2, bins, bounds, label=set_2_name, color=[0.1, 0.8, 0.1, 0.4])
        plt.hist(set_1-set_2, bins, bounds, label="Residuals", color=[0.1, 0.1, 0.8, 0.4])
        plt.legend()
        if ylabel: plt.ylabel(ylabel)
        if xlabel: plt.xlabel(xlabel)
        if title: plt.title(title)

        if show:
            plt.show()
    
    def plot_3hist(set_1, set_2, set_3, set_1_name="", set_2_name="", set_3_name="", title="", bounds=(-100,100), bins=100, show=True, ylabel="", xlabel=r"Faraday depth [$rad m^{-2}$]"):
        if set_1_name: plt.hist(set_1, bins, bounds, label=set_1_name, color=[0.8, 0.1, 0.1, 0.4])
        if set_2_name: plt.hist(set_2, bins, bounds, label=set_2_name, color=[0.1, 0.8, 0.1, 0.4])
        if set_3_name: plt.hist(set_3, bins, bounds, label=set_3_name, color=[0.1, 0.1, 0.8, 0.4])
        plt.legend()
        if ylabel: plt.ylabel(ylabel)
        if xlabel: plt.xlabel(xlabel)
        if title: plt.title(title)

        if show:
            plt.show()

    # NB: Plots must come in multiples of 3
    def plot_multiple_HVCs(snapshots, size=6, show=True):
        ny_plots = int(len(snapshots) / 3)

        plt.figure(figsize=(size*3, ny_plots*size))

        plt.rcParams.update({'font.size': 14})

        plt.tight_layout(pad=0)

        for s in range(len(snapshots)):
            snapshot = snapshots[s]
            rm_overlay = np.array([
                snapshot["RMs"]["pixel location x"],
                snapshot["RMs"]["pixel location y"],
                snapshot["RMs"]["RM"]
                ])
            plt.subplot(ny_plots, 3, s+1)
            plt.axis([0, snapshot['HI'].shape[0]-2, 0, snapshot['HI'].shape[1]-2])
            plt.tight_layout(w_pad=0, h_pad=1)
            plt.margins(tight=True)
            honours_plot.plot_fits_RM_overlay(rm_overlay, snapshot["HI"], show=False, index=snapshot["index"], pixel_corners=snapshot["HI_pixel_corners"])
        
        if show:
            plt.show()

    def plot_collated_RMs_multi(RMs, show=True):
        def repeated_setup(n):
            plt.subplot(3, 3, n)
            plt.tight_layout()
            plt.margins(tight=True)

        plt.figure(figsize=(13, 12))

        #plt.rcParams.update({'font.size': 14})

        repeated_setup(1)
        honours_plot.plot_colourspace_glat(RMs, show=False, show_colorbar=False, xlabel="Actual RMs", ylabel=r"Faraday depth [$rad m^{-2}$]" + "\n" + "Interpolated RMs", title="")

        repeated_setup(2)
        honours_plot.plot_colourspace_glat(RMs, show=False, show_colorbar=False, xlabel=r"Actual RMs" + "\n" + r"Faraday depth [$rad m^{-2}$]", ylabel=r"Bandpassed RMs", title="RM sample colour maps")

        repeated_setup(3)
        honours_plot.plot_colourspace_glat(RMs, show=False, show_colorbar=True, xlabel=r"Interpolated RMs", ylabel=r"Bandpassed RMs", title="")

        repeated_setup(4)
        honours_plot.plot_heatmap_single(RMs["RM"], RMs['ra_dec_obj'].galactic.b, r"Actual RMs", r"Galactic Longnitude [$deg$]", bins=50, show=False)

        repeated_setup(5)
        honours_plot.plot_heatmap_single(RMs["interpolation_raw"], RMs['ra_dec_obj'].galactic.b, r"Interpolated RMs" + "\n" + r"[$rad m^{-2}$]", "", bins=50, show=False)

        repeated_setup(6)
        honours_plot.plot_heatmap_single(RMs["interpolation_cor"], RMs['ra_dec_obj'].galactic.b, r"Bandpassed RMs", "", bins=50, show=False, title="Heatmap of RM-latitude")

        repeated_setup(7)
        honours_plot.plot_RM_histogram_single(RMs["RM"], RMs["interpolation_raw"], "Actual RMs", "Interpolated RMs", bounds=(-50,50), ylabel="Counts", xlabel="")

        repeated_setup(8)
        honours_plot.plot_RM_histogram_single(RMs["RM"], RMs["interpolation_cor"], "Actual RMs", "Bandpassed RMs", bounds=(-50,50), title="Histograms of RM comparison residuals")

        repeated_setup(9)
        honours_plot.plot_RM_histogram_single(RMs["interpolation_raw"], RMs["interpolation_cor"], "Interpolated RMs", "Bandpassed RMs", bounds=(-50,50), xlabel="")

        if show:
            plt.show()