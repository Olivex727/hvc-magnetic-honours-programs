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

        if rms_pos.any(): plt.scatter(rms_pos[0], rms_pos[1], s=rms_pos[2] * scale, color=(1, 1, 1, 0), edgecolors='red', linewidth=2)
        if rms_neg.any(): plt.scatter(rms_neg[0], rms_neg[1], s=-1 * rms_neg[2] * scale, color=(1, 1, 1, 0), edgecolors='blue', linewidth=2)

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
        plt.imshow(image.data, origin='lower')
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
        plt.imshow(image_cropped, origin='lower')
        if not pix_c == [0, 0]: plt.plot((pix_c[0]-pix_up[0]),(pix_c[1]-pix_down[1]), 'rx')
        if index >= 0: plt.title("HVC index "+str(index))

        if show:
            plt.show()

    def plot_fits(image, show=True, color_map="Greys", scale_function=(lambda x: x)):
        if color_map: colormap = plt.colormaps[color_map]
        else: colormap = None

        image = plt.imshow(scale_function(image.data), cmap=colormap, origin='lower')

        if show:
            if color_map: plt.colorbar(image, label=r"[${\log}_{10}(N_{HI}/cm^{-1})$]")
            plt.xticks([])
            plt.yticks([])
            plt.show()
        
        return image
    
    def plot_fits_RM_overlay(rms, image, scale=1, show=True, index=-1, pixel_corners=[], title_prefix=""):
        if not index.isdigit() or index >= 0: plt.title(title_prefix+(" " if title_prefix else "")+"RM field for HVC "+str(index))
        
        honours_plot.plot_RMs(rms, scale=scale)
        image = honours_plot.plot_fits(image, False)

        if pixel_corners: plt.scatter((pixel_corners[2][0]-pixel_corners[1][0]),(pixel_corners[2][1]-pixel_corners[0][1]), marker='x', color=(0, 0, 0, 1))
        
        if show:
            plt.colorbar(image, label=r"[${\log}_{10}(N_{HI}/cm^{-1})$]")
            plt.show()

        return image
    
    def plot_colourspace_glat(rms, show=True, scale=0.1, color_map="coolwarm", x_col='RM', y_col="interpolation_raw", show_colorbar=True, xlabel=r"Faraday depth (Actual) [$rad m^{-2}$)]", ylabel=r"Faraday depth (Interpolation) [$rad m^{-2}$]", title="Comparison of RMs", return_color=False):
        #plt.rcParams.update({'font.size': 13})

        b_list = rms['ra_dec_obj'].galactic.b

        colormap = plt.colormaps[color_map]

        scatter = plt.scatter(rms[x_col], rms[y_col], marker='s', s=scale, c=b_list, cmap=colormap)#colors)#, xerr=filtered['faraday_depth_err_radmm'], yerr=filtered['interpolation_err'], ecolor = "black")
        if show_colorbar: plt.colorbar(scatter, label=r"Galactic Latitude [$deg$]")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if show:
            plt.show()

        if return_color:
            return scatter

    def plot_heatmap_single(RMs, x_col='RM', y_col="interpolation_raw", xlabel="", ylabel="", title="", bins=20, rng=[[-80, 80], [-180, 180]], show=True):
        hist = plt.hist2d(RMs[x_col], RMs[y_col], bins=bins, range=rng)

        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)
        if title: plt.title(title)

        if show:
            plt.show()
        
        return hist
    
    def plot_H_alpha(Ha, bins=100, rng=(0,5), show=True):
        plt.hist(Ha, bins, rng, color=[0.8, 0.1, 0.8, 0.6])
        plt.title("RM sample associated H-alpha fluxes")
        plt.xlabel(r"Flux [$R$]")
        plt.ylabel("Counts")

        if show:
            plt.show()
    
    def plot_HI(H1, bins=100, rng=(0,5), show=True):
        plt.hist(H1, bins, rng, color=[0.1, 0.8, 0.8, 0.6])
        plt.title("RM sample associated HI fluxes")
        plt.xlabel(r"Column Density [$dex(cm^{-2})$]")
        plt.ylabel("Counts")

        if show:
            plt.show()

    def plot_RM_histogram_single(set_1, set_2, set_1_name="", set_2_name="", title="", bounds=(-100,100), bins=100, show=True, ylabel="", xlabel=r"Faraday depth [$rad m^{-2}$]", maximum=0, second=False, legend_size=14):

        if not second:
            y, x, _ = plt.hist(set_1-set_2, bins, bounds, label="Residuals", color=[0.1, 0.1, 0.8, 0.4], density=True)
            if not maximum: maximum = y.max()*1.05
        elif set_2_name:
            y, x, _ = plt.hist(set_2, bins, bounds, label=set_2_name, color=[0.1, 0.8, 0.1, 0.4], density=True)
            if not maximum: maximum = y.max()*1.05

        plt.ylim(0, maximum)
        if set_1_name:
            plt.hist(set_1, bins, bounds, label=set_1_name, color=[0.8, 0.1, 0.1, 0.4], density=True)
            plt.ylim(0, maximum)

        if set_2_name and not second:
            plt.hist(set_2, bins, bounds, label=set_2_name, color=[0.1, 0.8, 0.1, 0.4], density=True)
            plt.ylim(0, maximum)
        else:
            plt.hist(set_1-set_2, bins, bounds, label="Residuals", color=[0.1, 0.1, 0.8, 0.4], density=True)
            plt.ylim(0, maximum)

        plt.legend(fontsize=legend_size, loc='upper right')
        if ylabel: plt.ylabel(ylabel)
        if xlabel: plt.xlabel(xlabel)
        if title: plt.title(title)

        if show:
            plt.show()
    
        return maximum
    
    def plot_3hist(set_1, set_2, set_3, set_1_name="", set_2_name="", set_3_name="", title="", bounds=(-100,100), bins=100, show=True, ylabel="", xlabel=r"Faraday depth [$rad m^{-2}$]", legend_size=14):
        if set_1_name: plt.hist(set_1, bins, bounds, label=set_1_name, color=[0.8, 0.1, 0.1, 0.4])
        if set_2_name: plt.hist(set_2, bins, bounds, label=set_2_name, color=[0.1, 0.8, 0.1, 0.4])
        if set_3_name: plt.hist(set_3, bins, bounds, label=set_3_name, color=[0.1, 0.1, 0.8, 0.4])
        plt.legend(fontsize=legend_size, loc='upper right')
        if ylabel: plt.ylabel(ylabel)
        if xlabel: plt.xlabel(xlabel)
        if title: plt.title(title)

        if show:
            plt.show()

    # NB: Plots must come in multiples of 3
    def plot_multiple_HVCs(snapshots, scale=1, size=6, show=True, add_circles=False, average=False):
        ny_plots = int(len(snapshots) / 3)

        plt.figure(figsize=(size*3, ny_plots*size))

        plt.rcParams.update({'font.size': (1+size)*2})

        plt.tight_layout(pad=0)

        for s in range(len(snapshots)):
            snapshot = snapshots[s]

            if average:
                interpolation = np.mean(snapshot["RMs"]["RM"])
            else:
                interpolation = snapshot["RMs"]["interpolation_raw"]
            print(interpolation)
            rm_overlay = np.array([
                snapshot["RMs"]["pixel location x"],
                snapshot["RMs"]["pixel location y"],
                snapshot["RMs"]["RM"] - interpolation
                ])
            plt.subplot(ny_plots, 3, s+1)
            plt.axis([0, snapshot['HI'].shape[0]-2, 0, snapshot['HI'].shape[1]-2])
            plt.tight_layout(w_pad=0, h_pad=1)
            plt.margins(tight=True)
            honours_plot.plot_fits_RM_overlay(rm_overlay, snapshot["HI"], show=False, index=snapshot["HVC"]["Name"], pixel_corners=snapshot["HI_pixel_corners"], scale=scale)
            plt.xticks([])
            plt.yticks([])
            if add_circles:
                xlim = plt.xlim()
                ylim = plt.ylim()

                maximum = max(snapshot["HVC"]["dx"], snapshot["HVC"]["dy"])
                average = (snapshot["HVC"]["dx"]+snapshot["HVC"]["dy"])/2

                circle = plt.Circle((sum(xlim)/2, sum(ylim)/2), (average/maximum)*sum(xlim)/4, color='black', fill=False)
                plt.gca().add_patch(circle)
        
        if show:
            plt.show()

    def plot_multiple_HVCs_with_RM_sets(snapshots, scale=1, size=6, show=True):
        image = None

        ny_plots = len(snapshots)

        fig = plt.figure(figsize=((4+(size*3), ny_plots*size)))

        plt.rcParams.update({'font.size': (size)*2})

        plt.tight_layout(pad=0)

        keywords = ["RM", "interpolation_raw", "interpolation_cor"]
        titleword = ["Uncorrected", "Intp. corrected", "Fourier corrected"]

        for s in range(len(snapshots) * 3):
            snapshot = snapshots[int(s / 3)]
            rm_overlay = np.array([
                snapshot["RMs"]["pixel location x"],
                snapshot["RMs"]["pixel location y"],
                snapshot["RMs"]["RM"] - (snapshot["RMs"][keywords[s % 3]] if s % 3 != 0 else 0)
                ])
            plt.subplot(ny_plots, 3, s+1)
            plt.axis([0, snapshot['HI'].shape[0]-2, 0, snapshot['HI'].shape[1]-2])
            plt.tight_layout(w_pad=0, h_pad=1)
            plt.margins(tight=True)
            image = honours_plot.plot_fits_RM_overlay(rm_overlay, snapshot["HI"], show=False, index=snapshot["HVC"]["Name"], pixel_corners=snapshot["HI_pixel_corners"], scale=scale, title_prefix=titleword[s % 3])

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
        fig.colorbar(image, cax=cbar_ax, label=r"HI Intensity [${\log}_{10}(N_{HI}/cm^{-1})$]")

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

    def plot_cdfs(data1, cdf1, data2, cdf2, show=False, xlims=(-20,20)):
        plt.plot(data1, cdf1, label="Inner")
        plt.plot(data2, cdf2, label="Outer")
        plt.xlabel(r"Faraday Depth [$rad m^{-2}$]")
        plt.ylabel(r"Cumulative Proportion")
        plt.xlim(xlims)
        plt.ylim(0,1)
        plt.legend()
        if show:
            plt.show()

    def plot_cdf_lines(statx, statsgn, statv, y_inner, x_outer, limits):
        plt.axvline(x=statx/limits[1], ymin=y_inner if statsgn < 0 else y_inner-statv, ymax=y_inner+statv if statsgn < 0 else y_inner, c="red", label="Statisitc")

        plt.axhline(y=y_inner, xmin=(limits[0]+(x_outer if statsgn < 0 else statx))/(2*limits[0]), xmax=(limits[0]+(statx if statsgn < 0 else x_outer))/(2*limits[0]), c="black", label="Difference")

        plt.legend()

    def uncertainty_boxplot(master_rm, single=False, limit=20):
        rms = master_rm["B_virtual_unc [int]"].data * 1e6
        B_unc = rms[rms < 20]
        B_true_unc = np.std(master_rm["B_virtual [int]"].data * 1e6)

        if single:
            plt.figure(figsize=(6,1))
            plt.boxplot(B_unc, vert=False, showmeans=True, widths=0.6, sym="x")
            plt.axvline(B_true_unc, c='r', linestyle='--')
            plt.yticks([])
            plt.xlabel(r"Faraday Depth [$10^7 rad m^{-2}$]")
        else:
            plt.boxplot(B_unc, vert=False, showmeans=True, widths=0.6, sym="x")
            plt.axvline(B_true_unc, c='r', linestyle='--')
            plt.yticks([])

    def uncertainty_boxplots(master_rm_inner, master_rm_outer):
        fig = plt.figure(figsize=(6,2.5))
    
        fig.supxlabel(r"Magnetic Field Uncertainties [$\mu G$]")

        plt.tight_layout()

        plt.subplot(2, 1, 1)
        honours_plot.uncertainty_boxplot(master_rm_inner)
        plt.subplot(2, 1, 2)
        honours_plot.uncertainty_boxplot(master_rm_outer)
        plt.subplots_adjust(hspace=0.5, bottom=0.2)