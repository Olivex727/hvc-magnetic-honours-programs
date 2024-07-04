#!/usr/bin/env python3

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')

import numpy as np

from plotting import honours_plot as hplt
from collation import hvc_snapshot as snap, collation_tools as ct

from uncertainties import ufloat
from uncertainties import umath

from astropy import units as u
from astropy.units import astrophys as astru

import numpy as np
import scipy.stats as stats
from astropy.table import Table, vstack, hstack

import copy

import contextlib

def hvc_rms(i): return "../data_processed/hvc_rms/hvc_rms_index_"+str(i)

class hvc_looper:

    def add_survival_times_HVCs(collated_data, hvc_indicies=[]):
        return 0

    def add_magnetic_field_HVCs(collated_data, hvc_indicies=[], save_directory="../data_processed/", rm_load=True):
        return 0
    
    def uncertainty_subtract_HVCs(collated_data, hvc_indicies=[], load_directory="../data_processed/hvc_rms/", filter_significant=False, load_file="../data_processed/hvc_KS_tests/hvc_KS_average", save_file="../data_processed/results_pre"):
        master_hvcs = hvc_looper.load_HVCs(collated_data, hvc_indicies, load_directory)
        fwhm_table, _, _ = uncertainty_subtraction.subtract(master_hvcs)
        results = uncertainty_subtraction.uncertainty_readwrite(fwhm_table, filter_significant, load_file, save_file)
        return results
    
    def KStest_HVCs(collated_data, hvc_indicies=[], load_directory="../data_processed/hvc_rms/", save_file="../data_processed/hvc_KS_tests/hvc_KS", p_value=0.05, morph_type="average"):
        master_hvcs = hvc_looper.load_HVCs(collated_data, hvc_indicies, load_directory)
        results = KStest.KStest_HVCs(master_hvcs, p_value=p_value, morph_type="average")
        print("Converting")
        table_stat = Table(rows=results)
        if save_file:
            print("Saving data")
            ct.write_processed(table_stat, save_file+"_"+morph_type)
        print("Process Complete")
        return table_stat


    def add_magnetic_field_RMs(collated_data, hvc_indicies=[], save_directory="../data_processed/hvc_rms/", rm_load=True):
        rmbs = []
        print("=== CALCULATING HVC MAGNETIC FIELDS ===")
        print("Calculating HVC data")

        if not bool(hvc_indicies): hvc_indicies = list(range(len(collated_data["HVCs"])))
        l = len(hvc_indicies)

        for i in range(l):
            index = hvc_indicies[i]

            if rm_load: rm_load_file = hvc_rms(index)
            else: rm_load_file = ""

            with contextlib.redirect_stdout(None):
                rmb = magnetic_field_derivation.get_magnetic_field_points(snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False, rm_load_file=rm_load_file), index=index)

                if save_directory: ct.write_processed(rmb, save_directory+"hvc_rms_index_"+str(index)+"_with_B")

                rmbs.append(rmb)
            
            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
        
        print("Process complete")
        return rmbs
    
    def load_HVC_RMs(collated_data, hvc_indicies=[], directory="../data_processed/hvc_rms/", has_B=True):
        print("=== HVC RM LOADER ===")
        print("Taking HVC snapshots")

        rms = []

        if not bool(hvc_indicies): hvc_indicies = list(range(len(collated_data["HVCs"])))
        l = len(hvc_indicies)

        for i in range(l):
            index = hvc_indicies[i]

            with contextlib.redirect_stdout(None):
                rms.append(snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], rm_load_file=directory+"hvc_rms_index_"+str(index)+("_with_B" if has_B else ""))["RMs"])

            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
            
        print("Process complete")
        return rms
    
    def load_HVCs(collated_data, hvc_indicies=[], directory="../data_processed/hvc_rms/", has_B=True):
        print("=== HVC RM LOADER ===")
        print("Taking HVC snapshots")

        rms = []

        if not bool(hvc_indicies): hvc_indicies = list(range(len(collated_data["HVCs"])))
        l = len(hvc_indicies)

        for i in range(l):
            index = hvc_indicies[i]

            with contextlib.redirect_stdout(None):
                rms.append(snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], rm_load_file=directory+"hvc_rms_index_"+str(index)+("_with_B" if has_B else "")))

            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
            
        print("Process complete")
        return rms
    
    def save_HVC_RMs(collated_data, directory="../data_processed/hvc_rms/"):
        print("=== HVC RM SAVER ===")
        print("Taking HVC snapshots")
        l = len(collated_data["HVCs"])
        for index in range(l):
            with contextlib.redirect_stdout(None):
                snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], rm_save_file=directory+"hvc_rms_index_"+str(index))
            print(str(int((index+1)/l*100))+"% \r", sep="", end="", flush=True)
        print("Process complete")

    def plot_HVC_selection(hvc_indicies, collated_data, hvc_override=[], rm_load=True, scale=1, size=6, plot_cross_source=False, add_circles=False):
        snapshots = []
        print("=== GENERATING MULTIPLE HVC PLOTS ===")
        print("Calculating HVC data")
        l = len(hvc_indicies)
        override = bool(hvc_override)

        for i in range(len(hvc_indicies)):
            index = hvc_indicies[i]

            if rm_load: rm_load_file = hvc_rms(index)
            else: rm_load_file = ""

            with contextlib.redirect_stdout(None):
                if override:
                    snapshots.append(snap.take_snapshot(index, collated_data["RMs"], hvc_override, collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False, rm_load_file=rm_load_file))
                else:
                    snapshots.append(snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False, rm_load_file=rm_load_file))
            
            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
        
        print("Plotting HVC data")
        if plot_cross_source: hplt.plot_multiple_HVCs(snapshots, scale=scale, size=size, add_circles=add_circles)
        else: hplt.plot_multiple_HVCs_with_RM_sets(snapshots, scale=scale, size=size)

        print("Process complete")
        return snapshots
    
class magnetic_field_derivation:

    # DO NOT USE (YET)
    def get_magnetic_field_HVC(cropped_data, X=1):
        return magnetic_field_derivation.get_magnetic_field_points(cropped_data["RMs"])

    def get_magnetic_field_points(collated_data, X=1, index=0):
        print("=== CONVERTING RM TABLE ===")
        B_name = ["raw", "int"]
        B_list = [[], []]
        B_unc_list = [[], []]
        rm_table = collated_data["RMs"]
        hvc = collated_data["HVC"]

        l = len(rm_table)
        for index in range(len(rm_table)):
            B = magnetic_field_derivation.get_magnetic_field_point(rm_table[index], X, hvc["NH"], hvc["e_NH"])
            for i in range(len(B_name)):
                B_list[i].append(B[i].n)
                B_unc_list[i].append(B[i].s)
            print(str(int((index+1)/l*100))+"% \r", sep="", end="", flush=True)
            
        print("Process Complete")
        rm_table_new = copy.deepcopy(rm_table)

        for i in range(len(B_name)):
            rm_table_new.add_column(B_list[i], name="B_virtual ["+B_name[i]+"]")
            rm_table_new["B_virtual ["+B_name[i]+"]"].unit = u.G

            rm_table_new.add_column(B_unc_list[i], name="B_virtual_unc ["+B_name[i]+"]")
            rm_table_new["B_virtual_unc ["+B_name[i]+"]"].unit = u.G

        return rm_table_new

    # Returns magnetic field in gauss
    def get_magnetic_field_point(rm_point, X=1, H1=None, e_H1=None):
        return calculate.B_virt(
            H1,
            e_H1,
            rm_point["RM"],
            rm_point["RM_uncert"],
            rm_point["interpolation_raw"],
            rm_point["interpolation_unc"],
            X=X
            )

class calculate:

    # Returns in gauss
    def B_virt(H1, H1_err, rm, rm_unc, interp, interp_unc, intrinsic_unc=7, X=1):
        div = X * calculate.N_HI(H1, H1_err) * 1e6
        RMs = [
            3.8e18 * calculate.RM(rm, rm_unc, intrinsic_unc=intrinsic_unc)/div,
            3.8e18 * calculate.RM(rm, rm_unc, interp, interp_unc, intrinsic_unc=intrinsic_unc)/div
            ]

        return RMs
    
    # Returns in rad m-2
    def RM(rm, rm_unc, cor=0, cor_unc=0, intrinsic_unc=7):
        return ufloat(rm, rm_unc) - ufloat(cor, cor_unc) + calculate.intrinsic_RM_err(intrinsic_unc)
    
    def intrinsic_RM_err(rm_unc=7):
        return ufloat(0, rm_unc)

    def EM(H_alpha, H_alpha_err):
        return 2.75 * calculate.H_alpha(H_alpha, H_alpha_err) * calculate.T() ** 0.9
    
    # Returns temperature in 1e4 K
    def T():
        return ufloat(10, 2)

    # Returns H_alpha in Rayleighs
    def H_alpha(H_alpha, H_alpha_err):
        return ufloat(H_alpha, H_alpha_err)

    # Path length is in parsecs
    def path_length():
        return 1
    
    # HI Column Density
    def N_HI(H1, H1_err, individual_override=0):
        if individual_override:
            return ufloat(10 ** individual_override, 10 ** (0.5 * individual_override))
        else:
            return ufloat(H1, H1_err)
        
class KStest:

    def column_to_array(data):
        return data.data*1e6

    def split_RMs(RMs, centre, max_distance):
        mask = np.zeros(len(RMs), dtype=bool)
        for i in range(len(RMs)):
            rmi = RMs[i]
            mask[i] = rmi["ra_dec_obj"].separation(centre).value < max_distance

        RMs_inner = RMs[~mask]
        RMs_outer = RMs[mask]

        return RMs_inner, RMs_outer

    def morph_ring(hvc_snap, morph_type="average"):
        if morph_type == "average":
            return (hvc_snap["HVC"]["dx"] + hvc_snap["HVC"]["dy"])/2
        if morph_type == "minimum":
            return min(hvc_snap["HVC"]["dx"], hvc_snap["HVC"]["dy"])

    def make_cdfs(inner, outer):
        xs = np.linspace(-20, 20, 1000)

        inner_cdf = np.array(list(map(lambda x: KStest.edf(inner, x), xs)))
        outer_cdf = np.array(list(map(lambda x: KStest.edf(outer, x), xs)))

        return inner_cdf, outer_cdf, xs
    
    def edf(data, x):
        masked = data[data<x]
        return len(masked)/len(data)

    def KStest_single(snapshots, index = 0, show = False, dict_answer=True, p_value=0.05, morph_type="average", find_diff=True):
        hvc_snap = snapshots[index]
        inner_rms, outer_rms = KStest.split_RMs(hvc_snap["RMs"],hvc_snap["HVC"]["SkyCoord"], KStest.morph_ring(hvc_snap, morph_type=morph_type))
        inner = KStest.column_to_array(inner_rms["B_virtual [int]"])
        outer = KStest.column_to_array(outer_rms["B_virtual [int]"])

        print("Analysing HVC: " + hvc_snap["HVC"]["Name"])

        inner_cdf, outer_cdf, xs = KStest.make_cdfs(inner, outer)

        if show:
            hplt.plot_cdfs(xs, inner_cdf, xs, outer_cdf, show=False)

        ks_test = stats.ks_2samp(inner, outer)

        if find_diff:
            statx = ks_test.statistic_location
            statsgn = ks_test.statistic_sign
            statv = ks_test.statistic_sign

            y_inner = KStest.edf(inner, statx)

            y_outer = 0
            x_outer = 0
            for x in xs:
                y_outer = KStest.edf(outer, x)
                if y_outer >= y_inner:
                    x_outer = x
                    break

            diff = statx - x_outer

            if show:
                hplt.plot_cdf_lines(statx, statsgn, statv, y_inner, x_outer)

            if dict_answer:
                return {"Name":hvc_snap["HVC"]["Name"], "Statistic":ks_test.statistic, "p_value":ks_test.pvalue, "Statistic_x":ks_test.statistic_location, "Statistic_sgn":ks_test.statistic_sign, "Statistic_diff":diff, "Significant": ks_test.pvalue < p_value}

        if not dict_answer:
            return ks_test
        else:
            return {"Name":hvc_snap["HVC"]["Name"], "Statistic":ks_test.statistic, "p_value":ks_test.pvalue, "Statistic_x":ks_test.statistic_location, "Statistic_sgn":ks_test.statistic_sign, "Significant": ks_test.pvalue < p_value}
        
    def KStest_HVCs(snapshots, show=False, dict_answer=True, p_value=0.05, morph_type="average"):
        print("=== HVC KS TESTING ===")
        KSlist = []
        l = len(snapshots)
        print("Performing KS Tests")
        for i in range(len(snapshots)):
            with contextlib.redirect_stdout(None):
                KSlist.append(KStest.KStest_single(snapshots, index=i, show=show, dict_answer=dict_answer, p_value=p_value, morph_type=morph_type))
            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
        return KSlist
    
class uncertainty_subtraction:

    def subtract(snapshots):
        master_rm_inner, master_rm_outer, inners, outers = uncertainty_subtraction.get_stacked_sets(snapshots)
        return uncertainty_subtraction.uncertainty_subtract(inners, outers)
    
    def get_stacked_sets(snapshots):
        inners = []
        outers = []
        for hvc_snap in snapshots:
            inner_rms, outer_rms = KStest.split_RMs(hvc_snap["RMs"],hvc_snap["HVC"]["SkyCoord"], KStest.morph_ring(hvc_snap))
            inners.append(inner_rms)
            outers.append(outer_rms)

        master_rm_inner = vstack(inners)
        master_rm_outer = vstack(outers)

        return master_rm_inner, master_rm_outer, inners, outers
    
    def uncertainty_table(table_list):
        uncert = []
        for rms in table_list:
            m_list = rms["B_virtual_unc [int]"].data * 1e6
            o_list = rms["B_virtual [int]"].data * 1e6

            meas = np.mean(m_list)
            obsv = np.std(o_list)

            uncert.append({"Sigma [meas]":meas, "Sigma [obsv]":obsv, "Sigma [true]": np.sqrt(obsv**2 - meas**2)})

        uncert_table = Table(uncert)

        return uncert_table

    def uncertainty_subtract(inners, outers):
        inner_sigma = uncertainty_subtraction.uncertainty_table(inners)
        outer_sigma = uncertainty_subtraction.uncertainty_table(outers)

        sub = inner_sigma["Sigma [true]"]-outer_sigma["Sigma [true]"]
        fwhm = 2 * np.sqrt(2 * np.log(2)) * np.array(sub)

        fwhm_table = copy.deepcopy(inner_sigma)

        fwhm_table.add_column(sub, name="Sigma [diff]")
        fwhm_table.add_column(fwhm, name="FWHM")

        return fwhm_table, inner_sigma, outer_sigma
    
    def uncertainty_readwrite(uncert_table, filter_significant=False, load_file="../data_processed/hvc_KS_tests/hvc_KS_average", save_file="../data_processed/results_pre"):
        ks = ct.read_processed(load_file)
        hks = hstack([ks, uncert_table])

        if filter_significant:
            hks = hks[hks["Significant"]]
            hks = hks[~np.isnan(hks["Sigma [diff]"])]

        if save_file: ct.write_processed(hks, save_file)

        return hks
    
class postprocess_analysis:

    def generate_suvival_graph():
        return 0
    
    def get_survival_timescale():
        return 0