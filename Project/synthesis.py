#!/usr/bin/env python3

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')

import numpy as np

from plotting import honours_plot as hplt
from collation import hvc_snapshot as snap

from uncertainties import ufloat
from uncertainties import umath

from astropy import units as u
from astropy.units import astrophys as astru

import copy

import contextlib

def hvc_rms(i): return "../data_processed/hvc_rms/hvc_rms_index_"+str(i)

class hvc_looper:

    def add_magnetic_field_HVCs(collated_data, hvc_indicies=[]):
        return 0
    
    def add_survival_times_HVCs(collated_data, hvc_indicies=[]):
        return 0
    
    def save_HVC_RMs(collated_data, directory="../data_processed/hvc_rms/"):
        print("=== HVC RM SAVER ===")
        print("Taking HVC snapshots")
        l = len(collated_data["HVCs"])
        for index in range(l):
            with contextlib.redirect_stdout(None):
                snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], rm_save_file=directory+"hvc_rms_index_"+str(index))
            print(str(int((index+1)/l*100))+"% \r", sep="", end="", flush=True)
        print("Process complete")


    def plot_HVC_selection(hvc_indicies, collated_data, hvc_override=[], rm_load=True):
        snapshots = []
        print("=== GENERATING MULTIPLE HVC PLOTS ===")
        print("Calculating HVC data")
        l = len(hvc_indicies)
        override = bool(hvc_override)

        for i in range(len(hvc_indicies)):
            index = hvc_indicies[i]

            if rm_load: rm_load_file = hvc_rms(i)
            else: rm_load_file = ""

            with contextlib.redirect_stdout(None):
                if override:
                    snapshots.append(snap.take_snapshot(index, collated_data["RMs"], hvc_override, collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False, rm_load_file=rm_load_file))
                else:
                    snapshots.append(snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False, rm_load_file=rm_load_file))
            
            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
        
        print("Plotting HVC data")
        hplt.plot_multiple_HVCs(snapshots)

        #return snapshots
    
class magnetic_field_derivation:

    def get_magnetic_field_HVC(cropped_data):
        return magnetic_field_derivation.get_magnetic_field_points(cropped_data["RMs"])

    def get_magnetic_field_points(rm_table, path_length=1):
        print("=== CONVERTING RM TABLE ===")
        B_name = ["raw", "int", "cor"]
        B_list = [[], [], []]
        B_unc_list = [[], [], []]
        l = len(rm_table)
        for index in range(len(rm_table)):
            B = magnetic_field_derivation.get_magnetic_field_point(rm_table[index], path_length)
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
    def get_magnetic_field_point(rm_point, path_length=1):
        return calculate.B_virt(
            rm_point["H-alpha flux"],
            rm_point["H-alpha flux [Error]"],
            rm_point["RM"],
            rm_point["RM_uncert"],
            rm_point["interpolation_raw"],
            rm_point["interpolation_unc"],
            rm_point["interpolation_cor"]
            )

class calculate:

    # Returns in gauss
    def B_virt(H_alpha, H_alpha_err, rm, rm_unc, interp, interp_unc, interp_cor, intrinsic_unc=7):
        div = 0.81 * (calculate.path_length() * calculate.EM(H_alpha, H_alpha_err)) ** 0.5
        RMs = [
            1e-6 * calculate.RM(rm, rm_unc, intrinsic_unc=intrinsic_unc)/div,
            1e-6 * calculate.RM(rm, rm_unc, interp, interp_unc, intrinsic_unc=intrinsic_unc)/div,
            1e-6 * calculate.RM(rm, rm_unc, interp_cor, interp_unc, intrinsic_unc=intrinsic_unc)/div
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
    
class postprocess_analysis:

    def generate_suvival_graph():
        return 0
    
    def get_survival_timescale():
        return 0