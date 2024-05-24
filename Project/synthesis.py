#!/usr/bin/env python3

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')

import numpy as np

from plotting import honours_plot as hplt
from collation import hvc_snapshot as snap

import contextlib

class hvc_looper:

    def add_magnetic_field_HVCs(collated_data, hvc_indicies=[]):
        return 0
    
    def add_survival_times_HVCs(collated_data, hvc_indicies=[]):
        return 0

    def plot_HVC_selection(hvc_indicies, collated_data, hvc_override=[]):
        snapshots = []
        print("=== GENERATING MULTIPLE HVC PLOTS ===")
        print("Calculating HVC data")
        l = len(hvc_indicies)
        override = bool(hvc_override)

        for i in range(len(hvc_indicies)):
            index = hvc_indicies[i]
            with contextlib.redirect_stdout(None):
                if override:
                    snapshots.append(snap.take_snapshot(index, collated_data["RMs"], hvc_override, collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False))
                else:
                    snapshots.append(snapshot = snap.take_snapshot(index, collated_data["RMs"], collated_data["HVCs"], collated_data["HI"], collated_data["H-alpha"], collated_data["interpolation"], plot=False))
            
            print(str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
        
        print("Plotting HVC data")
        hplt.plot_multiple_HVCs(snapshots)

        return snapshots
    
class magnetic_field_derivation:

    def get_magnetic_field_HVC(cropped_data):
        return 0

    def get_magnetic_field_points(RM_table, path_length):
        return 0

    def get_magnetic_field_point(RM_point, path_length):
        return 0
    
    def calculate_path_length():
        return 0
    
class postprocess_analysis:

    def generate_suvival_graph():
        return 0
    
    def get_survival_timescale():
        return 0