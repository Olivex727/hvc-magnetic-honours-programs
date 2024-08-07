# === SETUP === #

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/project')

from collation import collator, hvc_snapshot as snap, collation_tools as ct
from synthesis import hvc_looper as hvcl, magnetic_field_derivation as mfd, KStest, uncertainty_subtraction as us

collated_data = collator.data_whole_sky(False, load_data=["../data_processed/proc_rms","../data_processed/proc_hvcs"], h1_img="../data_catalog/hi4pi-hvc-nhi-car.fits", override_RMs=True)

master_hvcs = hvcl.load_HVCs(collated_data, hvc_indicies=[0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 13, 24])

results = ct.read_processed("../data_processed/results_weighted_RM")

# === KS-EDF UNCERTAINTIES === #

import numpy as np

master_rm_inner, master_rm_outer, inners, outers = us.get_stacked_sets(master_hvcs)
inner_sigma = us.uncertainty_table(inners)
outer_sigma = us.uncertainty_table(outers)

#uncertainty_KS = np.sqrt(inner_sigma["Sigma [meas]"]**2 + inner_sigma["Sigma [obsv]"]**2 + outer_sigma["Sigma [meas]"]**2 + outer_sigma["Sigma [obsv]"]**2)
uncertainty_KS = np.sqrt(inner_sigma["Sigma [meas]"]**2 + outer_sigma["Sigma [meas]"]**2)

results.add_column(uncertainty_KS, name="KS unc")

# === BOOTSTRAP FUNCTIONS === #

from astropy.table import Table
import copy

def bootstrap_selection(sample):
    sample = copy.deepcopy(sample)
    sample.remove_column("ra_dec_obj")
    rand = np.round((len(sample)-1) * np.random.rand(len(sample))).astype(int)
    t = copy.deepcopy(sample)
    t.remove_rows(list(range(len(sample))))

    for index in rand:
        t.add_row(list(sample[index]))
    
    return t

def bootstrap_sample_creation(sample, console_out=""):
    samples = []
    l = len(sample)
    for i in range(len(sample)):
        samples.append(bootstrap_selection(sample))
        print(console_out+"Creating samples: "+str(int((i+1)/l*100))+"% \r", sep="", end="", flush=True)
    return samples

def bootstrap_evaluation(samples, callback):
    sample_out = []
    l = len(samples)

    for i in range(len(samples)):
        sout = "Evaluating samples: "+str(int((i+1)/l*100))+"% "
        sample = samples[i]
        bootstrapped = bootstrap_sample_creation(sample, sout+"")
        response = np.array(list(map(callback, bootstrapped)))
        sample_out.append(response)
        print(sout+"\r", sep="", end="", flush=True)

    return np.array(sample_out)

def uncertainty_calculate(rms):
    m_list = rms["RM_uncert"].data
    o_list = rms["RM"].data

    meas = np.mean(m_list)
    obsv = np.std(o_list)

    return np.sqrt(obsv**2 - meas**2)

import warnings

def bootstrap_uncertaintes(master_hvcs):
    with warnings.catch_warnings(action="ignore"):
        print("=== BOOTSTRAPPING UNCERTAINTIES ===")
        print("Getting stacked sets")
        master_rm_inner, master_rm_outer, inners, outers = us.get_stacked_sets(master_hvcs)
        print("Calculating inner magnetic field uncertainties")
        inner_mag = bootstrap_evaluation(inners, uncertainty_calculate)
        inner_unc = []
        for l in inner_mag:
            inner_unc.append(np.nanstd(l))
        inner_unc = np.array(inner_unc)
        print() 
        print("Calculating outer magnetic field uncertainties")
        outer_mag = bootstrap_evaluation(outers, uncertainty_calculate)
        outer_unc = []
        for l in outer_mag:
            outer_unc.append(np.nanstd(l))
        outer_unc = np.array(outer_unc)
        print()
        print("Calculating final uncertainty")
        final = np.sqrt(inner_unc ** 2 + outer_unc ** 2)
        return final

field_set = bootstrap_uncertaintes(master_hvcs)

results.add_column(field_set, name="Sigma unc")

ct.write_processed(results, "../data_processed/results_raw_RM")