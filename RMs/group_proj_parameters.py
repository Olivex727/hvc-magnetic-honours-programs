#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
Functions to be read in by groups_basic_pipeline_output_fdf_analyser_ver10_collab_restructure.py
or later versions.
"""

from astropy.cosmology import Planck18 as assumed_comology

##########################################################################################
###################################### PARAMETERS ########################################
##########################################################################################

SN_cut = 10 #Pol signal:noise below which we turf RMs
n_member_in_group_cut = 3 #Min. number of group members identified by Tully to include in analysis
max_absRM_cut = 6000 #Maximum abs RM that we'll allow in our sample before simply cutting it
log10_Msol_cut_upper = 14 #Base 10 log of the most massive groups we'll consider ('legacy' sample spans ~12--16)
log10_Msol_cut_lower = 12.5 #Base 10 log of the least massive groups we'll consider ('legacy' sample spans ~12--16)
x_and_y_off_cut_groups = 3 #Degs; Group must be closer than this to tile centre along both the ra and decl directions to be counted in the analysis
x_and_y_off_cut_RMs = 2 #Only include RMs inside this ra and dec direction offset from mosaic centre in degs, so we dont include increased edge noise. See plots in 3rd group project notes (search "#money_plot 30 Jun 2023")
sbid_cut = 8000 #20000 SBID below which we will not consider
gal_lat_cut = 20 #Galactic 'b' [degs] plane-wards of which we'll not accept a tile 
NVSS_matching_radius_arcsec = 45 # separation in arcsec for considering that an NVSS source matches an ASKAP source
weight_scheme_RM_correction = 'None' #Weighting scheme for RM foreground correction. Currently if 'invvar', then inverse variance. If 'huber, uses robust estimator for mean. Otherwise, fefaults to equal weighting of all data points
num_RMs_to_average_for_foreground_correction = 40 # 20 #This is the number of RMs that we derive a foreground Galactic RM correction from for any given source
rm_correction_source_exclusion_radius_degs = 0.4 # 0.5 #Applicable to 'surrounding source' foreground RM correction method. For a given RM we want to foreground-correct, when we select nearby sources whose RMs we fit/average, then don't include any sources closer to this value. Important if the group environment produces spatially correlated RMs on scales smaller than this value, because then we might subtract some of the RM signal off. 
rm_outlier_limit_radius_degs = 0.6 #What radius INSIDE which do we consider sources for RM outlier exclusion stats? 0.56 ~= 1 square deg.
h = assumed_comology.H0.value/100 #lil H
n_groups_probed_inside_2TR_cut = 1 #If an RM penetrates greater than this number of group halos inside the splashback radius, then mask it out.
