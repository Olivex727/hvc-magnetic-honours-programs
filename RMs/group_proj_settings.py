#! /Users/canders/miniconda3/envs/possum-mamba-env/bin/python python

"""
Routine settings to be read in by groups_basic_pipeline_output_fdf_analyser_ver10_collab_restructure.py
or later versions.
"""

##########################################################################################
###################################### SETTINGS ##########################################
##########################################################################################

basedir = '/Users/canders/Data/Groups/cutting_edge_data_2/'

include_CASDA_raw_RMs = True #Include the raw RMs downloaded from CASDA (which have no QA effective 9/01/2023)?
include_POSSUM_pipeline = False #Include the POSSUM pipeline data (which currently appear to have some QA issues too effective 9/01/2023)?
include_SPICE_RACS = False #Include the SPICE-RACS data?
do_not_chuck_sources = True #If true, will not use 'continue' in RM import loop to get rid of sources, but will instead apply a flag
mask_SR_cat_table = True #Mask SPICE-RACS table to get the 'good RM' dataset (speak to Alec)?
mask_results_table = True #Apply masks to the main RM results table, to e.g. apply mass or membership number cuts on group halos, etc.
do_NVSS_xmatching = False #X-match sources with NVSS RM catalogue?
save_data_out = False #Save out various numpy arrays, etc
save_final_data_out = True #Save out final table data
do_diagnostic_plots = True #Some additional plots, breaking down data into SBIDs to check for outliers/unusual behaviour, etc.
outlier_rejection_RM_correction = False #Before calculating the mean RM for gal foreground removal, should we reject outlier RMs?
add_control_spatial_offsets = False #To establish whether any effect is real, it helps to randomly offset the positions of the RMs w.r.t to the groups to check any result cannot come about randomly. 
factor_in_cosmology = True #The splashback radius cited by Tully assumed H_0==10 km/s Mpc. This switch recalculates R2T based on the nest mass and incorporating other possible H_0 values (set below).
do_search_dominant_galaxy_common_name = False #If True, query Simbad for the common name (NGC or Messier) of the brightest galaxy in the group
do_eliminate_RMs_in_bad_positions = True #If true, eliminate RMs that come from, or lie near, regions with v. bright and/or resolved sources
main_gal_centric_tully = True #Do we want to say the group is centred on the brightest galaxy, or on the mass-weighted mean position of the galaxies in the group?
add_noise_to_RMs_to_enforce_spatial_uniformity = True #The mosaic sensitivity drops off at the edges. Increase noise by the relevant factor for all RMs across the field to compensate?
limit_to_central_uniform_sensitivity_region = False #If true, will throw out all RMs not from central 2x2 degree regions of mosaic, where sensitivity is uniform. Since we select sources based on measured pol S/N, this is overkill --- no evidence of increased RM variance at mosaic edges for pol S/N>=X (see logbook 'Gal groups in POSSUM work #3' and tag #tag_no_lmc_goodness)
plot_RM_vs_sep_for_indiv_groups = False #Do RM vs. sep plots for individual groups?
do_hutschenreuter_foreground_RM_correction = False #Do the Hutschenreuter foreground correction?
