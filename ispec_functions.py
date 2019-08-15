#!/usr/bin/env python
#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import os
import sys
import numpy as np
import logging
import multiprocessing

################################################################################
#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/ted/iSpec/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


#--- Change LOG level ----------------------------------------------------------
#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))
################################################################################

def sendoff(y):
    logging.info("Processed spectrum saved as %s" % y)
    
def largest(a, b):
    if a > b:
        return a
    else:
        return b
    
def read_write_spectrum(spec_id):
    #--- Reading spectra -----------------------------------------------------------
    logging.info("Reading spectra")
    star_spectrum = ispec.read_spectrum(spec_id + "_datafile.txt")
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    ispec.write_spectrum(star_spectrum, spec_id + "_processed.fits")

def determine_radial_velocity_with_mask(spec_id):
    mu_cas_spectrum = ispec.read_spectrum(spec_id)
    #--- Radial Velocity determination with linelist mask --------------------------
    logging.info("Radial velocity determination with linelist mask...")
    # - Read atomic data
    mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Arcturus.372_926nm/mask.lst""
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Sun.372_926nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.A0.350_1095nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.F0.360_698nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.G2.375_679nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K0.378_679nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K5.378_680nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.M5.400_687nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/VALD.Sun.300_1100nm/mask.lst"
    ccf_mask = ispec.read_cross_correlation_mask(mask_file)

    models, ccf = ispec.cross_correlate_with_mask(mu_cas_spectrum, ccf_mask, \
                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                            velocity_step=1.0, mask_depth=0.01, \
                            fourier=False)

    # Number of models represent the number of components
    components = len(models)
    
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    
    # Applying correction
    logging.info("Applying radial velocity correction...")
    #rv = -96.40 # km/s
    mu_cas_spectrum = ispec.correct_velocity(mu_cas_spectrum, rv)
    ispec.write_spectrum(mu_cas_spectrum, spec_id)

def degrade_resolution(spec_id, initial_res, final_res):
    star_spectrum = ispec.read_spectrum(spec_id)
    #--- Resolution degradation ----------------------------------------------------
    logging.info("Degrading resolution...")
    from_resolution = initial_res
    to_resolution = final_res
    convolved_star_spectrum = ispec.convolve_spectrum(star_spectrum, to_resolution, \
                                                    from_resolution=from_resolution)
    ispec.write_spectrum(convolved_star_spectrum, spec_id)

def resample_spectrum(spec_id, min_wl, max_wl):
    star_spectrum = ispec.read_spectrum(spec_id)
    #--- Resampling  --------------------------------------------------------------
    logging.info("Resampling...")
    wavelengths = np.arange(min_wl, max_wl, 0.001)
    #resampled_star_spectrum = ispec.resample_spectrum(star_spectrum, wavelengths, method="bessel", zero_edges=True)
    resampled_star_spectrum = ispec.resample_spectrum(star_spectrum, wavelengths, method="linear", zero_edges=True)
    ispec.write_spectrum(resampled_star_spectrum, spec_id)

def normalize_whole_spectrum(spec_id, med_wave, max_wave):
    """
    Use the whole spectrum, strategy 'median+max'
    """
    star_spectrum = ispec.read_spectrum(spec_id)

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
    nknots = None # Automatic: 1 spline every 5 nm
    from_resolution = 40000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=med_wave
    max_wave_range=max_wave

    star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
    ispec.write_spectrum(normalized_star_spectrum, spec_id)

def find_fe_lines(spec_id):
    star_spectrum = ispec.read_spectrum(spec_id)
    
    #--- Continuum fit -------------------------------------------------------------
    logging.info("Fitting fixed continuum...")
    star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Find linemasks ------------------------------------------------------------
    logging.info("Finding line masks...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

    #telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

    # Read
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun
    #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
    #vel_telluric = 17.79  # km/s
    telluric_linelist = None
    vel_telluric = None

    resolution = 80000
    smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, resolution)
    min_depth = 0.05
    max_depth = 1.00
    star_linemasks = ispec.find_linemasks(star_spectrum, star_continuum_model, \
                            atomic_linelist=atomic_linelist, \
                            max_atomic_wave_diff = 0.005, \
                            telluric_linelist=telluric_linelist, \
                            vel_telluric=vel_telluric, \
                            minimum_depth=min_depth, maximum_depth=max_depth, \
                            smoothed_spectrum=star_spectrum, \
                            check_derivatives=False, \
                            discard_gaussian=False, discard_voigt=True, \
                            closest_match=False)
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (star_linemasks['wave_nm'] == 0)
    star_linemasks = star_linemasks[~rejected_by_atomic_line_not_found]

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (star_linemasks['ew'] == 0)
    star_linemasks = star_linemasks[~rejected_by_zero_ew]

    # Select only iron lines
    iron = star_linemasks['element'] == "Fe 1"
    iron = np.logical_or(iron, star_linemasks['element'] == "Fe 2")
    iron_star_linemasks = star_linemasks[iron]

    # Write iron regions with only masks limits and note:
    ispec.write_line_regions(iron_star_linemasks, "fe_linemasks.txt")

def find_linemasks(spec_id, species):
    code = "spectrum"
    star_spectrum = ispec.read_spectrum(spec_id)
    #--- Continuum fit -------------------------------------------------------------
    logging.info("Fitting fixed continuum...")
    star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")

    #--- Find linemasks ------------------------------------------------------------
    logging.info("Finding line masks...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

    #telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

    # Read
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun
    #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
    #vel_telluric = 17.79  # km/s
    telluric_linelist = None
    vel_telluric = None

    resolution = 80000
    smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, resolution)
    min_depth = 0.05
    max_depth = 1.00
    star_linemasks = ispec.find_linemasks(star_spectrum, star_continuum_model, \
                            atomic_linelist=atomic_linelist, \
                            max_atomic_wave_diff = 0.005, \
                            telluric_linelist=telluric_linelist, \
                            vel_telluric=vel_telluric, \
                            minimum_depth=min_depth, maximum_depth=max_depth, \
                            smoothed_spectrum=star_spectrum, \
                            check_derivatives=False, \
                            discard_gaussian=False, discard_voigt=True, \
                            closest_match=False)
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (star_linemasks['wave_nm'] == 0)
    star_linemasks = star_linemasks[~rejected_by_atomic_line_not_found]

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (star_linemasks['ew'] == 0)
    star_linemasks = star_linemasks[~rejected_by_zero_ew]

    # Select only lines of desired species
    elmt = star_linemasks['element'] == species
    elmt_star_linemasks = star_linemasks[elmt]
    
    file_id = species.replace(" ", "_") + "_linemasks.txt"

    # Write line regions with only masks limits and note:
    ispec.write_line_regions(elmt_star_linemasks, file_id)

def estimate_snr_from_flux(spec_id):
    star_spectrum = ispec.read_spectrum(spec_id)
    ## WARNING: To compare SNR estimation between different spectra, they should
    ##          be homogeneously sampled (consider a uniform re-sampling)
    #--- Estimate SNR from flux ----------------------------------------------------
    logging.info("Estimating SNR from fluxes...")
    num_points = 10
    estimated_snr = ispec.estimate_snr(star_spectrum['flux'], num_points=num_points)
    return estimated_snr

def estimate_snr_from_err(spec_id):
    star_spectrum = ispec.read_spectrum(spec_id)
    #--- Estimate SNR from errors --------------------------------------------------
    logging.info("Estimating SNR from errors...")
    efilter = star_spectrum['err'] > 0
    filtered_star_spectrum = star_spectrum[efilter]
    if len(filtered_star_spectrum) > 1:
        estimated_snr = np.median(filtered_star_spectrum['flux'] / filtered_star_spectrum['err'])
    else:
        # All the errors are set to zero and we cannot calculate SNR using them
        estimated_snr = 0
    return estimated_snr

def clean_telluric_regions(spec_id, clean_percent):
    star_spectrum = ispec.read_spectrum(spec_id)
    #--- Telluric velocity shift determination from spectrum --------------------------
    logging.info("Telluric velocity shift determination...")
    # - Telluric
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    models, ccf = ispec.cross_correlate_with_mask(star_spectrum, telluric_linelist, \
                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                            velocity_step=0.5, mask_depth=0.01, \
                            fourier = False,
                            only_one_peak = True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s

    #--- Clean regions that may be affected by tellurics ---------------------------
    logging.info("Cleaning tellurics...")

    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    # - Filter regions that may be affected by telluric lines
    #bv = 0.0
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], clean_percent)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(star_spectrum['waveobs'], \
                                telluric_linelist[dfilter], min_velocity=-bv+min_vel, \
                                max_velocity=-bv+max_vel)
    clean_star_spectrum = star_spectrum[~tfilter]
    ispec.write_spectrum(clean_star_spectrum, spec_id)

def determine_parameters(spec_id, lines_id):
    code = "spectrum"
    star_spectrum = ispec.read_spectrum(spec_id)
    # Use a fixed value because the spectrum is already normalized
    star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    #--- Model spectra ----------------------------------------------------------
    # Parameters
    initial_teff = 5771.0
    initial_logg = 4.44
    initial_MH = 0.00
    initial_vmic = 4.21 #ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
    initial_vmac = 4.21 #ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
    initial_vsini = 1.6
    initial_limb_darkening_coeff = 0.6
    initial_R = 40000
    initial_vrad = 0
    max_iterations = 6

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
    free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini"]

    # Free individual element abundance
    free_abundances = None
    linelist_free_loggf = None

    # Line regions
    line_regions = ispec.read_line_regions(lines_id)
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/spectrum_synth_turbospectrum_synth_sme_synth_moog_synth_synthe_synth_good_for_params_all_extended.txt")
    # Select only some lines to speed up the execution (in a real analysis it is better not to do this)
    #line_regions = line_regions[np.logical_or(line_regions['note'] == 'Ti 1', line_regions['note'] == 'Ti 2')]
    #line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)
    # Read segments if we have them or...
    #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    # ... or we can create the segments on the fly:
    segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
            linemasks=line_regions, \
            enhance_abundances=True, \
            use_errors = True, \
            vmic_from_empirical_relation = False, \
            vmac_from_empirical_relation = True, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code)
    ##--- Save results -------------------------------------------------------------
    logging.info("Saving results...")
    #dump_file = "example_results_synth_%s.dump" % (code)
    #ispec.save_results(dump_file, (params, errors, abundances_found, loggf_found, status, stats_linemasks))
    # If we need to restore the results from another script:
    #params, errors, abundances_found, loggf_found, status, stats_linemasks = ispec.restore_results(dump_file)

    #logging.info("Saving synthetic spectrum...")
    #synth_filename = "example_modeled_synth_%s.fits" % (code)
    #ispec.write_spectrum(modeled_synth_spectrum, synth_filename)
    
    parameters = "\n" + str(params)
    parameter_errors = "\n" + str(errors)
    file_id = spec_id + "_info.txt"
    
    with open(file_id, "w") as savefile:
        savefile.write(parameters)
        savefile.write(parameter_errors)

def determine_abundances(spec_id, species, params):
    multiprocessing.current_process().daemon=False
    code="spectrum"
    star_spectrum = ispec.read_spectrum(spec_id)
    
    # Use a fixed value because the spectrum is already normalized
    star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Model spectra ----------------------------------------------------------
    # Parameters
    initial_teff = params[1]
    initial_logg = params[5]
    initial_MH = params[4]
    initial_vmic = params[2]
    initial_vmac = params[0]
    initial_vsini = params[3]
    initial_limb_darkening_coeff = params[7]
    initial_R = params[6]
    initial_vrad = 0
    max_iterations = 6

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)



    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)


    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
    #free_params = ["vrad"]
    free_params = []

    # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

    element_name = species[:-2]
    free_abundances = ispec.create_free_abundances_structure([element_name], chemical_elements, solar_abundances)
    free_abundances['Abund'] += initial_MH # Scale to metallicity

    linelist_free_loggf = None
    
    linemask_id = species.replace(" ", "_") + "_linemasks.txt"

    # Line regions
    line_regions = ispec.read_line_regions(linemask_id)
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/limited_but_with_missing_elements_spectrum_synth_good_for_abundances_all_extended.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/limited_but_with_missing_elements_turobspectrum_synth_good_for_abundances_all_extended.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/limited_but_with_missing_elements_sme_synth_good_for_abundances_all_extended.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/limited_but_with_missing_elements_moog_synth_good_for_abundances_all_extended.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/limited_but_with_missing_elements_synthe_synth_good_for_abundances_all_extended.txt")
    # Select only the lines to get abundances from
    #line_regions = line_regions[np.logical_or(line_regions['note'] == element_name+' 1', line_regions['note'] == element_name+' 2')]
    #line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)

    # Read segments if we have them or...
    #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    # ... or we can create the segments on the fly:
    segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
            ispec.model_spectrum(star_spectrum, star_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
            linemasks=line_regions, \
            enhance_abundances=True, \
            use_errors = True, \
            vmic_from_empirical_relation = False, \
            vmac_from_empirical_relation = False, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code)

    abundance = str(abundances_found)
    results = abundance.split()
    del results[0:4]
    del results[1:3]
    del results[2:5]
    for i in range(0, len(results)):
        results[i] = float(results[i].replace(",",""))
    
    return results

teff_list = []
logg_list = []
vmic_list = []

def log_teff(result):
    teff_list.append(result)
    
def log_logg(result):
    logg_list.append(result)
    
def log_vmic(result):
    vmic_list.append(result)

def determine_error(spec_id, species, params, abundance, scatt_err):
    del teff_list[:]
    del logg_list[:]
    del vmic_list[:]
    
    plus_teff = list(params)
    plus_teff[1] += 100
    less_teff = list(params)
    less_teff[1] -= 100
    
    plus_logg = list(params)
    plus_logg[5] += 0.1
    less_logg = list(params)
    less_logg[5] -= 0.1
    
    plus_vmic = list(params)
    plus_vmic[2] += 0.5
    less_vmic = list(params)
    if params[2] > 0.5:
        less_vmic[2] -= 0.5
    
    #Enter a numerical argument here to restrict the number of cores used for multiprocessing.
    p = multiprocessing.Pool()
    
    p.apply_async(determine_abundances, [spec_id, species, plus_teff], callback = log_teff)
    p.apply_async(determine_abundances, [spec_id, species, less_teff], callback = log_teff)
    p.apply_async(determine_abundances, [spec_id, species, plus_logg], callback = log_logg)
    p.apply_async(determine_abundances, [spec_id, species, less_logg], callback = log_logg)
    p.apply_async(determine_abundances, [spec_id, species, plus_vmic], callback = log_vmic)
    p.apply_async(determine_abundances, [spec_id, species, less_vmic], callback = log_vmic)
    
    p.close()
    p.join()

    d_teff_p = abs(teff_list[0][0] - abundance)
    d_teff_m = abs(teff_list[1][0] - abundance)
    teff_error = largest(d_teff_p, d_teff_m)
    
    d_logg_p = abs(logg_list[0][0] - abundance)
    d_logg_m = abs(logg_list[1][0] - abundance)
    logg_error = largest(d_logg_p, d_logg_m)
    
    d_vmic_p = abs(vmic_list[0][0] - abundance)
    d_vmic_m = abs(vmic_list[1][0] - abundance)
    vmic_error = largest(d_vmic_p, d_vmic_m)
    
    error = (scatt_err**2 + teff_error**2 + logg_error**2 + vmic_error**2)**(0.5)
    return error
    
    