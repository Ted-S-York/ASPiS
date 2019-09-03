import ispec_functions as isp

#---------------------------USER INPUT & OPTIONS------------------------------#

#Format must match spectrum filename
spectrum_id = raw_input("Please define the spectrum to be processed: ")
info_id = spectrum_id + "_info.txt"

#Degree of telluric line preservation, 0 = complete removal, 100 = no removal
telluric_degree = 0

#Start- and endpoints for resolution degradation
start_resolution = 80000
final_resolution = 40000

#Median and maximum wave range for continuum normalisation, must be floats
med_wave = 3.0
max_wave = 4.0

#Minimum and maximum values for spectrum wavelength in nm, used in resampling
lambda_min = 424.76
lambda_max = 763.21

#----------------------------SPECTRUM PROCESSING------------------------------#

#Loads text file and converts to fits
isp.read_write_spectrum(spectrum_id)
spectrum_id += "_processed.fits"

#Tellluric and velocity corrections
isp.clean_telluric_regions(spectrum_id, telluric_degree)
isp.determine_radial_velocity_with_mask(spectrum_id)

#SNR estimations
snr_flux = isp.estimate_snr_from_flux(spectrum_id)
snr_err = isp.estimate_snr_from_err(spectrum_id)

#Create text file for additional spectrum info
with open(info_id, "w") as file:
    file.write('\nSNR from errors = %s' % snr_err)
    file.write('\nSNR from fluxes = %s' % snr_flux)

#Degrade resolution, normalise continuum and resample
#Number arguments are starting and ending resolution
isp.degrade_resolution(spectrum_id, start_resolution, final_resolution)
#Arguments are median and maximum wave range, must be floats 
isp.normalize_whole_spectrum(spectrum_id, med_wave, max_wave)
#Arguments are wavelength range of final spectrum
isp.resample_spectrum(spectrum_id, lambda_min, lambda_max)

isp.sendoff(spectrum_id)




