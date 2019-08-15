import ispec_functions as isp

#Change to match desired target, format must match spectrum filename
spectrum_id = raw_input("Please define the spectrum to be processed: ")
info_id = spectrum_id + "_info.txt"
title = spectrum_id

#Loads text file and converts to fits
isp.read_write_spectrum(spectrum_id)
spectrum_id += "_processed.fits"

#Tellluric and velocity corrections
#Numerical argument is degree of removal of Tellurics, 0 = complete removal, 100 = no removal
isp.clean_telluric_regions(spectrum_id, 0)
isp.determine_radial_velocity_with_mask(spectrum_id)

#SNR estimations
snr_err = isp.estimate_snr_from_flux(spectrum_id)
snr_flux = isp.estimate_snr_from_err(spectrum_id)

#Create text file for additional spectrum info
with open(info_id, "w") as file:
    file.write(title)
    file.write('\nSNR from errors = %s' % snr_err)
    file.write('\nSNR from fluxes = %s' % snr_flux)

#Degrade resolution, normalise continuum and resample
#Number arguments are starting and ending resolution
isp.degrade_resolution(spectrum_id, 80000, 40000)
#Arguments are median and maximum wave range, must be floats 
isp.normalize_whole_spectrum(spectrum_id, 3.0, 4.0)
#Arguments are wavelength range of final spectrum
isp.resample_spectrum(spectrum_id, 424.76, 763.21)

isp.sendoff(spectrum_id)




