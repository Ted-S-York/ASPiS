# ASPiS
Automated Spectrum Processing with iSpec

ABOUT

ASPiS is a set of python scripts for automation of the analysis of stellar spectra using iSpec (https://www.blancocuaresma.com/s/iSpec). It is designed to be used with spectra processed via the MEGARA reduction pipeline.

REQUIREMENTS

ASPiS requires python 2.7 with the same dependencies as iSpec (https://www.blancocuaresma.com/s/iSpec/manual/installation/linux/anacondapython). ASPiS was developed on Ubuntu and should be compatible with Debian-based operating systems, compatibility with other OS may be possible but is not guaranteed.

INSTALLATION AND USAGE

ASPiS can be installed by simply downloading the included files and placing them in a directory of choice, this does NOT need to be the same as the iSpec install directory. Once they are in place, the user must change the ispec_dir variable at the top of ispec_functions.py to match the iSpec directory path.

SPECTRUM_AVERAGER.PY:

This script takes the full reduced spectral data from MEGARA and produces a single median spectrum that can be parsed by iSpec. This requires a MEGARA output file that must be titled 'fulldata.mat' and should contain the wavelength (fullwave), flux (fullint) and Julian date (jd) information for all observations. A copy of fulldata.mat must be placed in the ASPiS directory for spectrum_averager.py to function. 

If desired, the data can be 'binned' according to pulsation phase for variable stars. To enable this, change the variable use_spectrum_binning to True and then input the duration of variability desired in days for the period variable and define the desired number of bins via the no_bins variable. It is recommended to use at least 4 bins, though higher numbers will increase computation time.

SPECTRUM_PREP.PY:

This script applies various corrections to the spectrum datafile produced by spectrum_averager.py. It will also produce a .fits file that will be used for analysis, and a info .txt file that results will be appended to. A .txt file of the correct format outputted by spectrum_averager must be present in the ASPiS directory. 

The degree of telluric line removal is specified by the variable telluric_degree, by default this is set to 0 for complete removal but can be any value up to 100, which represents no removal. The initial and final resolutions for resolution degradation are defined by the variables start_resolution and final_resolution respectively. The variables med_wave and max_wave determine the median and maximum wave ranges for continuum normalisation and lambda_min and lambda_max determine the wavelength range for the resampled spectrum in nm.

SPECTRUM_ANALYSER.PY:

This script takes the processed .fit file and analyses it to determine atmospheric parameters and stellar abundances. This requires an appropriate .fit file and accompanying info.txt file to be located in the ASPiS directory. 

Lists of detected linemasks can be saved during analysis if desired, otherwise the linemask files can be deleted after use by setting the purge_linemasks variable to True. The purge_log variable will enable or disable the deletion of the log file produced by the logging library at the end of each run. A minimum threshold value of detected linemasks can be set for abundance analysis using the line_threshold variable, any species that fails to meet this requirement will not be analysed to decrease computation time. The list of species to be analysed is defined in species_list; target species must be included here as comma-separated strings of the form "Fe 2", "Ca 1" etc.

Atmospheric parameters, chemical abundances and their errors will be appended to the info.txt file. The computation of abundance error is multiprocessed and by default will make use of as many cores as are available on the machine. If it is desired to restrict this for any reason, add the maximum number of cores as an argument to line 567 of ispec_functions.py. Linelists, model atmospheres etc. may also be selected by modifying ispec_functions.py.
