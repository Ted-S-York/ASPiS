import ispec_functions as isp
import os
import re

#---------------------------USER INPUT & OPTIONS------------------------------#

target_id = raw_input('Please define the target to analyse: ')
spectrum_id = target_id + '_processed.fits'
info_id = target_id + '_info.txt'

#if True, logging file will be deleted after run.
purge_log = True

#if True, linemask files will be deleted after each run.
purge_linemasks = True

#Species with fewer lines than this will not be analysed.
line_threshold = 1

#Species to be analysed should be listed here in the format "Fe 1" etc.
species_list = ["C 1","N 1","O 1","Na 1","Na 2","Mg 1","Mg 2","Si 1","Si 2",
                "S 1","S 2","Ca 1","Ca 2","Sc 1","Sc 2","Ti 1","Ti 2","V 1",
                "V 2","Cr 1","Cr 2","Mn 1","Mn 2","Fe 1","Fe 2","Co 1","Co 2",
                "Ni 1","Ni 2","Cu 1","Cu 2","Zn 1","Zn 2","Sr 1","Sr 2","Y 1",
                "Zr 1","Zr 2","Ba 1","Ba 2","La 1","La 2","Ce 1","Ce 2","Pr 1",
                "Pr 2","Nd 1","Nd 2","Sm 1","Sm 2","Eu 1","Eu 2","Gd 1", "Gd 2"]

#-------------------------ATMOSPHERIC PARAMETERS------------------------------#

#Find iron lines and determine atmospheric parameters
isp.find_fe_lines(spectrum_id)
isp.determine_parameters(spectrum_id, 'fe_linemasks.txt')

if purge_linemasks:
    os.remove("fe_linemasks.txt")

#----------------------WRITING ATM. PARAMS TO INFO FILE----------------------#

temp_id = spectrum_id + "_info.txt"

with open(temp_id, "r") as file:
    filedata = file.read()

filedata = filedata.replace(',', '\n')

with open(info_id, "a") as file:
    file.write("\n#-----ATMOSPHERIC PARAMETERS AND ERRORS-----#")
    file.write(filedata)
    file.write("\n#------CHEMICAL ABUNDANCES AND ERRORS-------#")
   
os.remove(temp_id)

#-----------------------READ PARAMETERS FROM INFO FILE-----------------------#

with open(info_id, "r") as file :
    params_data = file.read()

params = re.findall(r"[-+]?\d*\.\d+|\d+", params_data)
del params[0:2]
del params[8:]

for i in range(0, len(params)):
    params[i] = float(params[i])
       
#-------------------------------ANALYSE SPECIES------------------------------#

for x in range(0, len(species_list)):
    print("Beginning analysis for %s" % species_list[x])
    #Find linemasks of target species.
    isp.find_linemasks(spectrum_id, species_list[x])
    #Create linemasks file.
    file_name = species_list[x].replace(" ", "_") + "_linemasks.txt"
    #Species with fewer linemasks than the threshold will not be analysed.
    lines_detected = (sum(1 for line in open(file_name))-1)
    if lines_detected < line_threshold:
        print("Species does not have sufficient spectral lines, discarding...")
        os.remove(file_name)
    else:
        print("Detected %d spectral lines" % lines_detected)
        #Determine abundance.
        initial_abnd = isp.determine_abundances(spectrum_id, species_list[x], params)
        abundance = initial_abnd[0]
        scatt_err = initial_abnd[1]
        #DETERMINING ERRORS
        error = isp.determine_error(spectrum_id, species_list[x], params, abundance, scatt_err)
        #Save result + error to file.
        with open(info_id, "a") as file:
            line1 = "\n"+species_list[x]+": "+str(abundance)+" +/- "+str(error)+" N = "+str(lines_detected)
            file.write(line1)
        #Clean up linemasks file if desired.
        if purge_linemasks:
            os.remove(file_name)

#-----------------------------------CLEANUP-----------------------------------#

if purge_log:
    print("Purging log file...")
    os.remove("ispec.log")
    print("Log purged")
    