import scipy.io as sio
import numpy as np
from math import sqrt

#Loads matlab cell arrays for x and y axes of spectrum
print('Loading fulldata.mat...')
mat_contents = sio.loadmat('fulldata.mat')
fullwave = mat_contents['fullwave']
fullint = mat_contents['fullint']
jul_date = mat_contents['jd']

#-----------------------------USER INPUT SECTION------------------------------#
starname = raw_input('Please define the title to be used for this spectrum: ')

use_spectrum_binning = False

if use_spectrum_binning:
    #Input target pulsation period in days
    period = 1
    #Define desired number of spectral bins, 4 is a recommended minimum.
    no_bins = 4
    zero_point = jul_date[0, 0]

#-----------------------------------------------------------------------------#

#Get size of matlab arrays
print('Getting size of matlab arrays...')
cols = fullwave.shape[1]
rows = len(fullint)

#Creates blank lists to be overwritten by matlab data
x = []
if use_spectrum_binning:
    print('Creating bins...')
    y = [[] for n in range(no_bins)]
#    err = [[] for n in range(no_bins)]
    err = [np.zeros(cols) for n in range(no_bins)]
else:
    y = []
    err = []
    
#Iterate across cell arrays, write x values to list
print('Writing wavelength values...')
for i in range(cols):
    x.append(fullwave[0, i]/10)
    
#Calculate and write median y values to list WITHOUT binning
print('Writing flux values and errors...')
if not use_spectrum_binning:
    for j in range(cols):
        tempsum = []
        diffsum = 0
        for i in range(rows):
            tempsum.append(fullint[i, j])
        y.append(np.median(tempsum))
        for k in range(rows):
            diffsum += (fullint[k, j]-y[j])**2
        sigma = sqrt(diffsum/(rows-1))
        err.append(sigma/sqrt(rows))
    
    #Write unbinned textfile
    print('Saving datafile...')
    data = np.array([x, y, err])
    data = data.T
    
    with open(starname+'_datafile.txt', 'w+') as datafile_id:
        datafile_id.write("waveobs\tflux\terr\n")
        np.savetxt(datafile_id, data, fmt=['%s', '%s', '%s'], delimiter='\t')
    
#Calculate and write median y values to list WITH binning
else:
    #Analyse jd to sort observations into bins
    for j in range(cols):
        tempsum = [[] for n in range(no_bins)]
        for i in range(rows):
            #Determine phase of observation
            phase = jul_date[0, i] - zero_point
            #Find corresponding bin
            for n in range(no_bins):
                #Case for zero point observation
                if phase == 0:
                    active_bin = 0
                #Case for phases that are multiples of the period
                elif phase%period == 0 and n == no_bins - 1:
                    active_bin = n
                #Case for all other phases
                elif phase%period > n*period/(no_bins) and phase%period <= (n+1)*period/(no_bins):
                    active_bin = n
            tempsum[active_bin].append(fullint[i, j])
        for m in range(no_bins):
            #Write median int value to bin
            y[m].append(np.median(tempsum[m]))
            #Calculate standard error
#            diffsum = 0
#            for k in range(len(tempsum[m])):
#                diffsum += (tempsum[m][k]-y[m])**2
            #Write standard error to bin
#            sigma = np.sqrt(diffsum/(len(tempsum[m])-1))
#            err[m].append(sigma/np.sqrt(len(tempsum[m])))
                    
    #Write binned datafiles
    print('Saving datafiles...')
    for n in range(no_bins):
        data = np.array([x, y[n], err[n]])
        data = data.T
        title = starname + str(n+1) + '_datafile.txt'
        with open(title, 'w+') as datafile_id:
            datafile_id.write("waveobs\tflux\terr\n")
            np.savetxt(datafile_id, data, fmt=['%s', '%s', '%s'], delimiter='\t')

print('Complete')

    


