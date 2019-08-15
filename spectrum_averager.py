import scipy.io as sio
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

# Loads matlab cell arrays for x and y axes of spectrum
mat_contents = sio.loadmat('fulldata.mat')
fullwave = mat_contents['fullwave']
fullint = mat_contents['fullint']

# ----------USER INPUT SECTION----------
starname = 'test'
# Must be manually changed to match dimensions of matlab arrays
cols = 169317
rows = 15
# Display spectrum after calculation?
show_spectrum = 1
# --------------------------------------

# Creates blank lists to be overwritten by matlab data
x = np.zeros(cols)
y = np.zeros(cols)
err = np.zeros(cols)

# Iterate across cell arrays, write x values to list
for j in range(0, cols):
    # Divide by 10 to convert A to nm
    x[j] = (fullwave[0, j]/10)
    tempsum = 0
    diffsum = 0
    # Write averaged y values to list
    for i in range(0, rows):
        tempsum = tempsum + fullint[i, j]
    y[j] = tempsum/rows
    # Write standard error to list
    for k in range(0, rows):
        diffsum = diffsum + (fullint[k, j]-y[j])**2
    sigma = sqrt(diffsum/(rows-1))
    err[j] = sigma/sqrt(rows)

# Create transposed array of data for writing to text file
data = np.array([x, y, err])
data = data.T

# Write tab-spaced text file with headers
with open(starname+'_datafile.txt', 'w+') as datafile_id:
    datafile_id.write("waveobs\tflux\terr\n")
    np.savetxt(datafile_id, data, fmt=['%s', '%s', '%s'], delimiter='\t')
    
if show_spectrum == 1:
    plt.plot(x, y)
    plt.title(starname)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux')
    plt.show()
