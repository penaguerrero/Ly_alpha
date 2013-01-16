import os
import glob
import numpy
import sys
import table
import re
import string
from scipy.optimize import leastsq
from spectrum import spectrum
from spectrum import pyplot
from configobj import ConfigObj

def selection(A, cgs, lower, upper):
    # This functions selects according to a wavelength range
    xlist = []
    ylist = []
    for i in range(0, len(A)):
        if A[i] >= lower and A[i] <= upper:
            xlist.append(A[i])
            ylist.append(cgs[i])     
    A_selected = numpy.array(xlist)
    cgs_selected = numpy.array(ylist)
    return(A_selected, cgs_selected)

def get_continuum(arr):
    # create the full baseline spectrum...
    # arr = wavelength and flux array
    # baselinepars are the continuum parameters for a polynomial of 2nd degree:
    # baseline[0]*x^2 + baseline[1]*x + baseline[2]
    # and in case of a power law:
    # baselinepars[0] = scaling constant
    # baseline[1] = exponent
    # baseline[2] = re-scaling factor for the array
    n = len(arr[0])
    lnx = numpy.log(arr[0])
    lny = numpy.log(arr[1])
    multip = [lnx * lny]
    # initial guess for baseline or continuum parameters
    baselinepars = [1, 1, 1]
    baselinepars[1] = (n*sum(multip) - sum(lnx)*sum(lny)) / (n*sum(lnx)**2 - (sum(lnx))**2)
    b = baselinepars[1] 
    a = (sum(lny) - b * sum(lnx)) / n
    baselinepars[0] = numpy.exp(a)
#    if powerlaw:
#        arr[0].powerlaw = True
    cont_y = (0.95*baselinepars[0]*(arr[0]/baselinepars[2])**(-baselinepars[1])).squeeze()
    continuum = numpy.array([arr[0], cont_y])
    return(continuum)
#    else:
        #arr[0].powerlaw = False
#    return numpy.poly1d(baselinepars)(arr[0])

def area_under_curve(x1, x2): # Using Riemann Sum
    if x1 > x2:
        print('The calculated area will be negative.')
    # Compute step of integration
    delta_x = ((x2-x1)/1000)
    j = abs((x2-x1)/delta_x)
    i = int(j)
    print('Number of smaller areas to be summed =', i)
    # initialize
    n = 0       # number of steps
    tot_A = 0.0 # total area
    x = x1      # starting point
    # Begin numerical integration
    while n < i:
        delta_A = x**2 * delta_x
        tot_A = tot_A + delta_A
        n = n + 1
    return tot_A
        
# Path to files
Bstarpath = '../Tlusty/solar/BGuvspec_v10/'
Bstarpath2 = '../Tlusty/solar/BGuvspec_v2/'

# Detailed spectrum files
#specs = glob.glob(Bstarpath+'*.uv.7') # microturbulent velocity of 10
#specs = glob.glob(Bstarpath2+'*.uv.7')
specs = ['BL30000g425v2.uv.7', 'BL30000g425v2.uv.7']

# Detailed continuum files
#conts = glob.glob(Bstarpath+'*.uv.17') # microturbulent velocity of 10
#conts = glob.glob(Bstarpath2+'*.uv.17')
conts = ['BL30000g425v2.uv.17', 'BL30000g425v2.uv.17']

path_results = '../results'

# MEASUREMENTS
# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701

# The selection range
lowav = 1100.0
upwav = 1350.0

'''
### CREATING TEXT FILES OF SELECTED RANGES
total_files = len(specs)
print("%d files" % (total_files))
for i in range(0, total_files):
    #bstar_result_path = os.path.join(path_results, 'Norm_Bstarv2_' + repr(i) + '.txt')
    bstar_result_path = os.path.join(path_results, 'Norm_BL30000g425v2.txt')
    A, cgs = numpy.loadtxt(specs[i], dtype=float, unpack=True)
    A_cont, cgs_cont = numpy.loadtxt(conts[i], dtype=float, unpack=True)
    print("%d of %d: [ %s, %s ] ... " % (i+1,total_files, os.path.basename(specs[i]), os.path.basename(conts[i])))
    # Selecting the range for Ly_alpha
    wav, flx = selection(A, cgs, lowav, upwav)       
    wav_cont, flx_cont = selection(A_cont, cgs_cont, lowav, upwav)
    # Rebinning
    spec_arr = numpy.array([wav, flx])
    cont_arr = numpy.array([wav_cont, flx_cont])    
    rebin_spec, rebin_cont, cont_factor, spec_factor = spectrum.get_factors_and_rebin(spec_arr, cont_arr)
    # Dividing spectra over continuum
    norm_flx = rebin_spec[1] / rebin_cont[1]
    # Writting to text files
    txtout = open(bstar_result_path, 'w+')
    for j in range(0, len(norm_flx)):
        txtout.write((rebin_spec[0][j].__repr__() + " " + norm_flx[j].__repr__()) + "\n")
    txtout.close()
    print("\n")
    
'''
### Calling the selected and converted files
norm_stars = glob.glob(path_results+'/Norm_*.txt')
for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    wav, flx = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    spec_data = numpy.array([wav, flx])
    eqw = area_under_curve(1200, 1232)
    print('Ly-alpha EQW =', eqw)
    '''
    ## Getting the continuum to measure EQW
    conti = get_continuum(spec_data)
    ## THIS IS ONLY FOR THE NORMALIZED SPECTRA, WHERE THE THEORETICAL CONTINUUM APPLIES
    continuum_norm = spectrum.theo_cont(wav, scale_factor=1.0)
    # Ploting just to check
    # the line to mark Ly_alpha
    lyalpha_arr_wav = []
    for j in wav:
        lyalpha_arr_wav.append(Lyalpha)            
    # Plot limits
    low = 1160
    up = 1260
    # Figure
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.suptitle(norm_stars[i])
    pyplot.plot(wav, flx, 'b', conti[0], conti[1], 'magenta')
    pyplot.plot(lyalpha_arr_wav, flx, 'r--')
    pyplot.show()
    '''


        