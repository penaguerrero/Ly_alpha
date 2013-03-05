import os
import glob
import numpy
import string
from spectrum import spectrum
from spectrum import pyplot

'''
This program produces text files of a size that splot can work with.

### I decided to do separately turbulent velocities v10 and v2 because they have 97 and 163 files, respectively, so if an error occurs the whole
program stops AND because originally, I made them separately and I want to be able to compare with previous results. Names at first were just

'''

# directory of normalized names
normBstars = []
  
all_eqws = []

# Path to files
#Bstarpath_v10andv2 = glob.glob('../Tlusty/solar/B*/')    # this is to get all the files with different turbulent velocities
Bstarpath = '../Tlusty/solar/BGuvspec_v10/'
Bstarpath2 = '../Tlusty/solar/BGuvspec_v2/'

# Detailed spectrum files
#specs = glob.glob(Bstarpath+'*.uv.7') # microturbulent velocity of 10
specs = glob.glob(Bstarpath2+'*.uv.7')
#specs = glob.glob(Bstarpath2+'BG28000g275v2.uv.7')


path_results = '../results/BTlustyStars/'

# MEASUREMENTS
# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701


### CREATING TEXT FILES OF SELECTED RANGES
total_files = len(specs)
print("%d files" % (total_files))

### Calling the selected and converted files
norm_stars = glob.glob(path_results+'Norm_Bstar*v2*')
for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    spec_data = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    continuum_norm = spectrum.theo_cont(spec_data[0], scale_factor=1.0)
    eqw = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width=10.0)
    print('*** Ly-alpha EQW = %f' % eqw)
    all_eqws.append(eqw)
    
    base_name = os.path.basename(norm_stars[i])
    base_name2 = string.split(base_name, sep='_')
    Teff, base_name3 = string.split(base_name2[4], sep='g')
    base_name4 = base_name3.replace(".txt", "")
    logg = base_name4.replace("g", "")
    print(Teff, logg)
    expected_EQW, expected_lolim, expected_uplim = spectrum.EQW_initial_guess(spec_data, continuum_norm, Lyalpha, Teff, guessed_EQW=1)
    print ('expected_EQW', expected_EQW)

    #eqw_iter, lo, up = spectrum.EQW_iter(spec_data, continuum_norm, Lyalpha, guessed_width=2.0)
    #print('eqw_iter = %f, lo = %f, up = %f' % (eqw_iter, lo, up))
    
    print ' '

    # Ploting just to check
    # the line to mark Ly_alpha
    lyalpha_arr_w = []
    for j in spec_data[0]:
        lyalpha_arr_w.append(Lyalpha)            
    # Plot limits
    low = 1160
    up = 1260
    # Figure
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.ylim(-0.1, 1.1)
    pyplot.suptitle(norm_stars[i])
    pyplot.plot(spec_data[0], spec_data[1], 'b', continuum_norm[0], continuum_norm[1], 'magenta', lyalpha_arr_w, spec_data[1], 'r--')
    pyplot.show()

    #exit()

'''
# Printing file of normalized B star and corresponding eqw:
eqwf = open(path_results+'eqwsBstars_v10.txt', 'w+')
print >> eqwf, 'B star', 'Ly-alpha EQW'
print >> eqwf, '*** SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
    print >> eqwf,  normBstars[i], all_eqws[i]
eqwf.close()
 
'''


        