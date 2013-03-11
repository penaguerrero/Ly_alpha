import os
import glob
import numpy
import string
import table
from spectrum import spectrum
from spectrum import pyplot

'''
This program produces text files of a size that splot can work with.

### I decided to do separately turbulent velocities v10 and v2 because they have 97 and 163 files, respectively, so if an error occurs the whole
program stops AND because originally, I made them separately and I want to be able to compare with previous results. Names at first were just

'''

# directory of normalized names
normBstars = []
normBstars.append('Theoretical Tlusty B star')
  
all_eqws = []
all_eqws.append('Fixed 10A eqw')
VG_eqw = []
VG_eqw.append('Valls-Gabaud (1993) eqw')
expected_eqwstoVG = []
expected_eqwstoVG.append('stable approx Valls-Gabaud93 eqw')
half_eqwstimes2 = []
half_eqwstimes2.append('half eqw*2 at 1250')
half_eqwstimes2_avoiding_N = []
half_eqwstimes2_avoiding_N.append('half eqw*2 at 1237')

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
norm_stars = glob.glob(path_results+'Norm_Bstar*')
#norm_stars = path_results+'Norm_Bstar_136_v2_27000g450.txt'

for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    spec_data = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    continuum_norm = spectrum.theo_cont(spec_data[0], scale_factor=1.0)
    eqw = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width=10.0)
    print('*** Ly-alpha EQW = %f' % eqw)
    all_eqws.append(eqw)
    normBstars.append(os.path.basename(norm_stars[i]))
    
    base_name = os.path.basename(norm_stars[i])
    base_name2 = string.split(base_name, sep='_')
    Teff, base_name3 = string.split(base_name2[4], sep='g')
    base_name4 = base_name3.replace(".txt", "")
    logg = base_name4.replace("g", "")
    print(Teff, logg)
    VG_eqw_object, expected_EQW, expected_lolim, expected_uplim = spectrum.EQW_initial_guessVG(spec_data, continuum_norm, Lyalpha, Teff, guessed_EQW=1)
    print ('B star theoretical EQW = %f' % expected_EQW)
    VG_eqw.append(VG_eqw_object)
    expected_eqwstoVG.append(expected_EQW)

    half_eqw_times2_obj = spectrum.half_EQW_times2(spec_data, continuum_norm, Lyalpha, wave_limit=1250, right_side=True)
    print('half_eqw_times2 = %f' % half_eqw_times2_obj)
    half_eqwstimes2.append(half_eqw_times2_obj)
    
    half_eqw_times2_avoiding_N_obj = spectrum.half_EQW_times2(spec_data, continuum_norm, Lyalpha, wave_limit=1237, right_side=True)
    print('half_eqw_times2_avoiding_N_obj = %f' % half_eqw_times2_avoiding_N_obj)
    half_eqwstimes2_avoiding_N.append(half_eqw_times2_avoiding_N_obj)
    print ' '
    
'''
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


# Printing file of normalized B star and corresponding eqw:
eqwf = open(path_results+'eqwsBstars_v10.txt', 'w+')
print >> eqwf, 'B star', 'Ly-alpha EQW'
print >> eqwf, '*** SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
    print >> eqwf,  normBstars[i], all_eqws[i]
eqwf.close()
 

# Printing a file with the various eqws
print ('normBstars, all_eqws, half_eqwstimes2, half_eqwstimes2_avoiding_N, VG_eqw, expected_eqwstoVG',
       numpy.shape(normBstars), numpy.shape(all_eqws), numpy.shape(half_eqwstimes2), numpy.shape(half_eqwstimes2_avoiding_N),
       numpy.shape(VG_eqw), numpy.shape(expected_eqwstoVG))

comp_eqws = open(path_results+'eqws_COMPARISON.txt', 'w+')
table_comp = []
print >> comp_eqws, 'SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
#    print >> comp_eqws, normBstars[i], all_eqws[i], half_eqwstimes2[i], half_eqwstimes2_avoiding_N[i], VG_eqw[i], expected_eqwstoVG[i]
    temp = [normBstars[i], repr(all_eqws[i]), repr(half_eqwstimes2[i]), repr(half_eqwstimes2_avoiding_N[i]), 
            repr(VG_eqw[i]), repr(expected_eqwstoVG[i])] 
    table_comp.append(temp)
table.pprint_table(comp_eqws, table_comp)
comp_eqws.close()

'''


        