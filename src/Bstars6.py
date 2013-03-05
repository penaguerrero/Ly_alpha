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

#for Bstarpath in Bstarpath_v10andv2:
#print(Bstarpath)

# Detailed spectrum files
#specs = glob.glob(Bstarpath+'*.uv.7') # microturbulent velocity of 10
specs = glob.glob(Bstarpath2+'*.uv.7')
#specs = glob.glob(Bstarpath2+'BG28000g275v2.uv.7')

# Detailed continuum files
#conts = glob.glob(Bstarpath+'*.uv.17') # microturbulent velocity of 10
conts = glob.glob(Bstarpath2+'*.uv.17')
#conts = glob.glob(Bstarpath2+'BG28000g275v2.uv.17')

path_results = '../results/BTlustyStars/'

# MEASUREMENTS
# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701
# width of Lyman-alpha
width = 10.0

# The selection range
lowav = 1100.0
upwav = 1301.0


### CREATING TEXT FILES OF SELECTED RANGES
total_files = len(specs)
print("%d files" % (total_files))
v2 = 'v2'
v10 = 'v10'
    
for i in range(0, total_files):
    base_name = os.path.basename(specs[i])
    base_name2 = string.split(base_name, sep='.')
    base_name3 = base_name2[0].replace("BG", "")
    templogg = base_name3.replace("v2", "")
    #print(templogg)
    if v2 in specs[i]:
        bstar_result_path = os.path.join(path_results, 'Norm_Bstar_' + repr(i) + '_v2_' + templogg +'.txt')
        #bstar_result_path = os.path.join(path_results, 'Norm_BL30000g425v2.txt')
    elif v10 in specs[i]:
        bstar_result_path = os.path.join(path_results, 'Norm_Bstar_' + repr(i) + '_v10_' + templogg + '.txt')
        #bstar_result_path = os.path.join(path_results, 'Norm_BL30000g425v10.txt')
    #A, cgs = numpy.loadtxt(specs[i], dtype=float, unpack=True)
    #A_cont, cgs_cont = numpy.loadtxt(conts[i], dtype=float, unpack=True)
    #print("%d of %d: [ %s, %s ] ... " % (i+1,total_files, os.path.basename(specs[i]), os.path.basename(conts[i])))
    
    aka = 'Bstar_'+repr(i)+'.v2.'+templogg
    normBstars.append(aka)
    '''
    # Selecting the range for Ly_alpha
    wav, flx = spectrum.selection(A, cgs, lowav, upwav)       
    wav_cont, flx_cont = spectrum.selection(A_cont, cgs_cont, lowav, upwav)

    # Rebinning
    spec_arr = numpy.array([wav, flx])
    cont_arr = numpy.array([wav_cont, flx_cont])    
    print('shape of spec_arr: %s    ----   shape of cont_arr: %s' % (repr(spec_arr.shape), repr(cont_arr.shape)))
    # Rebinning data to delta_lambda=0.456, which is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)
    desired_delta_lambda = 0.465
    _, initial_rows = cont_arr.shape
    rebin_spec, rebin_cont = spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, spec_arr, cont_arr, initial_rows-1)
    print('shape of rebin_spec: %s    ----   shape of rebin_cont: %s' % (repr(rebin_spec.shape), repr(rebin_cont.shape)))

    # Dividing spectra over continuum
    norm_flx = rebin_spec[1] / rebin_cont[1]
    
    # Writting to text files
    txtout = open(bstar_result_path, 'w+')
    for j in range(0, len(norm_flx)):
        if rebin_spec[0][j] != 0.0:
            txtout.write((rebin_spec[0][j].__repr__() + " " + norm_flx[j].__repr__()) + "\n")
    txtout.close()
    print("\n")
    '''    

### Calling the selected and converted files
norm_stars = glob.glob(path_results+'Norm_Bstar*v2*')
for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    spec_data = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    ## Getting the continuum to measure EQW
    #conti = get_continuum(spec_data)
    ## THIS IS ONLY FOR THE NORMALIZED SPECTRA, WHERE THE THEORETICAL CONTINUUM APPLIES
    continuum_norm = spectrum.theo_cont(spec_data[0], scale_factor=1.0)
    #print('data', spec_data)
    ### Getting EQW
    eqw = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width)
    print('*** Ly-alpha EQW = %f' % eqw)
    all_eqws.append(eqw)
    
    print(norm_stars[i])
    #teqw = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, variable_width)
    #print('TEST Ly-alpha EQW = %f' % teqw)

    eqw_iter, lo, up = spectrum.EQW_iter(spec_data, continuum_norm, Lyalpha, guessed_width=2.0)
    print('eqw_iter = %f, lo = %f, up = %f' % (eqw_iter, lo, up))
    
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


        