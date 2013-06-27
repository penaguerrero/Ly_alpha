import os
import glob
import numpy
import string
import table
from spectrum import spectrum
from matplotlib import pyplot

'''
The focus of this code is for B stars with different metallicities to see hoe the EQW changes.

This program produces text files of a size that splot can work with, and that is the right resolution: delta_lambda=0.465

It also measures the EQW at different fixed widths and prints a nice table.
'''

# Path to files
BstarsZs = glob.glob('../Tlusty/DiffZv2/*')
path_results = '../results/Bstars30KdiffZ/rebinned_specs/'
path_plots = '../results/plots/'

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701
# width of Lyman-alpha
width = 10.0

# The selection range
lowav = 1100.0
upwav = 1301.0

'''
for Bstarpath in BstarsZs:
    print('In directory: %s' % Bstarpath)
    # Detailed spectrum files
    specs = glob.glob(Bstarpath+'/*.uv.7')
    #specs = glob.glob(Bstarpath2+'BG28000g275v2.uv.7')
    
    # Detailed continuum files
    conts = glob.glob(Bstarpath+'/*.uv.17') # microturbulent velocity of 2 km/s
    #conts = glob.glob(Bstarpath2+'BG28000g275v2.uv.17')        
    
    ### CREATING TEXT FILES OF SELECTED RANGES
    total_files = len(specs)
    print("%d files" % (total_files))
    
    #### 
    #### This part of the code creates the normalized text files with the desired resolution.
    ####
    
    ### CREATING TEXT FILES OF SELECTED RANGES
    total_files = len(specs)
    print("%d files" % (total_files))
    
    for i in range(0, total_files):
        base_name = os.path.basename(specs[i])
        print(base_name)
        base_name2 = string.split(base_name, sep='.')
        base_name3 = string.split(base_name2[0], sep='3')
        if base_name3[0] == 'BC':
            Zkey = 'Z2'
        elif base_name3[0] == 'BG':
            Zkey = 'Z1'
        elif base_name3[0] == 'BL':
            Zkey = 'Z05'
        elif base_name3[0] == 'BS':
            Zkey = 'Z5th'
        elif base_name3[0] == 'BT':
            Zkey = 'Z10th'
        Bstar = base_name2[0]+Zkey
        templogg = '30000g425'
        
        bstar_result_path = os.path.join(path_results, 'Norm_Bstar_' + Zkey + '_v2_' + templogg +'.txt')
        A, cgs = numpy.loadtxt(specs[i], dtype=float, unpack=True)
        A_cont, cgs_cont = numpy.loadtxt(conts[i], dtype=float, unpack=True)
        print("%d of %d: [ %s, %s ] ... " % (i+1, total_files, os.path.basename(specs[i]), os.path.basename(conts[i])))
            
        # Selecting the range for Ly_alpha
        wav, flx = spectrum.selection(A, cgs, lowav, upwav)       
        wav_cont, flx_cont = spectrum.selection(A_cont, cgs_cont, lowav, upwav)
    
        # Rebinning
        spec_arr = numpy.array([wav, flx])
        cont_arr = numpy.array([wav_cont, flx_cont])    
        print('shape of spec_arr: %s    ----   shape of cont_arr: %s' % (repr(spec_arr.shape), repr(cont_arr.shape)))
        # Rebinning data to delta_lambda=0.465, which is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)
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
     
print('Done!')
'''


'''
### This part of the code is calling the text files already created (the REBINNED files) in order to do the measurements: EQWs 
'''
    
# directory of normalized names
normBstars = []
normBstars.append('Theoretical Tlusty B star')
  
eqws_10 = []
eqws_10.append('Ly-alpha +- 10A eqw')
eqws_20 = []
eqws_20.append('Ly-alpha +- 20A eqw')
eqws_30 = []
eqws_30.append('Ly-alpha +- 30A eqw')
half_eqwstimes2_avoiding_N = []
half_eqwstimes2_avoiding_N.append('half eqw*2 at 1235')
half_eqwstimes2_avoiding_Si = []
half_eqwstimes2_avoiding_Si.append('(LyA+20)*2')

### Calling the selected and converted files
norm_stars = glob.glob(path_results+'Norm_Bstar*30000g425.txt')
#norm_stars = path_results+'Norm_Bstar_Z2_v2_30000g425.txt'

for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    spec_data = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    continuum_norm = spectrum.theo_cont(spec_data[0], scale_factor=1.0)
    
    eqw1 = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width=10.0*2)
    print('*** Ly-alpha EQW fixed at 10 A = %f' % eqw1)
    eqws_10.append(eqw1)
    normBstars.append(os.path.basename(norm_stars[i]))
    
    eqw2 = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width=20.0*2)
    print('*** Ly-alpha EQW fixed at 20 A = %f' % eqw2)
    eqws_20.append(eqw2)
    
    eqw3 = spectrum.EQW_line_fixed(spec_data, continuum_norm, Lyalpha, width=30.0*2)
    print('*** Ly-alpha EQW fixed at 30 A = %f' % eqw3)
    eqws_30.append(eqw3)    
    
    base_name = os.path.basename(norm_stars[i])
    base_name2 = string.split(base_name, sep='_')
    Teff, base_name3 = string.split(base_name2[4], sep='g')
    base_name4 = base_name3.replace(".txt", "")
    logg = base_name4.replace("g", "")
    print(Teff, logg)
       
    half_eqw_times2_avoiding_N_obj = spectrum.half_EQW_times2(spec_data, continuum_norm, Lyalpha, wave_limit=1235, right_side=True)
    print('half_eqw_times2_avoiding_N_obj = %f' % half_eqw_times2_avoiding_N_obj)
    half_eqwstimes2_avoiding_N.append(half_eqw_times2_avoiding_N_obj)

    half_eqw_times2_avoiding_Si_obj = spectrum.half_EQW_times2(spec_data, continuum_norm, Lyalpha, wave_limit=(Lyalpha+20.0), right_side=True)
    print('half_eqwstimes2_avoiding_Si_obj = %f' % half_eqw_times2_avoiding_Si_obj)
    half_eqwstimes2_avoiding_Si.append(half_eqw_times2_avoiding_Si_obj)
    print ' '
    
    
    # Ploting just to check
    # the line to mark Ly_alpha
    lya, _ = spectrum.find_nearest(spec_data[0], Lyalpha)
    flx_lya = spectrum.findXinY(spec_data[1], spec_data[0], lya)
    if flx_lya > 1.0:
        y_Lyalpha = [flx_lya+0.1, flx_lya+0.2]
    else:
        y_Lyalpha = [0.8, 1.05]
    lyalpha_arr_norm = []
    for i in range(0, len(y_Lyalpha)):
        lyalpha_arr_norm.append(Lyalpha)
    lyalpha_arr = numpy.array([lyalpha_arr_norm, y_Lyalpha])
    # Plot limits
    low = 1160
    up = 1260
    # Figure
    pyplot.figure(1, figsize=(10, 10))
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.ylim(-0.1, 1.1)
    #pyplot.suptitle(norm_stars[i])
    pyplot.plot(spec_data[0], spec_data[1], 'b', lyalpha_arr[0], lyalpha_arr[1], 'r--')#, continuum_norm[0], continuum_norm[1], 'magenta')
    epsfile = os.path.join(path_plots, "Bstar_"+base_name+".eps")
    pyplot.savefig(epsfile)    
    pyplot.show()

    print("Do you want to save this plot?  [y/N , meaning NO is default... :P ]")
    save_plt = raw_input()
    pyplot.ioff()
    if save_plt == 'n' or save_plt == '':
        print('Plot not saved...   :( ')
        os.remove(epsfile)
    else:
        print('Star %s was saved.' % base_name)



'''
# Printing a file with the various eqws
comp_eqws = open(path_results+'kk.txt', 'w+')
table_comp = []
print >> comp_eqws, 'For an early B star (all taken from Tlusty website) with:'
print >> comp_eqws, 'SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
#    print >> comp_eqws, normBstars[i], eqws_10[i], eqws_20[i], eqws_30[i], half_eqwstimes2_avoiding_N[i]
    temp = [normBstars[i], repr(eqws_10[i]), repr(eqws_20[i]), repr(eqws_30[i]), 
            repr(half_eqwstimes2_avoiding_N[i]), repr(half_eqwstimes2_avoiding_Si[i])] 
    table_comp.append(temp)
table.pprint_table(comp_eqws, table_comp)
comp_eqws.close()

print('Done!')
'''

