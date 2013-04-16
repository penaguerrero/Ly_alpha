import os
import glob
import numpy
import string
import table
from spectrum import spectrum
from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator
from pylab import *


'''
This program produces text files of a size that splot can work with, and that is the right resolution: delta_lambda=0.465

### I decided to do separately turbulent velocities v10 and v2 because they have 97 and 163 files, respectively, so if an error occurs the whole
program stops AND because originally, I made them separately and I want to be able to compare with previous results. Names at first were just

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



path_results = '../results/BTlustyStars/'
path_plots = '../results/plots/'

# MEASUREMENTS
# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.3 #1215.6701
# rest wavelengths from Leitherer et al(2011), ApJ, 141, 37
# the numbers are different because they look prettier in the plot
cIII_mie = 1172.0 #1175.53
siIII = 1202.3 #1206.50
nV = 1242.80
cIII = 1247.38
sII = 1250.58


### CREATING TEXT FILES OF SELECTED RANGES
#total_files = len(specs)
#print("%d files" % (total_files))

### Calling the selected and converted files
norm_stars = glob.glob(path_results+'Norm_Bstar*')
#norm_stars = path_results+'Norm_Bstar_136_v2_27000g450.txt'

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
    '''
    lyalpha_arr_w = []
    for j in spec_data[0]:
        lyalpha_arr_w.append(Lyalpha)            
    '''

    # the line to mark Ly_alpha
    lya, _ = spectrum.find_nearest(spec_data[0], Lyalpha)
    flx_lya = spectrum.findXinY(spec_data[1], spec_data[0], lya)
    if flx_lya > 1.0:
        y_Lyalpha = [flx_lya+0.1, flx_lya+0.2]
        N = len(y_Lyalpha)
        upy = y_Lyalpha[N-1] + 0.15
    else:
        y_Lyalpha = [1.0, 1.2]
        N = len(y_Lyalpha)
        upy = y_Lyalpha[N-1] + 0.15
    lyalpha_arr_norm = []
    cIII_mie_list = []
    siIII_list = []
    nV_list =[]
    cIII_list = []
    sII_list = []
    for i in range(0, N):
        lyalpha_arr_norm.append(Lyalpha)
        cIII_mie_list.append(cIII_mie)
        siIII_list.append(siIII)
        nV_list.append(nV)
        cIII_list.append(cIII)
        sII_list.append(sII)
    lyalpha_arr = numpy.array([lyalpha_arr_norm, y_Lyalpha])
    
    star_name = base_name.replace(".txt", "")
    # Figure
    # Plot limits
    low = 1160
    up = 1260
    lo_y = -0.05
    up_y = upy
    # making all fonts biger
    font = {#'family' : 'Vera Sans',
            'weight' : 'regular',
            'size'   : 17}
    matplotlib.rc('font', **font)
    
    pyplot.figure(1, figsize=(10, 10))
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Normalized Flux')# [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.ylim(lo_y,up_y)
    #pyplot.suptitle(norm_stars[i])
    pyplot.plot(spec_data[0], spec_data[1], 'b', lyalpha_arr[0], lyalpha_arr[1], 'r--')#, continuum_norm[0], continuum_norm[1], 'magenta')
    # remove the firt tick so that they do not overlap
    pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
    # mark other important lines
    pyplot.plot(cIII_mie_list, lyalpha_arr[1], 'r--', siIII_list, lyalpha_arr[1], 'r--', nV_list, lyalpha_arr[1], 'r--')
    pyplot.plot(cIII_list, lyalpha_arr[1], 'r--')#, sII_list, lyalpha_arr[1], 'r--')
    pyplot.text(cIII_mie-3, y_Lyalpha[N-1]+0.03, 'C III')
    pyplot.text(siIII-3, y_Lyalpha[N-1]+0.03, 'Si III')
    pyplot.text(Lyalpha-4.5, y_Lyalpha[N-1]+0.03, 'Ly-alpha')
    pyplot.text(nV-4, y_Lyalpha[N-1]+0.03, 'N V')
    pyplot.text(cIII-1, y_Lyalpha[N-1]+0.03, 'C III')
    #pyplot.text(sII-1, y_Lyalpha[N-1]+0.03, 'S II')
    #epsfile = os.path.join(path_plots, star_name+".eps")
    #pyplot.savefig(epsfile)    
    pyplot.show()
    '''
    print("Do you want to save this plot?  [y/N , meaning NO is default... :P ]")
    save_plt = raw_input()
    pyplot.ioff()
    if save_plt == 'n' or save_plt == '':
        print('Plot not saved...   :( ')
        os.remove(epsfile)
    else:
        print('Star %s was saved.' % base_name)
    '''


'''
# Printing file of normalized B star and corresponding eqw:
eqwf = open(path_results+'eqwsBstars_v10.txt', 'w+')
print >> eqwf, 'B star', 'Ly-alpha EQW'
print >> eqwf, '*** SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
    print >> eqwf,  normBstars[i], all_eqws[i]
eqwf.close()

# Printing a file with the various eqws
print ('normBstars, all_eqws, half_eqwstimes2, half_eqwstimes2_avoiding_N, VG_eqw, expected_eqwstoVG',
       numpy.shape(normBstars), numpy.shape(eqws_10), numpy.shape(eqws_20), numpy.shape(eqws_30), numpy.shape(half_eqwstimes2_avoiding_N))

comp_eqws = open(path_results+'eqws_COMPARISON3.txt', 'w+')
table_comp = []
print >> comp_eqws, 'SIGN CONVENTION: positive = emission,  negative = absorption'
for i in range(0, len(normBstars)):
#    print >> comp_eqws, normBstars[i], eqws_10[i], eqws_20[i], eqws_30[i], half_eqwstimes2_avoiding_N[i]
    temp = [normBstars[i], repr(eqws_10[i]), repr(eqws_20[i]), repr(eqws_30[i]), 
            repr(half_eqwstimes2_avoiding_N[i]), repr(half_eqwstimes2_avoiding_Si[i])] 
    table_comp.append(temp)
table.pprint_table(comp_eqws, table_comp)
comp_eqws.close()

'''


        