import os
import glob
import numpy
import string
from spectrum import spectrum
from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator
from pylab import *


'''
This program simply produces nice plots of the selected stars. 

### It reads the txt files creaded with running Bstars6.py (the nice table with the various EQWs is in Bstars8.py)

'''

# MEASUREMENTS
# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.3 #1215.6701
# rest wavelengths from Leitherer et al(2011), ApJ, 141, 37
# the numbers are different because they look prettier in the plot
cIII_mie = 1175.53
siIII = 1206.50
nV = 1242.80
cIII = 1247.38
sII = 1250.58

  
half_eqwstimes2_avoiding_Si = []
half_eqwstimes2_avoiding_Si.append('(LyA+20)*2')

path_results = '../results/BTlustyStars/'
path_plots = '../results/plots/'

# list of the stars to be ploted
norm_stars = [path_results+'Norm_Bstar_0_v10_15000g175.txt', path_results+'Norm_Bstar_38_v2_18000g225.txt']
# OR when I want to see the plot of ALL stars use the following line:
#norm_stars = glob.glob(path_results+'Norm_Bstar*')

for i in range(0, len(norm_stars)):
    print('Working with: %s' % (norm_stars[i]))
    spec_data = numpy.loadtxt(norm_stars[i], dtype=float, unpack=True)
    continuum_norm = spectrum.theo_cont(spec_data[0], scale_factor=1.0)
        
    base_name = os.path.basename(norm_stars[i])
    base_name2 = string.split(base_name, sep='_')
    Teff, base_name3 = string.split(base_name2[4], sep='g')
    base_name4 = base_name3.replace(".txt", "")
    logg = base_name4.replace("g", "")
    print(Teff, logg)
       
    half_eqw_times2_avoiding_Si_obj = spectrum.half_EQW_times2(spec_data, continuum_norm, Lyalpha, wave_limit=(Lyalpha+20.0), right_side=True)
    print('half_eqwstimes2_avoiding_Si_obj = %f' % half_eqw_times2_avoiding_Si_obj)
    half_eqwstimes2_avoiding_Si.append(half_eqw_times2_avoiding_Si_obj)
    print ' '
    
    # the line to mark the lines
    lya, _ = spectrum.find_nearest(spec_data[0], Lyalpha)
    flx_lya = spectrum.findXinY(spec_data[1], spec_data[0], lya)
    if flx_lya > 1.0:
        y_Lyalpha = [flx_lya+0.1, flx_lya+0.2]
        N = len(y_Lyalpha)
        upy = y_Lyalpha[N-1] + 0.1
    else:
        y_Lyalpha = [1.0, 1.2]
        N = len(y_Lyalpha)
        upy = y_Lyalpha[N-1] + 0.1
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
    font = {#'family' : 'regular',
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
    pyplot.text(siIII-4.3, y_Lyalpha[N-1]+0.03, 'Si III')
    pyplot.text(Lyalpha-4.5, y_Lyalpha[N-1]+0.03, 'Ly-alpha')
    pyplot.text(nV-4, y_Lyalpha[N-1]+0.03, 'N V')
    pyplot.text(cIII-0.5, y_Lyalpha[N-1]+0.03, 'C III')
    #pyplot.text(sII-1, y_Lyalpha[N-1]+0.03, 'S II')
    epsfile = os.path.join(path_plots, star_name+".eps")
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
    
print "Done!"

        