import numpy
import glob
#import re
import os
import shutil
import string
#import table
import science
from matplotlib import pyplot


'''
This program calls the plane-parallel CMFGEN files, converts data into two column text files with units of Angstroms and cgs, respectively, and measures EWs
storing them in a text file.

*** Everything is a numpy array.
'''
  
#### NOW DO THE WORK
# This is where the directory is
machine = 'gallo'
path_CMFGEN = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/plane_par_mod_26apr11/'
path_results = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/'
dataCMFGEN = science.utility.DataDir(path_CMFGEN)
data_results = science.utility.DataDir(path_results)
if not os.path.exists(data_results.path):
    print("Path: %s does not exist!" % (data_results.path))
    exit(1)
    
# Now use those files to plot
lower_wave = 1190.0
upper_wave = 1300.0

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701

# Find all the stuff in directory
ogrid_dir = dataCMFGEN.contents()

# Generate 1d files of flux and wavelength
templogg_list = []

for i in range(0, len(ogrid_dir)):
    temp_basename = string.split(os.path.basename(ogrid_dir[i]), sep='T')
    # equivalently, I could write the following for python 3
    #thebase = os.path.basename(ogrid_dir[i]).split('T')
    templogg_list.append(temp_basename[1])

'''
 THIS LOOP IS TO BE RUNNED WHEN NEEDING TO CONVERT CMFGEN FILES INTO 2 COLUMN FILES OF ANGSTROMS AND CGS UNITS. 
    ****   PROBLEM: I could not find a way to change the name of the 1d files of Jy and Hz, so this file is overwritten every time.
'''
get_text_files = False

if get_text_files == True:
    for i in range(0, len(templogg_list)):
        print('temperature and log g: %s' % templogg_list[i])
        new_path_CMFGEN = path_CMFGEN+'T'+templogg_list[i]
        print(new_path_CMFGEN)
        science.CMFGEN.into2cols(new_path_CMFGEN, path_results, lower_wave, upper_wave, templogg_list[i])


# Generate the lists of the fin and cont files
temps = []
loggs = []
fin_files = []
converted = '_Acgs.txt'
# there might be other files in the results directory, so search for the right files to measure EWs
for item in data_results.contents():
    if converted in item:
        fin_files.append(item)        
        kk1 = string.split(item, sep='_')
        kk2 = kk1[0].replace("T", "")
        kk3 = string.split(kk2, sep='logg')
        teff = int(kk3[0])
        logg = float(kk3[1])*0.01
        temps.append(teff)
        loggs.append(logg)
        
# And now create the continuums
cont_files = []
points_used_for_line_continuums = []

# Load the text files
for i in range(len(fin_files)):
    A, cgs = numpy.loadtxt(os.path.join(data_results.path,fin_files[i]), dtype=numpy.float64, unpack=True)
    print('Fin file loaded. Here is the temperature and log g:  %i  %0.2f' % (temps[i], loggs[i]))
    # Determine the average resolving Power of the whole spectrum
    # CAUTION: this takes a few minutes, skip this loop if wanting faster run time,
    # make     determine_avgR = False
    determine_avgR = False
    if determine_avgR == True:
        full_data = numpy.array([A, cgs])
        R_list = []
        for line in A:
            R = science.spectrum.resolving_power(line, full_data)
            R_list.append(R)
        avgR = int(sum(R_list)) / len(A)
        print('Average Resolving power = %i' % avgR)
    
    # Selecting the range for Ly_alpha
    wav, flx = science.spectrum.selection(A, cgs, lower_wave, upper_wave)
    selected_data = numpy.array([wav, flx])
    
    # Determine resolving power at LyAlpha and delta_lambda
    R_lyalpha = science.spectrum.resolving_power(Lyalpha, selected_data)
    original_delta_lambda = Lyalpha / float(R_lyalpha)
    #print('Initial Resolving Power = %i  ---  Initial delta_lambda = %f' % (R_lyalpha, original_delta_lambda))
    
    # Rebin to desired delta_lambda and calculate the local continuum as a linear approximation
    '''delta_lambda=0.456 is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)'''
    desired_delta_lambda = 0.465
    
    # This is to calculate the continuum before rebinning and then rebin data and continuum
    #continuum_selected_data2 = science.CMFGEN.find_linear_continuum(selected_data, temps[i])    
    #rebinned_selec_data2, rebinned_continuum_data = science.spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, selected_data, continuum_selected_data2)
    # equivalently, I can simply rebin the data array and then calculate a continuum
    rebinned_selec_data, smoothing_R_factor = science.spectrum.rebin_one_arr_to_desired_resolution(desired_delta_lambda, Lyalpha, selected_data, guessed_rows=500)
    #print ('smoothing_R_factor = %f' % smoothing_R_factor)
    rebinned_continuum_data = science.CMFGEN.find_linear_continuum(rebinned_selec_data, temps[i])        
    # plot limits
    low = 1190
    up = 1300
    pyplot.figure(1, figsize=(10, 7))
    '''
    # Plot everything to compare
    pyplot.title('Raw and Rebinned data')
    pyplot.suptitle(temps[i])
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.plot(selected_data[0], selected_data[1], 'b', 
                rebinned_selec_data[0], rebinned_selec_data[1], 'k',
                continuum_selected_data[0], continuum_selected_data[1], 'm:',
                rebinned_continuum_data[0], rebinned_continuum_data[1], 'r--',
                rebinned_continuum_data2[0], rebinned_continuum_data2[1], 'g:')
    pyplot.show()
    '''
    # Divide by the continuum 
    norm_rebinned_flx = rebinned_selec_data[1] / rebinned_continuum_data[1]
    norm_data = numpy.array([rebinned_selec_data[0], norm_rebinned_flx])
    # plot
    pyplot.title('Teff = '+repr(temps[i]))
    pyplot.suptitle('Rebinned and normalized data')
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Normalized Flux')
    pyplot.xlim(low, up)
    pyplot.plot(norm_data[0], norm_data[1], 'k')
    pyplot.show()
    

    
