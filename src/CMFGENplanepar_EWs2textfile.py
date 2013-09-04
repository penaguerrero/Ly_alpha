import numpy
import os
import string
import science
from matplotlib import pyplot


'''
This program calls the plane-parallel CMFGEN files, converts data into two column text files with units of Angstroms and cgs, respectively, and measures EWs
storing them in a text file.

It can work with ALL stars in the file or with just the selected few.
'''

##################################################################################################################################

''' These are all the switches than can be turned on or off in this code: '''
# This is the machine where I am working
#machine = 'gallo'
machine = 'pena'

# create the wavelength and intensity [cgs] text files from original CMFGEN files
make_text_files = False

# use all stars from folder
all_stars = True

# if all_stars is False this list will be used
stars_list = ['T16000logg300_OBSFLUX_1d_Acgs.txt', 'T22000logg375_OBSFLUX_1d_Acgs.txt', 'T25000logg375_OBSFLUX_1d_Acgs.txt']

# determine the average resolving power for the entire spectrum, 90-10,500 Angstroms
determine_avgR = False      

# To test the linear continuum fitting comparing the one obtained before rebinning versus after rebinning
compare_to_find_continuum_then_rebin = False

# These are the plots of either all stars or just the selected ones
want_to_see_plots = True

# Create the text file with the temperature, log g values, and equivalent widths.
# There are 2 EW measuremets: simple-EW = full integration set to 40A, and half-EW = integration from LyAlpha+20A * 2
create_textfile_of_EWs = True

##################################################################################################################################
  
#    NOW DO THE WORK

# This is where the directory is
path_CMFGEN = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/plane_par_mod_26apr11/'
path_results = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/'
dataCMFGEN = science.utility.DataDir(path_CMFGEN)
data_results = science.utility.DataDir(path_results)
if not os.path.exists(data_results.path):
    print("Path: %s does not exist!" % (data_results.path))
    exit(1)
    
# Now use those files to plot
lower_wave = 1100.0
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
if make_text_files == True:
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

# Choose if you want to do all stars or list the ones you want
if all_stars == True:
    stars = data_results.contents()
else:
    stars = stars_list

# there might be other files in the results directory, so search for the right files to measure EWs
for item in stars:
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

# List to store the EWs
simple_eqw_list = []
half_eqw_list = []

# Load the text files
for i in range(len(fin_files)):
    A, cgs = numpy.loadtxt(os.path.join(data_results.path,fin_files[i]), dtype=numpy.float64, unpack=True)
    print('Fin file loaded. Here is the temperature and log g:  %i  %0.2f' % (temps[i], loggs[i]))
    # Determine the average resolving Power of the whole spectrum
    # CAUTION: this takes a few minutes, skip this loop if wanting faster run time,
    # make     determine_avgR = False
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
    
    # Rebin the data array and then calculate a continuum
    rebinned_selec_data, smoothing_R_factor = science.spectrum.rebin_one_arr_to_desired_resolution(desired_delta_lambda, Lyalpha, selected_data, guessed_rows=500)
    #print ('smoothing_R_factor = %f' % smoothing_R_factor)
    continuum_reb_selec_data = science.CMFGEN.find_linear_continuum(rebinned_selec_data, temps[i])        
    # plot limits
    low = 1160
    up = 1260
    pyplot.figure(1, figsize=(10, 7))
    
    '''
    If wanting to comare the difference between determining the continuum before of after the rebinning turm
    compare_to_find_continuum_then_rebin to TRUE.
    '''
    if compare_to_find_continuum_then_rebin == True:
        # This is to calculate the continuum before rebinning and then rebin data and continuum
        continuum_selected_data2 = science.CMFGEN.find_linear_continuum(selected_data, temps[i])    
        rebinned_selec_data2, rebinned_continuum_data = science.spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, selected_data, continuum_selected_data2)
        # Plot everything to compare
        pyplot.title('Raw and Rebinned data')
        pyplot.suptitle(temps[i])
        pyplot.xlabel('Wavelength [$\AA$]')
        pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
        pyplot.xlim(low, up)
        pyplot.plot(selected_data[0], selected_data[1], 'b',    # the raw selected data
                    continuum_selected_data2[0], continuum_selected_data2[1], 'm:', # the continuum determined from the raw selected data
                    rebinned_selec_data[0], rebinned_selec_data[1], 'k',    # the rebinned selected data
                    continuum_reb_selec_data[0], continuum_reb_selec_data[1], 'r--',  # the continuum determined from the rebinned selected data
                    rebinned_continuum_data[0], rebinned_continuum_data[1], 'g:') # the rebinned coninuum determined from the raw selected data
        pyplot.show()
    
    # Divide by the continuum 
    norm_rebinned_flx = rebinned_selec_data[1] / continuum_reb_selec_data[1]
    norm_data = numpy.array([rebinned_selec_data[0], norm_rebinned_flx])
    # determine theoretical continuum
    continuum_theo = science.spectrum.theo_cont(norm_data[0], scale_factor=1.0)
    
    '''To see plots turn want_to_see_plots to TRUE'''
    if want_to_see_plots == True:
        # plot the normalized data
        pyplot.title('Teff = '+repr(temps[i])+'   log g = '+repr(loggs[i]))
        pyplot.suptitle('Rebinned and normalized data')
        pyplot.xlabel('Wavelength [$\AA$]')
        pyplot.ylabel('Normalized Flux')
        pyplot.xlim(low, up)
        pyplot.plot(norm_data[0], norm_data[1], 'k')
        pyplot.show()
    
    # Determine EWs
    simple_eqw = science.spectrum.EQW_line_fixed(norm_data, continuum_theo, Lyalpha, width=40.)
    print('Ly-alpha simple EQW = %f' % simple_eqw)
    simple_eqw_list.append(simple_eqw)
    
    eqw_limit = Lyalpha + 20.0
    half_eqw = science.spectrum.half_EQW_times2(norm_data, continuum_theo, Lyalpha, eqw_limit)
    print('Ly-alpha half-EQW = %f' % half_eqw)
    half_eqw_list.append(half_eqw)
    
# Print a file with the all EWs
if create_textfile_of_EWs == True:
    f = open(path_results+'planepar_EWs.txt', 'w+')
    print >> f, 'Teff    log g    simple-EW    half-EW'
    for i in range(len(temps)):
        s = ('%i    %0.2f    %0.3f    %0.3f\n' % (temps[i], loggs[i], simple_eqw_list[i], half_eqw_list[i]))
        f.write(s)
    f.close()

print 'Done!'

