import numpy
import glob
import re
import os
import string
from spectrum import spectrum
from spectrum import pyplot
from claus.clausfile import clausfile


'''
This program calls a CMFGEN file and converts data into two column text files with units of Angstroms and cgs, respectively.
*** Everything is a numpy array.
'''
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

# go into each temperature and get the lines and continuum files into a normal text file of angstroms and cgs
def CMFGENfile_into2cols(path_CMFGEN, path_results, lower_wave, upper_wave, templogg):
    # Find directory
    ogrid_dir = glob.glob(path_CMFGEN)
    print('Now in %s' % ogrid_dir)
    for single_dir in ogrid_dir:
        dir_name = os.path.join(single_dir+'/obs/')
        print'Going into directory: %s' % (dir_name)
        #files_to_use = ['obs_cont', 'obs_fin_10']
        files_to_use = ['obs_cont.fixed', 'obs_fin_10.fixed']
        for i in range(0, len(files_to_use)):
            file_in_use = dir_name + files_to_use[i]
            file_name = os.path.basename(file_in_use)
            print('loading %s' % (file_in_use))
            # Converting file into three columns: frequencies, janskies, error
            clausfile.clausfile(file_in_use, dest=os.path.join(path_results))
            # Converting into Angstroms and cgs to create a text file 
            sp = spectrum.spectrum(os.path.join(path_results+os.path.basename(file_in_use)+"_1d"))
            #sp.save2(file_in_use+"_Acgs.txt")
            output = (sp.A, sp.cgs)
            spectrum.write_1d(os.path.join(os.path.abspath(path_results), file_name+"_Acgs_"+templogg+".txt"), output)
      
  
#### NOW DO THE WORK
# This is where the directory is
path_CMFGEN = '../ogrid_24jun09/'
path_results = '../results/OstarsCMFGEN/'
if not os.path.exists(path_results):
    print("Path: %s does not exist!" % (path_results))
    exit(1)
# Now use those files to plot
lower_wave = 1190.0
upper_wave = 1300.0

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701

# Find directory
ogrid_dir = glob.glob('../ogrid_24jun09/*')

# Generate 1d files of flux and wavelength
templogg_list = []
for i in range(0, len(ogrid_dir)):
    temp_basename = string.split(os.path.basename(ogrid_dir[i]), sep='NT')
    templogg_list.append(temp_basename[1])

''' THIS LOOP IS TO BE RUNNED WHEN NEEDING TO CONVERT CMFGEN FILES INTO 2 COLUMN FILES OF ANGSTROMS AND CGS UNITS. 
    PROBLEM: I could not find a way to change the name of the 1d files of Jy and Hz, so this file is overwritten every time.''' 
'''
for i in range(0, len(templogg_list)):
    print('temperature and log g: %s' % templogg_list[i])
    new_path_CMFGEN = path_CMFGEN+'NT'+templogg_list[i]
    print(new_path_CMFGEN)
    CMFGENfile_into2cols(new_path_CMFGEN, path_results, lower_wave, upper_wave, templogg_list[i])
'''

all_eqws = []
eqw_from_line_all = []
eqw_from_madeup_line_all = []
all_aka = []
# Errors lists
avg_abs_err_to_real_continuum_all = []
avg_relative_err_to_real_continuum_all = []
avg_abs_err_to_real_continuum_all2 = []
avg_relative_err_to_real_continuum_all2 = []
all_chi_squared1 = []
all_chi_squared2 = []


# Generate the lists of the fin and cont files
fin_files = []
cont_files = []
fins = 'fin_10.fixed_Acgs'
conts = 'cont.fixed_Acgs'
results_dir = glob.glob('../results/OstarsCMFGEN/*')
for i in range(0, len(results_dir)):
    if fins in results_dir[i]:
        fin_files.append(results_dir[i])
    elif conts in results_dir[i]:
        cont_files.append(results_dir[i])
        
aka = 0
#a = open('ogrid_aka_dictionary.txt', 'w+')
points_used_for_line_continuums = []

# Load the text files
for f in range(0, len(fin_files)):
    A, cgs = numpy.loadtxt(fin_files[f], dtype=numpy.float64, unpack=True)    
    print('Loading %s' % fin_files[f])
    base_fin = os.path.basename(fin_files[f])
    print('Fin file loaded. Now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    wav, flx = selection(A, cgs, lower_wave, upper_wave)
    lines = numpy.array([wav, flx])
    #print('shape of lines_array', lines.shape)
    all_aka.append(aka)
    #print >> a, aka, fin_files[f]
    '''
    ### Text file of Resolving Power of every wavelength
    data_lines = numpy.array([A, cgs])
    resol = open(path_results+'CMFGENResolvingPower.txt', 'w+')
    resPow_avg = []
    for line in A:
        resPow = spectrum.resolving_power(line, data_lines)
        resPow_avg.append(resPow)
        print >> resol, line, resPow 
        print('wavelength and R', line, resPow)
    avgR = int(sum(resPow_avg)) / len(A)
    print >> resol, 'Average Reolving Power of CMFGEN lines file = '
    print >> resol, avgR
    resol.close()
    print('Average Reolving Power of CMFGEN lines file = %i' % (avgR))
    '''
    
    A_cont, cgs_cont = numpy.loadtxt(cont_files[f], dtype=numpy.float64, unpack=True)   
    base_cont = os.path.basename(cont_files[f])
    print('Loading %s' % cont_files[f]) 
    print('Continuum file loaded. Now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    wav_cont, flx_cont = selection(A_cont, cgs_cont, lower_wave, upper_wave)
    continuum = numpy.array([wav_cont, flx_cont])
    #print('shape of continuum_array', continuum.shape)
        
    # Rebinning to same length as lines array
    # since the lines array is much larger than the continuum array, this last one must be rebinned:
    reference_wavelength = Lyalpha
    # Making the continuum array the same size as the lines array.
    #rebinned_cont, rebinning_factor_cont = spectrum.rebinning_interpol(lines, continuum, reference_wavelength)
    #norm_flx = flx / rebinned_cont
    #norm_lines = numpy.array([wav, norm_flx])
    #norm_cont = spectrum.theo_cont(wav, scale_factor=1.0)
    rebinned_cont, rebinning_factor_cont = spectrum.rebin_one_arr_to_desired_rows(continuum, lines.shape[1])
    lines2, rebinned_cont = spectrum.correct_rebin(lines, rebinned_cont)
    norm_flx_initial_R = lines2[1] / rebinned_cont[1]
    norm_lines_initial_R = numpy.array([lines2[0], norm_flx_initial_R])
    norm_cont_initial_R = spectrum.theo_cont(lines2[0], scale_factor=1.0)

    # Measuring EQW
    # from original shape
    eqw_initial = spectrum.EQW_line_fixed(norm_lines_initial_R, norm_cont_initial_R, Lyalpha, 10.0)
    print ('*** Initial eqw with the EQW "fixed" function = %f' % (eqw_initial))
    
    # Determining Resolving Powers 
    initial_Resolution =  spectrum.resolving_power(Lyalpha, lines2)  # R=20,000  for CMFGEN as indicated in Palacios et al. 2010, A&A
    original_delta_lambda = Lyalpha / float(initial_Resolution)
    print('Initial Resolving Power = %i  ---  Initial delta_lambda = %f' % (initial_Resolution, original_delta_lambda))

    # Now rebinning according to a desired delta_lambda
    '''delta_lambda=0.456 is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)'''
    desired_delta_lambda = 0.465
    lines_rebinned, cont_rebinned = spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, lines, continuum, guessed_rows=500)
    #print('shape of lines_rebinned: %s    ----   shape of cont_rebinned: %s' % (repr(lines_rebinned.shape), repr(cont_rebinned.shape)))

    # Measuring EQW
    # from rebinned shape
    new_eqw = spectrum.EQW_line_fixed(lines_rebinned, cont_rebinned, Lyalpha, 10.0)
    print ('*** Rebinned array eqw_fixed = %f' % (new_eqw))
    all_eqws.append(new_eqw)
    
    # Normalization
    fixed_lines_rebinned, fixed_cont_rebinned = spectrum.correct_rebin(lines_rebinned, cont_rebinned)
    norm_flx = fixed_lines_rebinned[1,:] / fixed_cont_rebinned[1,:]
    
    # Arrays for eqw measurements
    norm_lines = numpy.array([fixed_lines_rebinned[0,:], norm_flx])
    cont_lines = spectrum.theo_cont(fixed_lines_rebinned[0,:], scale_factor=1.0)
    #print('******* norm_lines.shape', norm_lines.shape)

    # Ploting 
    # the line to mark Ly_alpha in rebinned spectra
    lyalpha_arr_reb = []
    for i in lines_rebinned[1]:
        lyalpha_arr_reb.append(Lyalpha)  
    # the line to mark Ly_alpha in rebinned normalized spectra
    lyalpha_arr_wav_rebinned = []
    for i in norm_flx:
        lyalpha_arr_wav_rebinned.append(Lyalpha)  
    # the line to mark Ly_alpha in non-rebinned spectra
    lyalpha_arr_wav = []
    for i in flx:
        lyalpha_arr_wav.append(Lyalpha)  
    # area under eqw
    #eqw_y = selection(lines_rebin[0], norm_flx, lolim, uplim)
    # Lines and Continuum
    # plot limits
    low = 1200
    up = 1230
    '''
    # Figure
    pyplot.figure(1, figsize=(10, 10))
    pyplot.title('AKA Ostar_'+repr(aka))
    pyplot.suptitle('Lines (blue): '+base_fin)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(lines_rebinned[0], lines_rebinned[1], 'b', lyalpha_arr_reb, lines_rebinned[1], 'r--')#, cont_rebinned[0], cont_rebinned[1], 'g')
    pyplot.show()
    
    pyplot.figure(3, figsize=(10, 10))
    pyplot.title('AKA Ostar_'+repr(aka))
    pyplot.suptitle('Normalization: '+fin_files[f]+' / '+cont_files[f])
    pyplot.xlim(low, up)
    #pyplot.ylim(0.4, 1.5)
    pyplot.plot(lines_rebinned[0], norm_flx, 'b', lyalpha_arr_wav_rebinned, norm_flx, 'r--')
    #pyplot.fill_between(eqw, eqw_y, 0,color='g')
    
    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('AKA Ostar_'+repr(aka))
    pyplot.suptitle('blue: '+base_fin+'  green: '+base_cont)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1000, 2000)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(A, cgs, 'b', lyalpha_arr_wav, flx, 'r--', wav_cont, flx_cont, 'g')
    '''
    # My continuum function, two-point linear equation: y = m*(x -x1) + y1 
    print('Enter the lower and upper X-axis to consider in the line:')
    x1 = 1274.63 #float(raw_input("Low limit: "))
    x2 = 1290.42 #float(raw_input("Upper limit: "))
    temp_x1, _ = spectrum.find_nearest(lines_rebinned[0], x1)
    temp_x2, _ = spectrum.find_nearest(lines_rebinned[0], x2)
    print('closest wavelengths: %f, %f' % (temp_x1, temp_x2))
    y1 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x1)
    y2 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x2)
    temp_obj = string.split(base_fin, sep='_')
    print (temp_obj[4])
    if (numpy.fabs(y2) > numpy.fabs(y1)) or (float(temp_obj[4]) > 35000.0) and (float(temp_obj[4]) < 48000.0):
        x1 = 1221.76
        x2 = 1298.58
        temp_x1, _ = spectrum.find_nearest(lines_rebinned[0], x1)
        temp_x2, _ = spectrum.find_nearest(lines_rebinned[0], x2)
        y1 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x1)
        y2 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x2)
    if (float(temp_obj[4]) > 48000.0):
        x1 = 1202.95
        x2 = 1298.58
        temp_x1, _ = spectrum.find_nearest(lines_rebinned[0], x1)
        temp_x2, _ = spectrum.find_nearest(lines_rebinned[0], x2)
        y1 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x1)
        y2 = spectrum.findXinY(lines_rebinned[1], lines_rebinned[0], temp_x2)
        
    m = (y2 - y1) / (x2 - x1)
    y_list = []
    for x in lines_rebinned[0]:
        y = m * (x - x1) + y1
        y_list.append(y)
    my_cont_arr = numpy.array([lines_rebinned[0], y_list])
    
    # My continuum without the flux array: : y = m*(x -x1) + y1 
    print('Enter the lower and upper X-axis to consider in the line:')
    tx1 = x1#float(raw_input("x1 pont: "))
    ty1 = y1+(y1/100.0)#float(raw_input("y1 pont: "))
    tx2 = x2#float(raw_input("x2 pont: "))
    ty2 = y2+(y2/100.0)#float(raw_input("y2 pont: "))
    tm = (ty2 - ty1) / (tx2 - tx1)
    ty_list = []
    for x in lines_rebinned[0]:
        ty = tm * (x - tx1) + ty1
        ty_list.append(ty)
    tmy_cont_arr = numpy.array([lines_rebinned[0], ty_list])
    
    '''        
    pyplot.figure(3, figsize=(10, 10))
    pyplot.title('AKA Ostar_'+repr(aka))
    pyplot.suptitle('My continuum fit to: '+base_fin+' and: '+base_cont)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(lines_rebinned[0], lines_rebinned[1], 'b', lyalpha_arr_reb, lines_rebinned[1], 'r--', 
                cont_rebinned[0], cont_rebinned[1], 'g', my_cont_arr[0], my_cont_arr[1], 'magenta',
                tmy_cont_arr[0], tmy_cont_arr[1], 'cyan')
    '''
    
    # Normalization to my continuum: line from array
    new_norm_flx = lines_rebinned[1] / my_cont_arr[1]
    new_norm_lines = numpy.array([lines_rebinned[0], new_norm_flx])
    
    # Normalization to my continuum: made-up line
    tnew_norm_flx = lines_rebinned[1] / tmy_cont_arr[1]
    tnew_norm_lines = numpy.array([lines_rebinned[0], tnew_norm_flx])
    
    # ERRORS
    mean = sum(my_cont_arr[1]) / len(my_cont_arr[1])
    var_list = []
    for i in range(0, len(my_cont_arr[1])):
        var = (my_cont_arr[1,i] - mean)**2
        var_list.append(var)
    variance =  sum(var_list) / len(my_cont_arr[1])
    chi2 = []
    for i in range(0, len(my_cont_arr[1])):
        diff_squared = ((my_cont_arr[1, i] - cont_rebinned[1, i])**2)
        chi2.append(diff_squared) 
    chi_squared = sum(chi2) / variance
    
    # Difference of my continuum from the real continuum
    relative_err_to_real_continuum_list = []
    abs_err_to_real_continuum_list = []
    for i in range(0, len(my_cont_arr[0])):
        relative_err_to_real_continuum = ((my_cont_arr[1,i] / cont_rebinned[1,i]) - 1.0) * 100
        relative_err_to_real_continuum_list.append(relative_err_to_real_continuum)
        abs_err_to_real_continuum = my_cont_arr[1,i] - cont_rebinned[1,i]
        abs_err_to_real_continuum_list.append(abs_err_to_real_continuum)
        
    avg_relative_err_to_real_continuum = sum(relative_err_to_real_continuum_list) / len(relative_err_to_real_continuum_list)
    avg_abs_err_to_real_continuum = sum(abs_err_to_real_continuum_list) / len(abs_err_to_real_continuum_list)

    # ERRORS for made-up line
    tmean = sum(tmy_cont_arr[1]) / len(tmy_cont_arr[1])
    tvar_list = []
    for i in range(0, len(tmy_cont_arr[1])):
        tvar = (tmy_cont_arr[1,i] - tmean)**2
        tvar_list.append(tvar)
    tvariance =  sum(tvar_list) / len(tmy_cont_arr[1])
    tchi2 = []
    for i in range(0, len(tmy_cont_arr[1])):
        tdiff_squared = ((tmy_cont_arr[1, i] - cont_rebinned[1, i])**2)
        tchi2.append(tdiff_squared) 
    tchi_squared = sum(tchi2) / tvariance
    
    # Difference of my made-up line continuum to the real continuum
    trelative_err_to_real_continuum_list = []
    tabs_err_to_real_continuum_list = []
    for i in range(0, len(tmy_cont_arr[0])):
        trelative_err_to_real_continuum = ((tmy_cont_arr[1,i] / cont_rebinned[1,i]) - 1.0) * 100
        trelative_err_to_real_continuum_list.append(trelative_err_to_real_continuum)
        tabs_err_to_real_continuum = tmy_cont_arr[1,i] - cont_rebinned[1,i]
        tabs_err_to_real_continuum_list.append(tabs_err_to_real_continuum)
        
    tavg_relative_err_to_real_continuum = sum(trelative_err_to_real_continuum_list) / len(trelative_err_to_real_continuum_list)
    tavg_abs_err_to_real_continuum = sum(tabs_err_to_real_continuum_list) / len(tabs_err_to_real_continuum_list)
    
    avg_abs_err_to_real_continuum_all.append(avg_abs_err_to_real_continuum)
    avg_relative_err_to_real_continuum_all.append(avg_relative_err_to_real_continuum)
    avg_abs_err_to_real_continuum_all2.append(tavg_abs_err_to_real_continuum)
    avg_relative_err_to_real_continuum_all2.append(tavg_relative_err_to_real_continuum)
    
    print('Average Absolute error of my normalization to the real one: %f' % (avg_abs_err_to_real_continuum))
    print('Average Absolute error of the other normalization to the real one: %f' % (tavg_abs_err_to_real_continuum))
    print('Average Relative error of my line to the continuum: %f' % (avg_relative_err_to_real_continuum))
    print('Average Relative error of made-up line to the continuum: %f' % (tavg_relative_err_to_real_continuum))
    print('Chi-squared of my fit = %f' % chi_squared)
    all_chi_squared1.append(chi_squared)
    print('Chi-squared of my fit with made up line = %f' % tchi_squared)
    all_chi_squared2.append(tchi_squared)
    
    coords = [x1, x2, y1, y2, ty1, ty2]
    points_used_for_line_continuums.append(coords)
    
    ''' 
    pyplot.figure(4, figsize=(10, 10))
    pyplot.title('AKA Ostar_'+repr(aka))
    pyplot.suptitle('Normalization: '+base_fin+' / '+base_cont)
    pyplot.xlim(low, up)
    #pyplot.ylim(0.4, 1.5)
    pyplot.plot(norm_lines[0], norm_lines[1], 'b', lyalpha_arr_wav_rebinned, norm_lines[1], 'r--',
                new_norm_lines[0], new_norm_lines[1], 'magenta', tnew_norm_lines[0], tnew_norm_lines[1], 'c--')
    #pyplot.fill_between(eqw, eqw_y, 0,color='g')
    pyplot.show()
    '''
    
    cont_eqw_from_line = spectrum.theo_cont(new_norm_lines[0])
    eqw_from_line = spectrum.EQW_line_fixed(new_norm_lines, cont_eqw_from_line, Lyalpha, 10.0)
    cont_eqw_from_madeup_line = spectrum.theo_cont(tnew_norm_lines[0])
    eqw_from_madeup_line = spectrum.EQW_line_fixed(tnew_norm_lines, cont_eqw_from_madeup_line, Lyalpha, 10.0)
    print ('*** My continuum line eqw_fixed = %f' % (eqw_from_line))
    print ('*** Made-up continuum line eqw_fixed = %f' % (eqw_from_madeup_line))
    eqw_from_line_all.append(eqw_from_line)
    eqw_from_madeup_line_all.append(eqw_from_madeup_line)
    
    aka = aka+1
    
    
#a.close()    
print('Done!')
#print numpy.shape(all_aka), numpy.shape(avg_abs_err_to_real_continuum_all), numpy.shape(avg_abs_err_to_real_continuum_all2)
#print numpy.shape(avg_relative_err_to_real_continuum_all), numpy.shape(avg_relative_err_to_real_continuum_all2)
#print numpy.shape(all_eqws), numpy.shape(eqw_from_line_all), numpy.shape(eqw_from_madeup_line_all), 
#print numpy.shape(all_chi_squared1), numpy.shape(all_chi_squared2), numpy.shape(points_used_for_line_continuums)


# EQW file with my two fitted lines
f = open('CMFGEN_continuum_tests.txt', 'w+')
print >> f, 'AKA, Absolute error between the lines as continuum --- Relative error between lines and continuum  --- EQWs: real, from line, from made-up line'
print >> f, 'Ostar_# - Line from flux array, Made up line - Line from flux array, Made up line - EQW[A]: real, from line, from made-up line - Chi^2 of my fit - Chi^2 made-up line'  
for i in range(0, len(all_aka)):
    print >> f, all_aka[i], avg_abs_err_to_real_continuum_all[i], avg_abs_err_to_real_continuum_all2[i], 
    avg_relative_err_to_real_continuum_all[i], avg_relative_err_to_real_continuum_all2[i], 
    all_eqws, eqw_from_line_all[i], eqw_from_madeup_line_all[i], 
    all_chi_squared1[i], all_chi_squared2[i], 
    points_used_for_line_continuums[i][0],
    points_used_for_line_continuums[i][1],
    points_used_for_line_continuums[i][2],
    points_used_for_line_continuums[i][3],
    points_used_for_line_continuums[i][4],
    points_used_for_line_continuums[i][5]
f.close()

# File with coords per object:
f = open('Ostar_coords.txt', 'w+')
print >> f, 'Ostar_#, Wavelengths, fluxes from arr, made-up fluxes'
print >> f, '          x1,   x2,      y1,   y2,        ty1,   ty2'
for i in range(0, len(all_aka)):
    print >> f, all_aka[i], points_used_for_line_continuums[i]
f.close()

