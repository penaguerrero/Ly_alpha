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
differences1_list = []
differences2_list = []
rel_error1 = []
rel_error2 = []
all_eqws = []

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

# Load the text files
for f in range(0, len(fin_files)):
    A, cgs = numpy.loadtxt(fin_files[f], dtype=numpy.float64, unpack=True)    
    print('Loading %s' % fin_files[f])
    print('Fin file loaded. Now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    wav, flx = selection(A, cgs, lower_wave, upper_wave)
    lines = numpy.array([wav, flx])
    print('shape of lines_array', lines.shape)
    
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
    print('Loading %s' % cont_files[f]) 
    print('Continuum file loaded. Now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    wav_cont, flx_cont = selection(A_cont, cgs_cont, lower_wave, upper_wave)
    continuum = numpy.array([wav_cont, flx_cont])
    print('shape of continuum_array', continuum.shape)
        
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
    norm_flx = lines2[1] / rebinned_cont[1]
    norm_lines = numpy.array([lines2[0], norm_flx])
    norm_cont = spectrum.theo_cont(lines2[0], scale_factor=1.0)

    # Measuring EQW
    # from original shape
    eqw_initial = spectrum.EQW_line_fixed(norm_lines, norm_cont, Lyalpha, 10.0)
    print ('Initial eqw with the EQW "fixed" function = %f' % (eqw_initial))
    
    # Determining Resolving Powers 
    initial_Resolution =  spectrum.resolving_power(Lyalpha, lines)  # R=20,000  for CMFGEN as indicated in Palacios et al. 2010, A&A
    original_delta_lambda = Lyalpha / float(initial_Resolution)
    print('Initial Resolving Power = %i  ---  Initial delta_lambda = %f' % (initial_Resolution, original_delta_lambda))

    # Now rebinning according to a desired delta_lambda
    '''delta_lambda=0.456 is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)'''
    desired_delta_lambda = 0.465
    lines_rebinned, cont_rebinned = spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, lines, continuum, guessed_rows=500)
    print('shape of lines_rebinned: %s    ----   shape of cont_rebinned: %s' % (repr(lines_rebinned.shape), repr(cont_rebinned.shape)))

    # Measuring EQW
    # from rebinned shape
    new_eqw = spectrum.EQW_line_fixed(lines_rebinned, cont_rebinned, Lyalpha, 10.0)
    print ('Rebinned array eqw_fixed = %f' % (new_eqw))
    all_eqws.append(new_eqw)
    
    # Normalization
    fixed_lines_rebinned, fixed_cont_rebinned = spectrum.correct_rebin(lines_rebinned, cont_rebinned)
    norm_flx = fixed_lines_rebinned[1,:] / fixed_cont_rebinned[1,:]
    
    # Arrays for eqw measurements
    norm_lines = numpy.array([lines_rebinned[0,:], norm_flx])
    cont_lines = spectrum.theo_cont(lines_rebinned[0,:], scale_factor=1.0)

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
    # Figure
    pyplot.figure(1, figsize=(10, 10))
    pyplot.title('Lines (blue) and Continuum (green)')
    pyplot.suptitle(fin_files[f])
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(lines_rebinned[0], lines_rebinned[1], 'b', lyalpha_arr_reb, lines_rebinned[1], 'r--')#, cont_rebinned[0], cont_rebinned[1], 'g')
    
    '''
    pyplot.figure(3, figsize=(10, 10))
    pyplot.title('Normalization')
    pyplot.suptitle(fin_files[f])
    pyplot.xlim(low, up)
    #pyplot.ylim(0.4, 1.5)
    pyplot.plot(lines_rebinned[0], norm_flx, 'b', lyalpha_arr_wav_rebinned, norm_flx, 'r--')
    #pyplot.fill_between(eqw, eqw_y, 0,color='g')

    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('Lines (blue) and Continuum (green)')
    pyplot.suptitle(fin_files[f])
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(wav, flx, 'b', lyalpha_arr_wav, flx, 'r--', wav_cont, flx_cont, 'g')
    '''
    pyplot.show()
    
    # My continuum function, two-point linear equation: y = m*(x -x1) + y1 
    print('Enter the lower and upper X-axis to consider in the line:')
    x1 = float(raw_input("Low limit: "))
    x2 = float(raw_input("Upper limit: "))
    temp_x1, _ = spectrum.find_nearest(lines_rebinned[0], x1)
    temp_x2, _ = spectrum.find_nearest(lines_rebinned[0], x2)
    print('closest wavelengths: %f, %f' % (temp_x1, temp_x2))
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
    ty1 = float(raw_input("y1 pont: "))
    tx2 = x2#float(raw_input("x2 pont: "))
    ty2 = float(raw_input("y2 pont: "))
    tm = (ty2 - ty1) / (tx2 - tx1)
    ty_list = []
    for x in lines_rebinned[0]:
        ty = tm * (x - tx1) + ty1
        ty_list.append(ty)
    tmy_cont_arr = numpy.array([lines_rebinned[0], ty_list])
        
    # Normalization to my continuum
    new_norm_flx = lines_rebinned[1] / my_cont_arr[1]
    
    pyplot.figure(1, figsize=(10, 10))
    pyplot.title('My continuum fit')
    pyplot.suptitle(fin_files[f])
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    #pyplot.ylim(0, 3e-8)
    pyplot.plot(lines_rebinned[0], lines_rebinned[1], 'b', lyalpha_arr_reb, lines_rebinned[1], 'r--', 
                cont_rebinned[0], cont_rebinned[1], 'g', lines_rebinned[0], y_list, 'magenta',
                lines_rebinned[0], ty_list, 'cyan')

    # Normalization to my continuum
    tnew_norm_flx = lines_rebinned[1] / tmy_cont_arr[1]
    
    # Relative Error
    for i in range(0, len(new_norm_flx)):
        absolute_error1 = abs(new_norm_flx[i] - norm_flx[i])
        relative_error1 = (absolute_error1) / norm_flx[i] * 100
        differences1_list.append(absolute_error1)
        rel_error1.append(relative_error1)
    for i in range(0, len(tnew_norm_flx)):        
        absolute_error2 = abs(tnew_norm_flx[i] - norm_flx[i])
        relative_error2 = (absolute_error2) / norm_flx[i] * 100
        differences2_list.append(absolute_error2)
        rel_error2.append(relative_error2)
    diff_with_mynorm_avg = (sum(differences1_list)/len(differences1_list))**0.5
    tdiff_with_mynorm_avg = (sum(differences2_list)/len(differences2_list))**0.5
    rel_error_avg1 = (sum(rel_error1)/len(rel_error1))**0.5
    rel_error_avg2 = (sum(rel_error2)/len(rel_error2))**0.5
    print('Average Absolute error of my normalization to the real one: %f' % (diff_with_mynorm_avg))
    print('Average Absolute error of the other normalization to the real one: %f' % (tdiff_with_mynorm_avg))
    print('Average Relative error of my normalization to the real one: %f' % (rel_error_avg1))
    print('Average Relative error of the other normalization to the real one: %f' % (rel_error_avg2))
    # Chi squared
    chi_squared1 = []
    chi_squared2 = []
    for re1 in rel_error1:
        temp_chi_squared1 = re1*re1
        chi_squared1.append(temp_chi_squared1)
    Avg_chi_squared1 = sum(chi_squared1)/len(chi_squared1)
    for re2 in rel_error2:
        temp_chi_squared2 = re2*re2
        chi_squared2.append(temp_chi_squared2)
    Avg_chi_squared2 = sum(chi_squared2)/len(chi_squared2)
    print('Average Chi squared of my fit = %f' % Avg_chi_squared1)
    print('Average Chi squared of my fitwith made up line = %f' % Avg_chi_squared2)

    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('Normalization')
    pyplot.suptitle(fin_files[f])
    pyplot.xlim(low, up)
    #pyplot.ylim(0.4, 1.5)
    pyplot.plot(lines_rebinned[0], norm_flx, 'b', lyalpha_arr_wav_rebinned, norm_flx, 'r--',
                my_cont_arr[0], new_norm_flx, 'magenta', tmy_cont_arr[0], tnew_norm_flx, 'c--')
    #pyplot.fill_between(eqw, eqw_y, 0,color='g')

    pyplot.show()

print('Absolute error between the lines as continuum --- Relative error between lines and continuum  --- EQWs')
print('Line from flux array, Made up line  ---  Line from flux array, Made up line  --->  EQW[A]')
print(differences1_list, differences2_list, rel_error1, rel_error2, all_eqws)




