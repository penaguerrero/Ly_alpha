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
            # This is my way
            #print'Currently reading: %s' % (file_name)
            #print("Converting...")
            #name_of_outfile = 'Ostar_'+single_dir+file_name+'.txt'
            #clausfile.convert.generate_file(path_results, name_of_outfile)   
            #print'Here is the 1D nice file: %s' % (name_of_outfile)  
            
            # Joe's way
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
#txt_files = [path_results+'*cont_1d', path_results+'*fin10_1d']
cont_files = glob.glob(os.path.join(os.path.abspath(path_results),'*cont*Acgs*'))
fin_files = glob.glob(os.path.join(os.path.abspath(path_results),'*_fin_10*Acgs*'))

# Load the text files
for f in range(0, len(fin_files)):
    A, cgs = numpy.loadtxt(fin_files[f], dtype=numpy.float64, unpack=True)    
    print('Lines file loaded, now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    #wav = A
    #flx = cgs    
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
    print('Continuum file loaded, now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    #wav_cont = A_cont
    #flx_cont = cgs_cont
    wav_cont, flx_cont = selection(A_cont, cgs_cont, lower_wave, upper_wave)
    continuum = numpy.array([wav_cont, flx_cont])
    print('shape of continuum_array', continuum.shape)
    
    # Rebinning
    lines_rebin, cont_rebin, new_cont_factor, new_lines_factor = spectrum.get_factors_and_rebin(lines, continuum, 250)
    Resolution =  spectrum.resolving_power(Lyalpha, lines)  # R=20,000  for CMFGEN as indicated in Palacios et al. 2010, A&A
    original_delta_lambda = Lyalpha / float(Resolution)
    closest_lyA, idxLyA = spectrum.find_nearest(lines_rebin[0,:], Lyalpha)
    diffLyA1 = lines_rebin[0, idxLyA+1] - closest_lyA
    diffLyA2 = closest_lyA - lines_rebin[0, idxLyA-1]
    new_delta_lambda = (diffLyA1 + diffLyA2) / 2.0
    new_resolution = Lyalpha / new_delta_lambda
    smoothing_factor = float(Resolution) / new_resolution
    print 'original_delta_lambda and new_delta_lambda : %f, %f' % (original_delta_lambda, new_delta_lambda)
    print 'Original resolution = %i  -----  new resolution = %i' % (Resolution, int(new_resolution))
    print 'SMOOTHING FACTOR = %f' % (smoothing_factor)
        
    # Normalization
    norm_flx = lines_rebin[1,:] / cont_rebin[1,:]
    
    # Arrays for eqw measurements
    norm_lines = numpy.array([lines_rebin[0,:], norm_flx])
    cont_lines = spectrum.theo_cont(lines_rebin[0,:], scale_factor=1.0)
    
    # Measuring EQW
    lolim = 1210.67
    uplim = 1220.67
    eqw_fixed = spectrum.EQW_line_fixed(norm_lines, cont_lines, Lyalpha, 10.0)
    eqw = spectrum.EQW(norm_lines, cont_lines, lolim, uplim)
    print ('This is the eqw with the "fixed" function = %f' % (eqw_fixed))
    print ('This is the eqw with the EQW function = %f' % (eqw))
    
    # TEST for rebinning according to a desired delta_lambda
    '''delta_lambda=0.456 is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)'''
    desired_delta_lambda = 0.465
    testlines_rebin, test_cont_rebin = spectrum.rebin_arrays_for_desired_resolution(desired_delta_lambda, Lyalpha, lines, continuum)
    print('shape of testlines_rebin: %s    ----   shape of test_cont_rebin: %s' % (repr(testlines_rebin.shape), repr(test_cont_rebin.shape)))
    #testing_eqw_iter, lo_testing_eqw_iter, up_testing_eqw_iter = spectrum.EQW_iter(testlines_rebin, test_cont_rebin, Lyalpha)
    #print('testing_eqw_iter = %f, lo_testing_eqw_iter = %f, up_testing_eqw_iter = %f' % (testing_eqw_iter, lo_testing_eqw_iter, up_testing_eqw_iter))
    testeqw_fixed = spectrum.EQW_line_fixed(testlines_rebin, test_cont_rebin, Lyalpha, 10.0)
    print ('testeqw_fixed = %f' % (testeqw_fixed))

    '''
    # TEST with rebinning through interpolation
    # since the lines array is much larger than the continuum array, this last one must be rebinned:
    reference_wavelength = Lyalpha
    # Making the continuum array the same size as the lines array.
    rebinned_cont, rebinning_factor_cont = spectrum.rebinning_interpol(lines, continuum, reference_wavelength)
    print 'This is the test for the rebinning function: factor = %f' % (rebinning_factor_cont)
    norm_test = flx / rebinned_cont[1,:]
    test_norm_lines = numpy.array([wav, norm_test])
    test_norm_cont = spectrum.theo_cont(wav, scale_factor=1.0)
    test1_eqw = spectrum.EQW_line_fixed(test_norm_lines, test_norm_cont, Lyalpha, 10.0)
    test2_eqw = spectrum.EQW(test_norm_lines, test_norm_cont, lolim, uplim)
    print ('Test eqw with the EQW "fixed" function = %f' % (test1_eqw))
    print ('Test eqw with the EQW function = %f' % (test2_eqw))

    # ANOTHER TEST with rebinning
    rebin_lines_factor = 5.0
    rebin_cont_factor = spectrum.finding_rebin_it_factor_for_small_arr(lines, continuum, rebin_lines_factor)
    rebin_lines_array = spectrum.rebin_it(lines, rebin_lines_factor)    
    rebin_cont_array = spectrum.rebin_it(continuum, rebin_cont_factor)
    print 'Now testing the rebin_it function: factor for lines = %f, facotr for continuum %f' % (rebin_lines_factor, rebin_cont_factor)
    norm_anothertest = flx / rebinned_cont[1,:]
    anothertest_norm_lines = numpy.array([wav, norm_anothertest])
    anothertest_norm_cont = spectrum.theo_cont(wav, scale_factor=1.0)
    anothertest1_eqw = spectrum.EQW_line_fixed(anothertest_norm_lines, anothertest_norm_cont, Lyalpha, 10.0)
    anothertest2_eqw = spectrum.EQW(anothertest_norm_lines, anothertest_norm_cont, 1210.67, 1220.67)
    print ('Test eqw with the EQW "fixed" function = %f' % (anothertest1_eqw))
    print ('Test eqw with the EQW function = %f' % (anothertest2_eqw))
    '''
    # Ploting 
    # the line to mark Ly_alpha
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
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.ylim(0, 3e-8)
    pyplot.plot(wav, flx, 'b', lyalpha_arr_wav, flx, 'r--', wav_cont, flx_cont, 'g')
    
    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('Normalization zoom function')
    pyplot.xlim(low, up)
    pyplot.ylim(0.4, 1.5)
    pyplot.plot(lines_rebin[0], norm_flx, 'b')
    #pyplot.fill_between(eqw, eqw_y, 0,color='g')
    
    '''    
    pyplot.figure(3, figsize=(10, 10))
    pyplot.title('Normalization interpolation')
    pyplot.xlim(low, up)
    pyplot.ylim(0.4, 1.5)
    pyplot.plot(wav, norm_test, 'b')
    
    pyplot.figure(4, figsize=(10, 10))
    pyplot.title('Normalization rebin_it')
    pyplot.xlim(low, up)
    pyplot.ylim(0.4, 1.5)
    pyplot.plot(wav, norm_anothertest, 'b')
    '''
    pyplot.show()





