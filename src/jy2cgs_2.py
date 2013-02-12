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
def CMFGENfile_into2cols(path_CMFGEN, path_results, lower_wave, upper_wave):
    # Find directory
    ogrid_dir = glob.glob(path_CMFGEN)
    for single_dir in ogrid_dir:
        dir_name = os.path.join(single_dir+'/obs/')
        print'Going into directory: %s' % (dir_name)
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
            sp = spectrum.spectrum(file_in_use+"_1d")
            #sp.save2(file_in_use+"_Acgs.txt")
            output = (sp.A, sp.cgs)
            spectrum.write_1d(os.path.join(os.path.abspath(path_results), file_name+"_Acgs.txt"), output)
        
#### NOW DO THE WORK
# This is where the directory is
path_CMFGEN = '../ogrid_24jun09/test/*'
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
ogrid_dir = glob.glob(path_CMFGEN)

# Generate 1d files of flux and wavelength
#CMFGENfile_into2cols(path_CMFGEN, path_results, lower_wave, upper_wave)
#txt_files = [path_results+'*cont_1d', path_results+'*fin10_1d']
cont_files = glob.glob(os.path.join(os.path.abspath(path_results),'*cont*Acgs*'))
fin_files = glob.glob(os.path.join(os.path.abspath(path_results),'*_fin_10*Acgs*'))
#print cont_files
#print fin_files
#results = [fin_files, cont_files]
#for onedfiles in results:
# Load the text files
for f in range(0, len(fin_files)):
    A, cgs = numpy.loadtxt(fin_files[f], dtype=numpy.float64, unpack=True)    
    print('File loaded, now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    #print(A.shape, cgs.shape)
    #wav = A
    #flx = cgs
    wav, flx = selection(A, cgs, lower_wave, upper_wave)
    lines = numpy.array([wav, flx])

    A_cont, cgs_cont = numpy.loadtxt(cont_files[f], dtype=numpy.float64, unpack=True)    
    print('File loaded, now selecting desired range: %f to %f' % (lower_wave, upper_wave))
    # Selecting the range for Ly_alpha
    #print(A.shape, cgs.shape)
    #wav_cont = A_cont
    #flx_cont = cgs_cont
    wav_cont, flx_cont = selection(A_cont, cgs_cont, lower_wave, upper_wave)
    continuum = numpy.array([wav_cont, flx_cont])
        
    # Rebinning
    lines_rebin, cont_rebin, new_cont_factor, new_lines_factor = spectrum.get_factors_and_rebin(lines, continuum, 100)
    
    # Normalization
    norm_flx = lines_rebin[1] / cont_rebin[1]
    
    # Arrays for eqw measurements
    norm_lines = numpy.array([lines_rebin[0,:], norm_flx])
    cont_lines = spectrum.theo_cont(lines_rebin[0,:], scale_factor=1.0)
    
    # Measuring EQW
    eqw_fixed = spectrum.EQW_lyA_fixed(norm_lines, cont_lines, 10.0)
    eqw = spectrum.EQW(norm_lines, cont_lines, 1210.6, 1220.6)
    print ('This is the eqw with the "fixed" function = %f' % (eqw_fixed))
    print ('This is the eqw with the EQW function = %f' % (eqw))
    
    # Ploting 
    # the line to mark Ly_alpha
    lyalpha_arr_wav = []
    for i in flx:
        lyalpha_arr_wav.append(Lyalpha)  
    
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
    #pyplot.suptitle(aka)
    pyplot.plot(wav, flx, 'b', lyalpha_arr_wav, flx, 'r--', wav_cont, flx_cont, 'g')
    
    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('Normalization')
    pyplot.plot(lines_rebin[0], norm_flx, 'b')
    
    pyplot.show()





