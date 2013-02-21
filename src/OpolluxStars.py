import numpy
import glob
from spectrum import pyplot
from spectrum import spectrum


'''
This program uses the already converted files of pollux (from CMFGEN format into two-column Angstrom-cgs file) 
to create a selection text file that needs to be converted in order for splot to read it.

'''



# The SEDs and their txt files
# the full path is: /Users/pena/Documents/AptanaStudio3/Ly_alpha... 

### The SED files are the ones that contain the data
### their corresponding text files contain detailed model information: metallicity, luminosity, etc.

#seds = ['../pollux/C_s27500g3.00z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-6.02Vinfty1530beta0.8finfty1vcl0.sed',
#        '../pollux/C_s40000g3.50z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-4.85Vinfty2644beta0.8finfty1vcl0.sed']
seds = glob.glob('../pollux/*.sed')

#    txts = ['../pollux/C_s27500g3.00z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-6.02Vinfty1530beta0.8finfty1vcl0.sed.txt',
#            '../pollux/C_s40000g3.50z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-4.85Vinfty2644beta0.8finfty1vcl0.sed.txt']
#txts = glob.glob('../pollux/*.sed.txt')

path_results = '../results/OpolluxStars/'

star_counter = 0
lower_cut = 1190.0
upper_cut = 1300.0

#a = open(path_results+'aka_dictionary.txt', 'w+')
#print >> a, 'AKA', 'Real name'

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701

for i in range(0, len(seds)):
    print('Loading %s' % (seds[i]))
    aka='A.K.A. star'+repr(star_counter)
    print(aka)    
    
    ## Printing master aka and real name file
    #print >> a, 'Ostar_'+repr(star_counter), seds[i]

    # Load the file
    A, cgs = numpy.loadtxt(seds[i], dtype=float, unpack=True)
    #print('File loaded, now selecting desired range: %f to %f' % (lower_cut, upper_cut))
    
    # Selecting the range for Ly_alpha
    wav, flx = spectrum.selection(A, cgs, lower_cut, upper_cut)
    data = numpy.array([wav, flx])
    
    # Rebinning data to delta_lambda=0.456, which is an average between the typical real observed dlta_lambdas of STIS (0.75) and COS (0.18)
    desired_delta_lambda = 0.465
    data_rebinned, smoothing_R_factor = spectrum.rebin_one_arr_to_desired_resolution(desired_delta_lambda, Lyalpha, data, guessed_rows=500)
    print('shape of data_rebinned: %s  ---   smoothing_R_factor = %f' % (repr(data_rebinned.shape), smoothing_R_factor))

    # Printing text file of selected spectrum
    #f = open(path_results+'selec_data_Ostar_'+repr(star_counter)+'.txt', 'w+')
    #for i in range(0, len(data_rebinned[0])):
    #    print >> f, data_rebinned[0, i], data_rebinned[1, i]
    #f.close()
    
    star_counter = star_counter + 1

    '''
    # Ploting 
    # the line to mark Ly_alpha
    lyalpha_arr_wav = []
    for i in flx:
        lyalpha_arr_wav.append(Lyalpha)  
    
    # plot limits
    low = 1200
    up = 1230
    # Figure
    pyplot.figure(1)
    pyplot.title('Non-rebinned data')
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.suptitle(aka)
    pyplot.plot(wav, flx, 'b')
    pyplot.plot(lyalpha_arr_wav, flx, 'r--')

    pyplot.figure(2)
    pyplot.title('Rebinned data')
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.suptitle(aka)
    pyplot.plot(data_rebinned[0], data_rebinned[1], 'b')
    pyplot.plot(lyalpha_arr_wav, flx, 'r--')
    
    pyplot.show()
    '''  
#a.close()
   
    
    
