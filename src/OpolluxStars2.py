import numpy
import glob
import os
import string
import table
from spectrum import pyplot
from spectrum import spectrum


'''
This program uses the already converted files of pollux (from CMFGEN format into two-column Angstrom-cgs file) 
to create a selection text file that needs to be converted in order for splot to read it.

'''


def find_corr_points(x1, x2, arr):
    ''' This function fits a straight line to the given points, using the array coordinates
    ### arr is an array of wavelength and flux '''
    temp_x1, _ = spectrum.find_nearest(arr[0], x1)
    temp_x2, _ = spectrum.find_nearest(arr[0], x2)
    print('closest wavelengths: %f, %f' % (temp_x1, temp_x2))
    y1 = spectrum.findXinY(arr[1], arr[0], temp_x1)
    y2 = spectrum.findXinY(arr[1], arr[0], temp_x2)
    return (y1, y2)

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

EQWs = []
points_used_for_line_continuums = []
#points_used_for_line_continuums.append('Wavelengths    Fluxes')
#points_used_for_line_continuums.append('x1  x2    y1  y2')
AKAs = []

for i in range(0, len(seds)):
    print('Loading %s' % (seds[i]))
    base_seds = os.path.basename(seds[i])
    base_name2 = string.split(base_seds, sep='z')
    base_name3 = string.split(base_name2[0], sep='g')
    Teff = base_name3[0].replace("C_s", "")
    logg = base_name3[1]
    print(Teff, logg)
    aka='Ostar'+repr(star_counter)+'_'+Teff+'g'+logg
    print(aka) 
    AKAs.append(aka)   
    
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
    
    # Ploting 
    # the line to mark Ly_alpha
    lyalpha_arr_wav = []
    for i in data_rebinned[1]:
        lyalpha_arr_wav.append(Lyalpha)  
    
    # plot limits
    low = 1190
    up = 1300
    # Figure
    '''
    pyplot.figure(1, figsize=(10, 10))
    pyplot.title('Rebinned data')
    pyplot.suptitle(aka)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.plot(data_rebinned[0], data_rebinned[1], 'b', lyalpha_arr_wav, data_rebinned[1], 'r--')
    pyplot.show()

    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('Non-rebinned data')
    pyplot.suptitle(aka)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.plot(wav, flx, 'b')
    pyplot.plot(lyalpha_arr_wav, flx, 'r--')
    '''

    # My continuum function, two-point linear equation: y = m*(x -x1) + y1 
    print('Enter the lower and upper X-axis to consider in the line:')
    x1 = 1274.63 #float(raw_input("Low limit: "))
    x2 = 1290.42 #float(raw_input("Upper limit: "))
    y1, y2 = find_corr_points(x1, x2, data_rebinned)
    temp_obj = string.split(base_seds, sep='_')
    temp_obj2 = string.split(temp_obj[1], sep='.')
    temp = temp_obj2[0].replace("s", "")
    temp2 = string.split(temp, sep='g')
    print (temp2[0])  
    if (y1 / y2) > 2.0 or (numpy.fabs(y2) > numpy.fabs(y1)):        
        x1 = 1271.61 #1283.15
        x2 = 1298.58
        y1, y2 = find_corr_points(x1, x2, data_rebinned)         
    if (float(temp2[0]) > 35000.0) and (float(temp2[0]) < 48000.0):
        x1 = 1221.76
        x2 = 1298.58
        y1, y2 = find_corr_points(x1, x2, data_rebinned)
        if (numpy.fabs(y2) > numpy.fabs(y1)):
            x1 = 1204.37
            x2 = 1298.58
            y1, y2 = find_corr_points(x1, x2, data_rebinned)         
            
    if (float(temp2[0]) > 48000.0):
        x1 = 1202.95
        x2 = 1298.58
        y1, y2 = find_corr_points(x1, x2, data_rebinned)

    m = (y2 - y1) / (x2 - x1)
    y_list = []
    for x in data_rebinned[0]:
        y = m * (x - x1) + y1
        y_list.append(y)
    my_cont_arr = numpy.array([data_rebinned[0], y_list])
    
    coords = [x1, x2, y1, y2]
    points_used_for_line_continuums.append(coords)    
    
    # Normalization to my continuum
    norm_flx = data_rebinned[1] / my_cont_arr[1]
    norm_lines = numpy.array([data_rebinned[0], norm_flx])
    
    '''
    lyalpha_arr_norm = []
    for i in norm_lines[1]:
        lyalpha_arr_norm.append(Lyalpha)  
     
    # Fitted plots 
    pyplot.figure(2, figsize=(10, 10))
    pyplot.title('My continuum fit')
    pyplot.suptitle(aka)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    pyplot.plot(data_rebinned[0], data_rebinned[1], 'b', lyalpha_arr_wav, data_rebinned[1], 'r--', my_cont_arr[0], my_cont_arr[1], 'g')

    pyplot.figure(3, figsize=(10, 10))
    pyplot.title('My continuum fit')
    pyplot.suptitle(aka)
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(1190, 1300)
    pyplot.plot(norm_lines[0], norm_lines[1], 'b', lyalpha_arr_norm, norm_lines[1], 'r--')

    pyplot.show()
    '''
    # Adding theoretical continuum
    continuum = spectrum.theo_cont(norm_lines[0], scale_factor=1.0)
    
    # EQW measurement
    eqw = spectrum.EQW_line_fixed(norm_lines, continuum, Lyalpha, 10.0)
    print('EQW = %f' % eqw) 
    EQWs.append(eqw)
    
    star_counter = star_counter + 1

#a.close()

Oeqws_file = open(path_results+'Ostar_eqws.txt', 'w+')
table_Ostars = []
print >> Oeqws_file, 'SIGN CONVENTION: positive = emission,  negative = absorption'

for i in range(0, len(EQWs)):
    print(AKAs[i]+'EQW = %f' % EQWs[i])
    temp = [AKAs[i], repr(EQWs[i])]
    table_Ostars.append(temp)
    
table.pprint_table(Oeqws_file, table_Ostars)
Oeqws_file.close()


    
