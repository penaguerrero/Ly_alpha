import numpy
import glob
from spectrum import pyplot
from spectrum import spectrum


'''
This program uses the already converted files of pollux (from CMFGEN format into two-column Angstrom-cgs file) 
to create a selection text file that needs to be converted in order for splot to read it.

'''

def selection(A, cgs, lower, upper):
    # This functions selects according to a wavelength range
    # A and cgs are numpy arrays
    xlist = []
    ylist = []
    for i in range(0, len(A)):
        if A[i] >= lower and A[i] <= upper:
            xlist.append(A[i])
            ylist.append(cgs[i])     
    A_selected = numpy.array(xlist)
    cgs_selected = numpy.array(ylist)
    return(A_selected, cgs_selected)


# The SEDs and their txt files
# the full path is: /Users/pena/Documents/AptanaStudio3/Ly_alpha... 

### The SED files are the ones that contain the data
### their corresponding text files contain detailed model information: metallicity, luminosity, etc.

seds = ['../pollux/C_s27500g3.00z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-6.02Vinfty1530beta0.8finfty1vcl0.sed',
        '../pollux/C_s40000g3.50z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-4.85Vinfty2644beta0.8finfty1vcl0.sed']
#seds = glob.glob('../pollux/*.sed')

#    txts = ['../pollux/C_s27500g3.00z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-6.02Vinfty1530beta0.8finfty1vcl0.sed.txt',
#            '../pollux/C_s40000g3.50z0.0t5.0_a0.00c0.00n0.00o0.00_Mdot-4.85Vinfty2644beta0.8finfty1vcl0.sed.txt']
#txts = glob.glob('../pollux/*.sed.txt')

path_results = '../results/'

star_counter = 0
lower_cut = 1190.0
upper_cut = 1300.0

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.6701

for i in range(0, len(seds)):
    print('Loading %s' % (seds[i]))
    aka='A.K.A. star'+repr(star_counter)
    print(aka)    
    star_counter = star_counter + 1
    
    # Load the file
    A, cgs = numpy.loadtxt(seds[i], dtype=float, unpack=True)
    print('File loaded, now selecting desired range: %f to %f' % (lower_cut, upper_cut))
    
    # Selecting the range for Ly_alpha
    wav, flx = spectrum.selection(A, cgs, lower_cut, upper_cut)
    
    # Ploting 
    # the line to mark Ly_alpha
    lyalpha_arr_wav = []
    for i in flx:
        lyalpha_arr_wav.append(Lyalpha)  
    
    # plot limits
    low = 1200
    up = 1230
    # Figure
    pyplot.xlabel('Wavelength [$\AA$]')
    pyplot.ylabel('Flux [ergs/s/cm$^2$/$\AA$]')
    pyplot.xlim(low, up)
    pyplot.suptitle(aka)
    pyplot.plot(wav, flx, 'b')
    pyplot.plot(lyalpha_arr_wav, flx, 'r--')

    pyplot.show()
      
    
    
    
