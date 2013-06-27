import numpy
#import glob
import os
#import string
from matplotlib import pyplot
from pylab import *
from spectrum import spectrum
from matplotlib.ticker import MaxNLocator


'''
This program simply produces nice plots of the selected stars. 

### It reads the text files of the selected wavelength range. These files were obtained with OpolluxStars2.py code, along with
### the nice table that prints the EWs is in OpolluxStars2.py

'''

def yarr_for_line(arr, line, above=False, below=False):
    ''' This function defines a little array to make the line that marks the ID of each line.
    line = the rest wavelength of the line of interest
    arr = two dimensional array of wavelength and flux where to find the line
    It RETURNS: the y array of each line
    '''
    if (above==False) and (below==False):
        print('both false')
        line_in_arr, _ = spectrum.find_nearest(arr[0], line)
        flx_line_in_arr = spectrum.findXinY(arr[1], arr[0], line_in_arr)
        if flx_line_in_arr > 1.03:
            y_arr_of_line_in_arr = [flx_line_in_arr+0.1, flx_line_in_arr+0.3]
        else:
            y_arr_of_line_in_arr = [0.6, 0.4]
        return(y_arr_of_line_in_arr)
    elif above==True:
        y_arr_of_line_in_arr = [1.2, 1.4]
        return(y_arr_of_line_in_arr)        
    elif below==True:
        y_arr_of_line_in_arr = [0.6, 0.4]
        return(y_arr_of_line_in_arr)        

def arr_ref_lines(y_line, line):
    ''' This function creates the array to mark in the plot the reference lines
    y_line = is the y array created by the yarr_for_line function
    line = the rest wavelength of the line of interest
    RETURNS: the array for plotting
    '''
    y_line_list = []
    for _ in range(0, len(y_line)):
        y_line_list.append(line)
    arr_ref_line = numpy.array([y_line_list, y_line])
    return(arr_ref_line)

def ref_line(arr, line, above, below):
    y_line = yarr_for_line(arr, line, above, below)
    y_line_arr = arr_ref_lines(y_line, line)
    return(y_line_arr)


path_results = '../results/OpolluxStars/'
path_plots = '../results/plots/'

# Rest wavelength of Ly_alpha from NIST
Lyalpha = 1215.3 #1215.6701
# rest wavelengths from Leitherer et al(2011), ApJ, 141, 37
# the numbers are different because they look prettier in the plot
cIII_mie = 1175.53
siIII = 1206.50
nV = 1242.80
cIII = 1247.38
sII = 1250.58

stars = ['NormOstar_36_40000g3.50.txt']

for s in stars:
    star = path_results+s
    base_name1 = (os.path.basename(star))
    base_name2 = base_name1.replace(".txt", "")
    # Load the file
    norm_lines = numpy.loadtxt(star, dtype=float, unpack=True)
    print('File loaded, now working with: %s' % (star))

    # Marking the relevant lines
    lyalpha_arr = ref_line(norm_lines, Lyalpha, above=False, below=False)
    cIIImie_arr = ref_line(norm_lines, cIII_mie, above=True, below=False)
    siIII_arr = ref_line(norm_lines, siIII, above=True, below=False)
    nV_arr = ref_line(norm_lines, nV, above=False, below=True)
    cIII_arr = ref_line(norm_lines, cIII, above=False, below=False)
    sII_arr = ref_line(norm_lines, sII, above=False, below=False)
    
    # the line that marks the where the EW was inegrated from
    int_limits = [Lyalpha-5.0, Lyalpha+5.0]
    y_int = [1.0, 1.0]
    integration = numpy.array([int_limits, y_int])
    
    # line of the continuum
    continuum_plot = []
    for i in norm_lines[0]:
        continuum_plot.append(1.0)

    low = 1160
    up = 1260
    lo_y = -0.05

    # making all fonts biger
    font = {#'family' : 'Vera Sans',
            'weight' : 'regular',
            'size'   : 17}
    matplotlib.rc('font', **font)
    
    N = len(lyalpha_arr[0])
    pyplot.figure(3, figsize=(10, 6))
    xlab = pyplot.xlabel('Wavelength [$\AA$]')
    ylab = pyplot.ylabel('Normalized Flux')# [ergs/s/cm$^2$/$\AA$]')
    # add some space between labels and axis
    xlab.set_position((0.5, 0.02))
    ylab.set_position((0.9, 0.5))
    pyplot.xlim(low, up)
    pyplot.plot(norm_lines[0], norm_lines[1], 'k', lyalpha_arr[0], lyalpha_arr[1], 'r--')
    pyplot.plot(norm_lines[0], continuum_plot, 'r:')         #continuum
    pyplot.plot(integration[0], integration[1], 'r', lw=2)  # integration interval
    # remove the firt tick so that they do not overlap
    pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
    # mark other important lines
    pyplot.plot(cIIImie_arr[0], cIIImie_arr[1], 'r--', siIII_arr[0], siIII_arr[1], 'r--', nV_arr[0], nV_arr[1], 'r--')
    pyplot.plot(cIII_arr[0], cIII_arr[1], 'r--')#, sII_arr[0], sII_arr[1], 'r--')
    pyplot.text(cIII_mie-2, cIIImie_arr[1][N-1]+0.03, 'C III')
    pyplot.text(siIII-2, siIII_arr[1][N-1]+0.03, 'Si III')
    # When the leyend looks better above
    #pyplot.text(Lyalpha-4, y_Lyalpha[N-1]+0.03, 'Ly-alpha')
    #pyplot.text(nV-2.5, y_Lyalpha[N-1]+0.03, 'N V')
    #pyplot.text(cIII-1, y_Lyalpha[N-1]+0.03, 'C III')
    #pyplot.text(sII-1, y_Lyalpha[N-1]+0.03, 'S II')
    # When the legend looks better below
    pyplot.text(Lyalpha-4, 0.4-0.09, 'Ly-alpha')
    pyplot.text(nV-3.5, 0.4-0.09, 'N V')
    pyplot.text(cIII-0.5, 0.4-0.09, 'C III')
    epsfile = os.path.join(path_plots, base_name2+".eps")
    pyplot.savefig(epsfile)
    pyplot.show()    
    
    print("Do you want to save this plot?  [y/N , meaning NO is default... :P ]")
    save_plt = raw_input()
    pyplot.ioff()
    if save_plt == 'n' or save_plt == '':
        print('Plot not saved...   :( ')
        os.remove(epsfile)
    else:
        print('Star %s was saved.' % star)
    
    # Adding theoretical continuum
    continuum = spectrum.theo_cont(norm_lines[0], scale_factor=1.0)
    
    # EQW measurement
    eqw = spectrum.EQW_line_fixed(norm_lines, continuum, Lyalpha, 10.0)
    print('EQW = %f' % eqw) 
    

print "Done!"
