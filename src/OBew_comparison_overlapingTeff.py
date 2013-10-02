import os
import numpy
import string
import science
from scipy import stats
from matplotlib import pyplot
from matplotlib import pylab
from matplotlib.ticker import MaxNLocator
from pylab import *

'''
This code compares the EWs of O and B stars from spherical CMFGEN, plane-parallel CMFGEN, and TLUSTY models. It plots 
EW as a function of temperature for an effective temperature.
'''

############################################################################################################

''' These are all the switches than can be turned on or off in this code: '''
# This is the machine where I am working
#machine = 'gallo'
machine = 'pena'

# Create the file with the EW differences of CMFGEN, TLUTY, and plane-parallel CMFGEN
create_diffs_file = True

# Make and show the comparison plot of EW versus Teff as a function of log g.
make_comparison_plot = True
show_plot = True    # If show_plot is active, the saved plot will be empty.
save_plot = True   # To save the plot show_plot must be off
plot_name = 'OvsB_EW_Teff29k_and_higher.eps'

############################################################################################################

# Get the paths to obtain the data
path_plots = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/plots/'
path_results = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/'
#table_OBew = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/diff_eqws_with_mtvel.txt'
table_OBew = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/diff_eqws_with_mtvel copy_mod.txt'
table_planepar = '/Users/'+machine+'/Documents/AptanaStudio3/Ly_alpha/results/CMFGENplane_par/planepar_EWs.txt'
data_stars = science.utility.DataDir(table_OBew)
data_plots = science.utility.DataDir(path_plots)
data_results = science.utility.DataDir(path_results)
data_planepar = science.utility.DataDir(table_planepar)
if not os.path.exists(data_results.path):
    print("Path: %s does not exist!" % (data_results.path))
    exit(1)

# obtain CMFGEN and Tlusty data
Teff, logg, mtvel, lya5, lyahalf = numpy.loadtxt(table_OBew, skiprows=1, usecols=(0,1,2,3,5), unpack=True, dtype=numpy.float64)
# get the names of the stars
names_column = open(table_OBew, 'r')
star_names = []
for line in names_column:
    line_list = string.split(line)
    star_name = line_list[8]
    star_names.append(star_name)
names_column.close()
star_names.pop(0) # not to include the colum header row

# This is all the data for the O and B stars from CMFGEN and TLUSTY
logg_round = numpy.around(logg, decimals=0)
data = numpy.array([Teff, logg_round, mtvel, lya5, lyahalf], dtype=numpy.float64)

# separate CMFGEN from TLUSTY
Btlusty = numpy.array([Teff[lyahalf!=0.0], logg[lyahalf!=0.0], lyahalf[lyahalf!=0.0]])
Ocmfgen = numpy.array([Teff[lya5!=0.0], logg[lya5!=0.0], lya5[lya5!=0.0]])

# obtain the CMFGEN plane-parallel data
ppTeff, pplogg, ppsimpleew, pphalfew = numpy.loadtxt(table_planepar, skiprows=1, unpack=True)
# all the data for plane-parallel CMFGEN
ppCMFGEN = numpy.array([ppTeff, pplogg, ppsimpleew, pphalfew])

EW_diff_tlusty_list = [] 
EW_diff_tlusty_list_simple = []
EW_diff_Ocmfgen_list = []
line_list = []

# Find overlaping temperatures between Tlusty and Plane-Parallel CMFGEN
for i in range(len(ppCMFGEN[0])):
    for j in range(len(Btlusty[0])):
        if ppCMFGEN[0][i] == Btlusty[0][j]:
            if ppCMFGEN[1][i] == Btlusty[1][j]:
                Oempty = 'No data'
                ew_overlap = numpy.fabs(Btlusty[2][j]) - numpy.fabs(ppCMFGEN[3][i])
                EW_diff_tlusty_list.append(ew_overlap)
                ew_overlap_simple = numpy.fabs(Btlusty[2][j]) - numpy.fabs(ppCMFGEN[2][i])
                EW_diff_tlusty_list_simple.append(ew_overlap_simple)
                line = ('%i    %0.2f    %s    %0.2f    %0.2f    %0.2f    %0.2f\n' 
                     % (ppCMFGEN[0][i], ppCMFGEN[1][i], Oempty, Btlusty[2][j], ppCMFGEN[3][i], ew_overlap, ew_overlap_simple))
                line_list.append(line)

# Find overlaping temperarures for Tlusty and CMFGEN, where Tlusty is hotter than 26,000
for i in range(len(Btlusty[0])):
    for j in range(len(Ocmfgen[0])):
        if Btlusty[0][i] == Ocmfgen[0][j]:
            if Btlusty[1][i] == Ocmfgen[1][j]:
                ppempty = 'No data'
                ew_overlap = numpy.fabs(Ocmfgen[2][j]) - numpy.fabs(Btlusty[2][i])
                EW_diff_Ocmfgen_list.append(ew_overlap)
                line = ('%i    %0.2f    %0.2f      %0.2f    %s    %0.2f    %s\n' 
                     % (Btlusty[0][i], Btlusty[1][i], Ocmfgen[2][j], Btlusty[2][i], ppempty, ew_overlap, ppempty))
                line_list.append(line)


# Find the maximum and minimum difference
max_diff_plane_with_B = max(EW_diff_tlusty_list)
min_diff_plane_with_B = min(EW_diff_tlusty_list)
mean_diff_plane_with_B = sum(EW_diff_tlusty_list) / len(EW_diff_tlusty_list)
EW_diff_tlusty_arr = numpy.array([EW_diff_tlusty_list])
EEW_diff_tlusty_round = numpy.around(EW_diff_tlusty_arr, decimals=0)
mode_diff_plane_with_B = stats.mode(EEW_diff_tlusty_round, axis=None)

max_EW_diff_tlusty_list_simple = max(EW_diff_tlusty_list_simple)
min_EW_diff_tlusty_list_simple = min(EW_diff_tlusty_list_simple)
mean_EW_diff_tlusty_list_simple = sum(EW_diff_tlusty_list_simple) / len(EW_diff_tlusty_list_simple)
EW_diff_tlusty_list_simple_arr = numpy.array([EW_diff_tlusty_list_simple])
EW_diff_tlusty_list_simple_round = numpy.around(EW_diff_tlusty_list_simple_arr, decimals=0)
mode_EW_diff_tlusty_list_simple = stats.mode(EW_diff_tlusty_list_simple_round, axis=None)

max_diff_B_with_O = max(EW_diff_Ocmfgen_list)
min_diff_B_with_O = min(EW_diff_Ocmfgen_list)
mean_diff_B_with_O = sum(EW_diff_Ocmfgen_list) / len(EW_diff_Ocmfgen_list)
EW_diff_Ocmfgen_arr = numpy.array([EW_diff_Ocmfgen_list])
EW_diff_Ocmfgen_round = numpy.around(EW_diff_Ocmfgen_arr, decimals=0)
mode_diff_B_with_O = stats.mode(EW_diff_Ocmfgen_round, axis=None)

#print stuff
print('max_diff_plane_with_B = %0.2f' % max_diff_plane_with_B)
print('min_diff_plane_with_B = %0.2f' % min_diff_plane_with_B)
print EEW_diff_tlusty_round
print(len(EW_diff_tlusty_list))
print('mode_diff_plane_with_B = ', mode_diff_plane_with_B)
print('mean_diff_plane_with_B = %0.2f' % mean_diff_plane_with_B)

print('max_EW_diff_tlusty_list_simple = %0.2f' % max_EW_diff_tlusty_list_simple)
print('min_EW_diff_tlusty_list_simple = %0.2f' % min_EW_diff_tlusty_list_simple)
print EW_diff_tlusty_list_simple_round
print(len(EW_diff_tlusty_list_simple))
print('mode_EW_diff_tlusty_list_simple = ', mode_EW_diff_tlusty_list_simple)
print('mean_EW_diff_tlusty_list_simple = %0.2f' % mean_EW_diff_tlusty_list_simple)

print('max_diff_B_with_O = %0.2f' % max_diff_B_with_O)
print('min_diff_B_with_O = %0.2f' % min_diff_B_with_O)
print EW_diff_Ocmfgen_round
print len(EW_diff_Ocmfgen_list)
print('mode_diff_B_with_O = ' , mode_diff_B_with_O)
print('mean_diff_B_with_O = %0.2f' % mean_diff_B_with_O)

if create_diffs_file == True:
    f = open(os.path.join(data_results.path,'diff_in_EW_TlustyCmfgenPlanePar.txt'), 'w+')
    print >> f, 'Teff    log g    O_ew        B_ew    PlanePar_ew    Abs_Diff*    Abs_Diff_simple**\n'
    print >> f, '   * Abs_Diff = |O_ew| - |B_ew|'
    print >> f, '    OR        = |B_ew| - |PlanePar_ew(with half-integration method)|\n'
    print >> f, '   ** Abs_Diff_simple = |B_ew| - |PlanePar_ew(with simple-integration method)|\n'
    print >> f,('Maximum difference of between B_ew and PlanePar_ew = %0.2f' % max_diff_plane_with_B)
    print >> f,('Minimum difference of between B_ew and PlanePar_ew = %0.2f' % min_diff_plane_with_B)
    print >> f,('Maximum difference of between O_ew and B_ew = %0.2f' % max_diff_B_with_O)
    print >> f,('Minimum difference of between O_ew and B_ew = %0.2f\n' % min_diff_B_with_O)
    for line in line_list:
        f.write(line)
    f.close()


if make_comparison_plot == True:
    # making all fonts biger
    font = {#'family' : 'regular',
            'weight' : 'regular',
            'size'   : 17}
    matplotlib.rc('font', **font)

    fig = pyplot.figure(1, figsize=(10,10))
    pyplot.title('Comparison of EWs')
    pyplot.xlabel('Effective Temperature [K]')
    pyplot.ylabel('Equivalent Width [$\AA$]')
    b, = pyplot.plot(Btlusty[0], Btlusty[2], 'bo')
    o, = pyplot.plot(Ocmfgen[0], Ocmfgen[2], 'k^')
    pp, = pyplot.plot(ppCMFGEN[0], ppCMFGEN[3], 'r*', markersize=10)
    # remove the firt tick so that they do not overlap
    pyplot.gca().yaxis.set_major_locator(MaxNLocator(9, prune='lower'))
    # make the leyend for the data sets
    codes = ["TLUSTY", "CMFGEN", "Plane-parallel CMFGEN"]
    data_sets = []
    data_sets.append(b)
    data_sets.append(o)
    data_sets.append(pp)    
    leg = pyplot.figure(1).legend(data_sets, codes, title='Codes', bbox_to_anchor=(0.87, 0.35), labelspacing=0.2)
    for t in leg.get_texts():
        t.set_fontsize(17)     
    '''
    # This is an equivalent way to make the legend
    pylab.legend((r'TLUSTY', r'CMFGEN', r'Plane-Parallel CMFGEN'), title='Codes', bbox_to_anchor=(0.87, 0.35))
    ltext = pylab.gca().get_legend().get_texts()
    pylab.setp(ltext[0], fontsize=17, color='b')
    pylab.setp(ltext[1], fontsize=17, color='k')
    pylab.setp(ltext[2], fontsize=17, color='r')
    '''
    if show_plot == False:
        pyplot.show()
    if save_plot ==  True:
        #save the plot
        epsfile = os.path.join(data_plots.path, plot_name)
        pyplot.savefig(epsfile)
    

print 'Code has finished!'


