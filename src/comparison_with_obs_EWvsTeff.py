import numpy
import os
from matplotlib import pyplot
from pylab import *
from matplotlib.ticker import MaxNLocator


'''
This code makes a nice plot of EW vs. effective temperature. It compares our EW determinations to those of Savage & Code (1970)
and Savage & Panek (1974).
'''

table_OBew = '../results/table_OB_desired_eqws.txt'
path_plots = '../results/plots/'
plot_name = 'comparison_with_obs_EWvsTeff3.eps'

# There are 3 types of data for this plot: CMFGEN, TLUSTY, and OBSERVATIONS
data_type = ["TW(CMFGEN)", "TW(TLUSTY)", "SC70 E(B-V)<0.05", "SC70 0.05<E(B-V)<0.15", "SC70 0.15<E(B-V)<0.25", "SC70 E(B-V)>0.25", "SP74"]

# This is the data from Svage & Code (70) = sc70
teff_sc70 = [15170., 24400, 30000., 18700., 17620., 26850., 18200., 25290, 21580., 13160., 31000., 33400., 30000., 29300., 32000., 
             24820., 21480., 31270., 29910., 26390., 25180., 20990., 12250., 23930., 17030., 25000., 23000., 20500., 15780., 23100., 
             22900., 17900., 16860., 26200., 16980., 19690., 21380., 21380., 25350., 31000., 28000., 24200., 33000., 22500., 26300., 
             24720., 20180., 14140.]
ew_sc70 = [28., 25., 14., 29., 25., 14., 20., 19., 15., 35., 20., 13., 31., 33., 17., 24., 21., 16., 25., 22., 10., 12., 43., 11., 22., 12.,
           20., 15., 36., 14., 17., 26., 22., 14., 26., 22., 30., 12., 18., 24., 33., 31., 26., 20., 12., 16., 29., 33.]
# E(B-V)
ebv = [0.04, 0.31, 0.10, 0.07, 0.02, 0.01, 0.09, 0.05, 0.01, 0.01, 0.09, 0.03, 0.13, 0.36, #0.21, 
       0.07, 0.05, 0.06, 0.06, 0.06, 0.04, 0.02, 0.00, 0.00, 0.02, 0.01, 0.03, 0.04, 0.04, 0.01, 
       0.05, 0.03, 0.01, 0.02, 0.01, 0.00, 0.05, 0.17, 0.21, 0.07, 0.19, 0.22, 0.40, 0.32, 0.02, 0.06, 
       0.03, 0.05, 0.01]

# This is the data from Svage & Panek (74) = sp74
teff_sp74 = [14400., 14500., 16470., 12800., 15650., 17940., 12100., 13160., 14200., 18000., 18260., 17200., 12250., 17030., 19525., 18280., 
             11500., 15780., 25000., 17900., 16860., 19690., 23440., 14890., 13400., 17620., 18890., 17850., 14140.]
ew_sp74 = [45.0, 21.0, 25.0, 65., 31., 22., 30., 42., 45., 21., 20., 19., 55., 25., 10., 18., 23., 45., 20., 30., 
           21., 19., 11., 35., 44., 21., 18., 17., 45.]

# Now put temperatures in increasing order for both sets of observations and
# first observations set
dtype = [('teff', float), ('ew', float), ('ebv', float)]
values = []
for i in range(0, len(teff_sc70)):
    pair = (teff_sc70[i], ew_sc70[i], ebv[i])
    values.append(pair)
obs_sc70 = numpy.array(values, dtype=dtype)
ordered_sc70 = numpy.sort(obs_sc70, order='teff')
# undo the tuple and convert EW to negatives according to convention: NEGATIVE = ABSORPTION
teff_sc70o = []
ew_sc70on = []
ebv_o = []
for i in range(0, len(ordered_sc70)):
    t = ordered_sc70[i][0]
    e = ordered_sc70[i][1] * (-1)
    eb = ordered_sc70[i][2]
    teff_sc70o.append(t)
    ew_sc70on.append(e)
    ebv_o.append(eb)

# second observations set
dtype = [('teff', float), ('ew', float)]
values = []
for i in range(0, len(teff_sp74)):
    pair = (teff_sp74[i], ew_sp74[i])
    values.append(pair)
obs_sp74 = numpy.array(values, dtype=dtype)
ordered_sp74 = numpy.sort(obs_sp74, order='teff')
# undo the tuple and convert EW to negatives according to convention: NEGATIVE = ABSORPTION
teff_sp74o = []
ew_sp74on = []
for i in range(0, len(ordered_sp74)):
    t = ordered_sp74[i][0]
    e = ordered_sp74[i][1] * (-1)
    teff_sp74o.append(t)
    ew_sp74on.append(e)


# OUR EW DETERMINATIONS, these are alreade ordered: increasing temperature
# Reading the text table as numpy arrays
data = numpy.loadtxt(table_OBew, skiprows=3, usecols=(0,1,2,3,4,5), unpack=True)
#data: 0=temperature, 1=logg, 2=LyAew+-5, 3=LyAew+-20, 4=(LyAew+20)*2, 5=difference between cols 3 and 4

# separate the data according to the model used.
teff_cmfgen = []
teff_tlusty = []
ew_cmfgen = []
ew_tlusty = []
for i in range(0, len(data[2])):
    t = data[0][i]
    ec = data[2][i]
    et = data[3][i]
    if ec != 0:
        teff_cmfgen.append(t)
        ew_cmfgen.append(ec)
    elif et != 0:
        teff_tlusty.append(t)
        ew_tlusty.append(et)


# making all fonts biger
font = {#'family' : 'Vera Sans',
        'weight' : 'regular',
        'size'   : 17}
matplotlib.rc('font', **font)

# plotting
colors = ['g', 'm', 'c', 'y']
fig = pyplot.figure(1, figsize=(10, 10))
xlab = pyplot.xlabel('Effective Temperature [K]')
ylab = pyplot.ylabel('Lyman-alpha Equivalent Width [$\AA$]')
# add some space between labels and axis
xlab.set_position((0.5, 0.02))
ylab.set_position((0.9, 0.5))
# change the symbols according to  cmfgen (empty with mfc='None') or tlusty (filled)
c1 = pyplot.plot(teff_cmfgen, ew_cmfgen, 'w^')
c2 = pyplot.plot(teff_tlusty, ew_tlusty, 'wo')
#c3 = pyplot.plot(teff_sc70o, ew_sc70on, 'r>')
c4 = pyplot.plot(teff_sp74o, ew_sp74on, 'b*', ms=7)
# IF E(B-V) > 0.1 choose a different color
c3T = []
c3e = []
c5T = []
c5e = []
c6T = []
c6e = []
c7T = []
c7e = []
for i in range(0, len(ebv_o)):
    if (ebv_o[i] <= 0.05) :
        c3T.append(teff_sc70o[i])
        c3e.append(ew_sc70on[i])
    if (ebv_o[i] > 0.05) and (ebv_o[i] < 0.15):
        c5T.append(teff_sc70o[i])
        c5e.append(ew_sc70on[i])
    elif ebv_o[i] > 0.14 and (ebv_o[i] < 0.25):
        c6T.append(teff_sc70o[i])
        c6e.append(ew_sc70on[i])
    elif ebv_o[i] >= 0.25:
        c7T.append(teff_sc70o[i])
        c7e.append(ew_sc70on[i])
#pyplot.plot(teff_cmfgen, ew_cmfgen, 'bo-', mfc='None')
#pyplot.plot(teff_cmfgen, ew_tlusty, 'rs-', mfc='None')

c3 = pyplot.plot(c3T, c3e, 'r>')
c5 = pyplot.plot(c5T, c5e, 'm+', ms=10)
c6 = pyplot.plot(c6T, c6e, 'g1', ms=10)
c7 = pyplot.plot(c7T, c7e, 'kd')
    
curves = []
curves.append(c1)
curves.append(c2)
curves.append(c3)
curves.append(c4)
curves.append(c5)
curves.append(c6)
curves.append(c7)
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
# make the leyend box
pyplot.text(43300, -22.5, ' ')
leg = pyplot.figure(1).legend(curves, data_type, bbox_to_anchor=(0.87, 0.4), labelspacing=0.2)
for t in leg.get_texts():
    t.set_fontsize(17)     
#save the plot
epsfile = os.path.join(path_plots, plot_name)
pyplot.savefig(epsfile)
pyplot.show()    

print("Do you want to save this plot?  [y/N , meaning NO is default... :P ]")
save_plt = raw_input()
pyplot.ioff()
if save_plt == 'n' or save_plt == '':
    print('Plot not saved...   :( ')
    os.remove(epsfile)
else:
    print('Plot %s was saved.' % plot_name)


print ('Done!')