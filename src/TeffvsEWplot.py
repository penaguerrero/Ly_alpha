import numpy
import os
from matplotlib import pyplot
from pylab import *
from matplotlib.ticker import MaxNLocator


'''
This code makes a nice plot of EW vs. effective temperature. It uses the input table used for Starburst999 to determine Ly-alpha EWs.
'''

sb99table = '../Starburst99/lymanalpha.txt'
path_plots = '../results/plots/'
plot_name = 'EWvsTeff.eps'

logg = [2.25, 3.25, 4.25]
temps = [15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 32.5, 35., 37.5, 40., 42.5, 45., 48.5]
teff = []
for t in temps:
    Te = t*1000.
    teff.append(Te)

data = numpy.loadtxt(sb99table, skiprows=3, usecols=(2, 6, 10), unpack=True)

# making all fonts biger
font = {#'family' : 'Vera Sans',
        'weight' : 'regular',
        'size'   : 17}
matplotlib.rc('font', **font)

# plotting
fig = pyplot.figure(1, figsize=(10, 10))
xlab = pyplot.xlabel('Effective Temperature [K]')
ylab = pyplot.ylabel('Lyman-alpha Equivalent Width [$\AA$]')
# add some space between labels and axis
xlab.set_position((0.5, 0.02))
ylab.set_position((0.9, 0.5))
# change the symbols according to  cmfgen or tlusty
for i in range(0, len(teff)):
    if teff[i] < 30000:
        pyplot.plot(teff[i], data[0][i], 'bo-')
        pyplot.plot(teff[i], data[1][i], 'rs-')
        pyplot.plot(teff[i], data[2][i], 'g^-')
    else:
        pyplot.plot(teff[i], data[0][i], 'bo-', mfc='None')
        pyplot.plot(teff[i], data[1][i], 'rs-', mfc='None')
        pyplot.plot(teff[i], data[2][i], 'g^-', mfc='None')

c1 = pyplot.plot(teff, data[0], 'b')
c2 = pyplot.plot(teff, data[1], 'r')
c3 = pyplot.plot(teff, data[2], 'g')
        
curves = []
curves.append(c1)
curves.append(c2)
curves.append(c3)
# remove the firt tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
# make the leyend box
pyplot.text(43300, -22.5, 'log g')
leg = pyplot.figure(1).legend(curves, logg, bbox_to_anchor=(0.87, 0.35), labelspacing=0.2)
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