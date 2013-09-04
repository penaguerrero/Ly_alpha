import numpy
import os
from matplotlib import pyplot
from pylab import *
from matplotlib.ticker import MaxNLocator


'''
This code makes a nice plot of EW vs. effective temperature. It compares our EW determinations to those of Savage & Code (1970)
and Savage & Panek (1974).

It also determines the value of the star's radius and log g.

'''

table_OBew = '../results/table_OB_desired_eqws.txt'
table_obsInfo = '../results/table_obsInfo.txt'
path_plots = '../results/plots/'
plot_name = 'comparison_with_obs_EWvsTeff3.eps'


def char_consecutive(s, find):
    """Count consecutive occurrences of a character in a string
    s: string
    find: character to find

    return: count of characters
    """
    found = False
    consecutive = 0
    for char in s:
        if char == find:
            found = True
        else:
            found = False

        if found:
            consecutive += 1
    return consecutive

def find_infO(Ostar):
    # look what is the spectral type
    for i in range(0, len(Osp)):
        if Osp[i] in Ostar:
            # now find the luminosity type
            # first determine how many times I and V appear 
            chars = char_consecutive(Ostar, 'I')
            condition = 0
            if chars == 2:
                condition = 2
            elif chars == 3:
                condition = 3
            elif chars == 1:
                chars = char_consecutive(Ostar, 'V') 
                if chars == 1:
                    condition = 4
                elif condition == 0:
                    condition = 1
            elif chars == 0:
                chars = char_consecutive(Ostar, 'V')
                if chars == 1:
                    condition = 5
            # now do the corresponding operation
            for j in range(0, len(lumtype)):
                if lumtype[j] in Ostar:
                    if condition == 1:
                        l = deltaOlum[i]*4 + Olum[i]
                        m = Omass[i] - deltaOmass[i]*4
                        t = Oteff[i] - deltaOteff[i]*4
                    elif condition == 2:
                        l = deltaOlum[i]*3 + Olum[i]
                        m = Omass[i] - deltaOmass[i]*3
                        t = Oteff[i] - deltaOteff[i]*3
                    elif condition == 3:
                        l = deltaOlum[i]*2 + Olum[i]
                        m = Omass[i] - deltaOmass[i]*2
                        t = Oteff[i] - deltaOteff[i]*2
                    elif condition == 4:
                        l = deltaOlum[i] + Olum[i]
                        m = Omass[i] - deltaOmass[i]
                        t = Oteff[i] - deltaOteff[i]
                    elif condition == 5:
                        l = Olum[i]
                        m = Omass[i]
                        t = Oteff[i]
    return(l, m, t)

def find_infB(Bstar):
    # look what is the spectral type
    for i in range(0, len(Bsp)):
        if Bsp[i] in Bstar:
            # now find the luminosity type
            # first determine how many times I and V appear 
            chars = char_consecutive(Bstar, 'I')
            condition = 0
            if chars == 2:
                condition = 2
            elif chars == 3:
                condition = 3
            elif chars == 1:
                chars = char_consecutive(Bstar, 'V') 
                if chars == 1:
                    condition = 4
                elif condition == 0:
                    condition = 1
            elif chars == 0:
                chars = char_consecutive(Bstar, 'V')
                if chars == 1:
                    condition = 5
            # now do the corresponding operation
            for j in range(0, len(lumtype)):
                if lumtype[j] in Bstar:
                    if condition == 1:
                        l = deltaBlum[i]*4 + Blum[i]
                        m = Bmass[i] + deltaBmass[i]*4
                        t = Bteff[i] - deltaBteff[i]*4
                    elif condition == 2:
                        l = deltaBlum[i]*3 + Blum[i]
                        m = Bmass[i] + deltaBmass[i]*3
                        t = Bteff[i] - deltaBteff[i]*3
                    elif condition == 3:
                        l = deltaBlum[i]*2 + Blum[i]
                        m = Bmass[i] + deltaBmass[i]*2
                        t = Bteff[i] - deltaBteff[i]*2
                    elif condition == 4:
                        l = deltaBlum[i] + Blum[i]
                        m = Bmass[i] + deltaBmass[i]
                        t = Bteff[i] - deltaBteff[i]
                    elif condition == 5:
                        l = Blum[i]
                        m = Bmass[i]
                        t = Bteff[i]
    return(l, m, t)


# There are 3 types of data for this plot: CMFGEN, TLUSTY, and OBSERVATIONS
data_type = ["TW(CMFGEN)", "TW(TLUSTY)", "SC70 E(B-V)<0.05", "SC70 0.05<E(B-V)<0.15", "SC70 0.15<E(B-V)<0.25", "SC70 E(B-V)>0.25", "SP74"]

# This is the data from Svage & Code (70) = sc70
# Spectral types
st_sc70 = ['B3IV', 'B1I', 'B0.5V', 'B3V', 'B3V', 'B0.5IV', 'B0.5V', 'B1V', 'B2III', 'B7III', 'O9.5II', 'B0V', 'B0IV', 'O9.5V', 'O9III',
           'B0I', 'B2IV', 'O9.5V', 'O9.5I', 'B0.5I', 'B1II', 'B2II', 'B7V', 'B1V', 'B3V', 'B1II', 'B2V', 'B1.5V', 'B5V',
           'B1V', 'B2V', 'B5IV', 'B3V', 'B2IV', 'B5V', 'B2.5V', 'B2.5V', 'B8I', 'B1V', 'B0V', 'B1V', 'B1III', 'O9.5V', 
           'B1.5V', 'B1V', 'B2IV', 'B3IV', 'B5V']

'''
# temparatures obtained from literature
teff_sc70 = [15170., 24400, 30000., 18700., 17620., 26850., 18200., 25290, 21580., 13160., 31000., 33400., 30000., 29300., 32000., 
             24820., 21480., 31270., 29910., 26390., 25180., 20990., 12250., 23930., 17030., 25000., 23000., 20500., 15780., 23100., 
             22900., 17900., 16860., 26200., 16980., 19690., 21380., 21380., 25350., 31000., 28000., 24200., 33000., 22500., 26300., 
             24720., 20180., 14140.]
'''
ew_sc70 = [28., 25., 14., 29., 25., 14., 20., 19., 15., 35., 20., 13., 31., 33., 17., 24., 21., 16., 25., 22., 10., 12., 43., 11., 22., 12.,
           20., 15., 36., 14., 17., 26., 22., 14., 26., 22., 30., 12., 18., 24., 33., 31., 26., 20., 12., 16., 29., 33.]
# E(B-V)
ebv = [0.04, 0.31, 0.10, 0.07, 0.02, 0.01, 0.09, 0.05, 0.01, 0.01, 0.09, 0.03, 0.13, 0.36, #0.21, 
       0.07, 0.05, 0.06, 0.06, 0.06, 0.04, 0.02, 0.00, 0.00, 0.02, 0.01, 0.03, 0.04, 0.04, 0.01, 
       0.05, 0.03, 0.01, 0.02, 0.01, 0.00, 0.05, 0.17, 0.21, 0.07, 0.19, 0.22, 0.40, 0.32, 0.02, 0.06, 
       0.03, 0.05, 0.01]

# This is the data from Svage & Panek (74) = sp74
'''
teff_sp74 = [14400., 14500., 16470., 12800., 15650., 17940., 12100., 13160., 14200., 18000., 18260., 17200., 12250., 17030., 19525., 18280., 
             11500., 15780., 25000., 17900., 16860., 19690., 23440., 14890., 13400., 17620., 18890., 17850., 14140.]
'''
ew_sp74 = [45.0, 21.0, 25.0, 65., 31., 22., 30., 42., 45., 21., 20., 19., 55., 25., 10., 18., 23., 45., 20., 30., 
           21., 19., 11., 35., 44., 21., 18., 17., 45.]
# Spectral Types 
st_sp74 = ['B9II', 'B3V', 'B3V', 'B8V', 'B5III', 'B3V', 'B8I', 'B7III', 'B7IV', 'B3IV', 'B2.5IV', 'B3IV', 'B7V', 'B3V', 'B2.5IV',
           'B2.5IV', 'B7III', 'B5V', 'B5IV', 'B3V', 'B2.5V', 'B2.5IV', 'B5IV', 'B6III', 'B3IV', 'B3IV', 'B2.5V', 'B7IV']

# From table 3.1 of Conti, P. S., Crowther, P. A., & Leitherer, C. 2008, "From luminous hot stars to starburst galaxies", cambridge 
#university press, Cambridge, UK
# for Spectral type V
Osp =['8', '9', '9.5']
Oteff = [35000., 33000., 31500.]
deltaOteff = [500., 375., 625.]
Omass = [31, 25, 23]
deltaOmass = [0.25, 0.25, 0.75]
Olum = [5.10, 4.91, 4.78]
deltaOlum = [0.1325, 0.1675, 0.2]
# for Spectral type V
Bsp = ["0", "0.5", "1", "1.5", "2", "2.5", "3", "5", "6", "7", "8", "9"]
Bteff = [29500, 28000, 26000, 24000, 21000, 19000, 17500, 15400, 14400, 13300, 12300, 11300]
deltaBteff = [500., 500., 1125., 875., 625., 625., 500., 475., 400., 300., 225., 150.]
Bmass = [19, 18, 14, 12, 9, 7.5, 6, 5, 4.5, 4, 4, 3.5]
deltaBmass = [0.25, 0.5, 1., 1.5, 2., 1.875, 2., 1.75, 1.625, 1.5, 1.25, 1.125]
Blum = [4.61, 4.48, 4.24, 3.96, 3.59, 3.34, 3.06, 2.79, 2.64, 2.50, 2.35, 2.20]
deltaBlum = [0.245, 0.275, 0.3125, 0.38, 0.46, 0.4975, 0.5525, 0.45825, 0.1025, 0.5475, 0.5325, 0.5175]

# interpolate to thet the luminosities
lumtype = ['I', 'II', 'III', 'IV', 'V']
lum_sc70 = []
mass_sc70 = []
teff_sc70 = []
lum_sp74 = []
mass_sp74 = []
teff_sp74 = []

for i in range(0, len(st_sc70)):
    if 'O' in st_sc70[i]:
        l, m, t = find_infO(st_sc70[i])
    elif 'B' in st_sc70[i]:
        l, m, t = find_infB(st_sc70[i])
    lum_sc70.append(l)
    mass_sc70.append(m)
    teff_sc70.append(t)

for i in range(0, len(st_sp74)):
    if 'O' in st_sp74[i]:
        l, m, t = find_infO(st_sp74[i])
    elif 'B' in st_sp74[i]:
        l, m, t = find_infB(st_sp74[i])
    lum_sp74.append(l)
    mass_sp74.append(m)
    teff_sp74.append(t)

'''
teff_sc70 = [17000, 21500, 28000, 17500, 17500, 27500, 28000, 26000, 19750, 12700, 29625, 29500, 29000, 31500, 30250,
             27500, 19125, 31500, 29000, 26000, 22625, 19125, 13300, 26000, 17500, 22625, 21000, 24000, 15400, 26000,
             21000, 14925, 17500, 20375, 15400, 19000, 19000, 11400, 26000, 29500, 26000, 23750, 31500, 24000, 26000, 
             20375, 17000, 15400]
teff_sp74 = [10850, 17500, 17500, 12300, 14450, 17500, 11400, 12700, 13000, 17000, 18500, 17000, 13300, 17500, 18375,
             18375, 12700, 15400, 17500, 15400, 17500, 19000, 18375, 14975, 13600, 17000, 17000, 19000, 13000]
'''
# Using this info calculate the value of R/Rsun
R_sc70 = []
g_sc70 = []
R_sp74 = []
g_sp74 = []
for i in range(0, len(lum_sc70)):
    logR = 0.5*numpy.log10(lum_sc70[i]) - 2*numpy.log10(teff_sc70[i]) + 7.52
    R = 10**(logR)
    logg = 4.44 + numpy.log10(mass_sc70[i]/(R**2))
    R_sc70.append(R)    # in solar units
    g_sc70.append(logg)
    print(st_sc70[i], lum_sc70[i], teff_sc70[i], R, logg)

for i in range(0, len(lum_sp74)):
    logR = 0.5*numpy.log10(lum_sp74[i]) - 2*numpy.log10(teff_sp74[i]) + 7.52
    R = 10**(logR)
    logg = 4.44 + numpy.log10(mass_sp74[i]/(R**2))
    R_sp74.append(R)    # in solar units
    g_sp74.append(logg)

# Write this in a table
txt = open(table_obsInfo, 'w+')
txt.write('Savage & Code 1970\n')
txt.write('SpTyp    Teff    L/Ls    M/Ms    R/Rsun    log(g)\n')
for i in range(0, len(R_sc70)):
    txt.write('%s    %i     %.3f     %.3f    %.3f    %.3f\n' % (st_sc70[i], teff_sc70[i], lum_sc70[i], mass_sc70[i], R_sc70[i], g_sc70[i]))
txt.write('Savage & Panek 1974\n')
for i in range(0, len(R_sp74)):
    txt.write('%s    %i     %.3f     %.3f    %.3f    %.3f\n' % (st_sp74[i], teff_sp74[i], lum_sp74[i], mass_sp74[i], R_sp74[i],g_sp74[i]))
txt.close()



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