import os
import glob
import numpy
import table
import matplotlib
matplotlib.use("macosx")
#from pprint import pprint
#from pylab import *
#from random import randrange
from matplotlib.ticker import MaxNLocator
from science import spectrum 
from matplotlib import pyplot

'''
This code makes the plots of EQW vs log(age) using the output files of Starburst99.

In this version of the plots and simply add a zoom-in in the stellar part. 

'''

def read_S99file(txt_file):
    data = numpy.loadtxt(txt_file, dtype='float', skiprows=7, unpack=True)
    return data

# Path to the files
path_results = '../results/SB99_round4/'
#path_results = '../results/SB99_round4/test/OLDversionSB99/txt_files/'
#path_results = '../results/SB99_round4/test/OLDsbNEWlya/txt_files/'
#path_results = '../results/SB99_round4/test/OLDsbMODlya/txt_files/'
#path_results = '../results/SB99_round4/test/NEWversionSB99/txt_files/'
path_plots = '../results/plots/'
S99 = glob.glob('../Starburst99/round4/*')
#S99 = glob.glob('../results/SB99_round4/test/OLDversionSB99/*')
#S99 = glob.glob('../results/SB99_round4/test/OLDsbNEWlya/*')
#S99 = glob.glob('../results/SB99_round4/test/OLDsbMODlya/*')
#S99 = glob.glob('../results/SB99_round4/test/NEWversionSB99/*')

constant_ste = []
constant_neb = []
constant_tot = []
fixed_ste = []
fixed_neb = []
fixed_tot = []
i = 0
IMFs = ['0.5',
        #'0.8',
        #'1.1',
        #'1.2',
        #'1.3',
        #'1.4',
        #'1.5',
        #'1.6',
        '1.7',
        #'1.8',
        #'1.9', 
        #'2.0', 
        #'2.1',
        #'2.2',
        '2.3',
        #'2.5',
        '2.6']

total_chosenfiles = []

colors = ['blue', 'red', 'green', 'magenta', 'k', '#594A2B', '#994A2B', '#2B494A', '#7E1BE0', '#1BE0D3', '#E0531B', '#94E01B', '#9868C4']
file_count = 13

'''
colors = []
def generate_colors(count):
    ###return value is list containing count number of colors ###
    colors = []
    for _ in range(file_count):
        colors.append("#%s" % "".join([hex(randrange(1, 254))[2:] for _ in range(3)]))
    return colors
colors = generate_colors(file_count)
'''

# This loop only choses the selected IMFs to annalyze the files and create the plots
for each_folder in S99:
    for j in range (0, len(IMFs)):
        if 'imf'+IMFs[j] in each_folder:
            total_chosenfiles.append(each_folder)
    
for foldr in total_chosenfiles:
    print("%d chosen files" % len(total_chosenfiles))
    print("File %i of %i" % (i+1, len(total_chosenfiles))) # for only a few IMFs
    print foldr
    # No need to increase i because it is being increased in the next line
    file_name = os.path.basename(foldr)
    # THIS IS FOR THE NEBULAR COMPONENT
    # Get luminosity of Lyman-alpha
    logLHalpha_file = '/Output/'+file_name+'.ewidth'
    txt_file = foldr+logLHalpha_file
    data = read_S99file(txt_file)
    lumLyA = (10**data[2]) * 8.9
    time_and_lumLyA = numpy.array([data[0], lumLyA])
    
    # Nebular eqw(LyA)
    logneb_contLyA_file = '/Output/'+file_name+'.spectrum'
    txt_file = foldr+logneb_contLyA_file
    data = read_S99file(txt_file)
    indexes=numpy.where( data == 1195.00 )
    age = []
    nebcontinuum =[]
    for indx in indexes[1]:
        age_at_indx = (data[0][indx])
        nebcont_at_indx = 10**(data[2][indx])
        age.append(age_at_indx)
        nebcontinuum.append(nebcont_at_indx)
    age_and_nebcontinuum = numpy.array([age, nebcontinuum])   
    
    # Finding corresponding LyA luminosities and times with their nebular parts
    table_age_and_nebcontinuum = []
    hdr = ['age_nebcontinuum', 'log age', 'LyA', 'nebcontinuum', 'neb_eqw']
    table_age_and_nebcontinuum.append(hdr)
    neb_age = []
    neb_eqw = []
    for j in range(0, len(age_and_nebcontinuum[0])):
        if age_and_nebcontinuum[0][j] in time_and_lumLyA[0]:
            lLyA = spectrum.findXinY(time_and_lumLyA[1], time_and_lumLyA[0], age_and_nebcontinuum[0][j])
            eqw = lLyA / age_and_nebcontinuum[1][j]   # dealing with anti-logs
            #eqw = lLyA - age_and_nebcont[1][j]   # dealing with logs
            neb_age.append(age_and_nebcontinuum[0][j])
            neb_eqw.append(eqw)
            temp = [repr(age_and_nebcontinuum[0][j]), repr(numpy.log10(age_and_nebcontinuum[0][j])), repr(lLyA), repr(age_and_nebcontinuum[1][j]), repr(eqw)]
            table_age_and_nebcontinuum.append(temp)
    neb_age_and_eqw = numpy.array([neb_age, neb_eqw])
    # Saving the table into file 
    table_out = open(path_results+'age_and_nebcontinuum_'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_age_and_nebcontinuum)
    table_out.close()
    print('file  %s   is now written!' % (path_results+'age_and_nebcontinuum_'+file_name+'.txt'))
    
    # Stellar eqw(LyA)
    ste_eqwLyA_file = '/Output/'+file_name+'.irfeature'
    txt_file = foldr+ste_eqwLyA_file
    data = read_S99file(txt_file)
    age = []
    eqwlya = []
    for j in range(0, len(data[0])):
        if data[0][j] in neb_age_and_eqw[0]:
            #age.append(data[0][j]/1e6)            # age in mega years    
            #age.append(numpy.log10(data[0][j]))   # age in log
            age.append(data[0][j])                 # age in scientific notation
            #eqwlya.append(data[7][j]*(-1))       # the negative sign is because EQWs are in absorption (ONLY used for initial lyas table)
            eqwlya.append(data[7][j])       # no need to use the negative sign (EQWs are in absorption), data is fixed
        
    # length check
    while len(eqwlya) > len(neb_age_and_eqw[1]):
        print ('Fixing dimensions of arrays so that neb_age_and_eqw = ste_age_and_eqw')
        eqwlya.pop(len(eqwlya)-1)
        age.pop(len(age)-1)
    #print('Fixed, lengths are equal:', len(eqwlya), len(neb_age_and_eqw[1]))
    ste_age_and_eqw = numpy.array([age, eqwlya])
    
    # Total EQW = nebular + stellar
    tot_eqw = []
    tot_age = []
    for j in range(0, len(ste_age_and_eqw[1])):
        tot = ste_age_and_eqw[1][j] + neb_age_and_eqw[1][j]
        tot_eqw.append(tot)
        tot_age.append(numpy.log10(ste_age_and_eqw[0][j]))
    tot_age_and_eqw = numpy.array([tot_age, tot_eqw])
        
    # Divide continuous from instantaneous bursts
    pair_ste = [tot_age, ste_age_and_eqw[1]]
    pair_neb = [tot_age, neb_age_and_eqw[1]]
    word = 'fixed'
    if word in foldr:
        fixed_ste.append(pair_ste)
        fixed_neb.append(pair_neb)
        fixed_tot.append(tot_age_and_eqw)
    else:
        constant_ste.append(pair_ste)
        constant_neb.append(pair_neb)
        constant_tot.append(tot_age_and_eqw)
    
    table_time_and_lumLyA = []
    hdr = ['time', 'lum LyA']#, 'age_nebcontinuum', 'nebcontinuum']
    table_time_and_lumLyA.append(hdr)
    for j in range(0, len(time_and_lumLyA[0])):
        temp = [repr(time_and_lumLyA[0][j]), repr(time_and_lumLyA[1][j])]#, repr(age_and_nebcontinuum[0][j]), repr(age_and_nebcontinuum[1][j])]
        table_time_and_lumLyA.append(temp)
    #table.pprint_table(sys.stdout, table_time_and_lumLyA)
    # Saving the table into file 
    table_out = open(path_results+'time_and_lumLyA_'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_time_and_lumLyA)
    table_out.close()
    print('file  %s   is now written!' % (path_results+'time_and_lumLyA_'+file_name+'.txt'))
        
    table_s99LyAeqw = [] 
    hdr = ['ste_age', 'ste_eqw', 'neb_age', 'neb_eqw', 'log age', 'tot_eqw']
    table_s99LyAeqw.append(hdr)
    for j in range(0, len(pair_ste[0])):
        temp = [repr(ste_age_and_eqw[0][j]), repr(ste_age_and_eqw[1][j]), repr(neb_age_and_eqw[0][j]), repr(neb_age_and_eqw[1][j]), repr(tot_age_and_eqw[0][j]), repr(tot_age_and_eqw[1][j])]
        table_s99LyAeqw.append(temp)
    #table.pprint_table(sys.stdout, table_s99LyAeqw)
    # Saving the table into file 
    table_out = open(path_results+'s99LyAeqw'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_s99LyAeqw)
    table_out.close()
    print('file  %s   is now written!' % (path_results+'s99LyAeqw'+file_name+'.txt'))
        
    i = i+1

# Give me the values of the EWs at the specified age and print it on the screen with a specific IMF 
# Salpeter IMF in this case, IMF = 2.3
ages = [5e6, 15e6]
EWs_file_list = glob.glob(os.path.join(path_results, 's99LyAeqw*_imf2.3.txt'))
cont_ew_ste = 0.0
cont_ew_neb = 0.0
cont_ew_tot = 0.0
fixed_ew_ste = 0.0
fixed_ew_neb = 0.0
fixed_ew_tot = 0.0
# the numbers in these list correspond to the ages list: 1st number at 1st age, ...
cont_ew_ste_list = []
cont_ew_neb_list = []
cont_ew_tot_list = []
fixed_ew_ste_list = []
fixed_ew_neb_list = []
fixed_ew_tot_list = []
for f in EWs_file_list:
    dataEWs = numpy.loadtxt(f, dtype=numpy.float64, skiprows=1, unpack=True)
    # dataEWs: 0=ste_age, 1=ste_eqw, 2=neb_age, 3=neb_eqw, 4=log age, 5=tot_eqw
    print 'This is file %s' % f
    if 'fixed' in f:
        for age in ages:
            for i in range(len(dataEWs[0])):
                if age == dataEWs[0][i]:
                    fixed_ew_ste = dataEWs[1][dataEWs[0] == age]
                    fixed_ew_neb = dataEWs[3][dataEWs[0] == age]
                    fixed_ew_tot = dataEWs[5][dataEWs[0] == age]
                else:
                    fixed_ew_ste = numpy.interp(age, dataEWs[0], dataEWs[1])
                    fixed_ew_neb = numpy.interp(age, dataEWs[0], dataEWs[3])
                    fixed_ew_tot = numpy.interp(age, dataEWs[0], dataEWs[5])
            fixed_ew_ste_list.append(fixed_ew_ste)
            fixed_ew_neb_list.append(fixed_ew_neb)
            fixed_ew_tot_list.append(fixed_ew_tot)
        print 'fixed_ew_ste_list, fixed_ew_neb_list, fixed_ew_tot_list', fixed_ew_ste_list, fixed_ew_neb_list, fixed_ew_tot_list
    elif 'continous' in f:
        for age in ages:
            for i in range(len(dataEWs[0])):
                if age == dataEWs[0][i]:
                    cont_ew_ste = dataEWs[1][dataEWs[0] == age]
                    cont_ew_neb = dataEWs[3][dataEWs[0] == age]
                    cont_ew_tot = dataEWs[5][dataEWs[0] == age]
                else:
                    cont_ew_ste = numpy.interp(age, dataEWs[0], dataEWs[1])
                    cont_ew_neb = numpy.interp(age, dataEWs[0], dataEWs[3])
                    cont_ew_tot = numpy.interp(age, dataEWs[0], dataEWs[5])
            cont_ew_ste_list.append(cont_ew_ste)
            cont_ew_neb_list.append(cont_ew_neb)
            cont_ew_tot_list.append(cont_ew_tot)
        print 'cont_ew_ste_list, cont_ew_neb_list, cont_ew_tot_list', cont_ew_ste_list, cont_ew_neb_list, cont_ew_tot_list

    
lower = 6.3
upper = 8
font = {#'family' : 'Vera Sans',
        'weight' : 'regular',
        'size'   : 15}
# OR ONE CAN USE FOR INDIVIDUAL AXIS:
#f1.tick_params(axis='x', labelsize=17)
#f1.tick_params(axis='y', labelsize=17)
matplotlib.rc('font', **font)

figs = pyplot.figure(1, figsize=(8, 12))
figs.subplots_adjust(hspace=0.10)
#figs.suptitle()
# Set common labels
figs.text(0.5, 0.06, 'log age [yr]', ha='center', va='center')
figs.text(0.035, 0.5, 'Ly-alpha EW [$\AA$]', ha='center', va='center', rotation='vertical')

f1 = figs.add_subplot(211)
f1.set_title('Instantaneous SFR')
#[label.set_visible(False) for label in f1.get_xticklabels()]
#step(f1.get_xticklabels(), visible=False)    # this is the same as the previous line but for some reason it is not working
f1.set_xlim(lower, upper)
f1.set_ylim(-60, 250)
f1.text(6.5, -30, 'stellar component')   
f1.text(6.33, 70, 'total=')   
f1.text(6.33, 50, 'nebular')   
f1.text(6.4, 35, '+')  
f1.text(6.36, 20, 'stellar')
f1.text(7.2, 10, 'nebular component')  
f1.text(7.6, 207 , 'IMF exp')
# adding the enhancement of the stellar component  
f2 = figs.add_subplot(212)
#f2.set_xlim(lower, upper) 
f2.set_xlim(6.5, 7.3)
f2.set_ylim(-20, -3.)
f2.text(6.75, -7, 'stellar component')
f2.text(7.0, -13.3, 'IMF exp = 0.5', rotation=-21)
f2.text(6.9, -14.95, 'IMF exp = 2.6', rotation=-27)
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

figs = pyplot.figure(2, figsize=(8, 12)) 
figs.subplots_adjust(hspace=0.10)
f3 = pyplot.subplot(211)
f3.set_title('Constant SFR') 
# Set common labels
figs.text(0.5, 0.06, 'log age [yr]', ha='center', va='center')
figs.text(0.035, 0.5, 'Ly-alpha EW [$\AA$]', ha='center', va='center', rotation='vertical')
f3.set_xlim(lower, upper)
#f3.set_ylim(-15, 170)
f3.text(6.6, 6, 'stellar component')   
f3.text(7.27, 37, 'total = nebular + stellar')   
f3.text(6.8, 210, 'nebular component')   
f3.text(7.65, 320, 'IMF exp')  
# adding the enhancement of the stellar component  
f4 = figs.add_subplot(212)#, sharex=f3)
#f4.set_xlim(lower, upper) 
f4.set_xlim(6.5, 7.3)
f4.set_ylim(-8, -1.0)
f4.text(6.67, -1.7, 'stellar component')   
f4.text(7.1, -2.7, 'IMF exp = 0.5')   
f4.text(7.0, -4.3, 'IMF exp = 1.7', rotation=-17)   
f4.text(6.9, -5.3, 'IMF exp = 2.3', rotation=-27)   
f4.text(6.8, -5.5, 'IMF exp = 2.6', rotation=-35)   
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

# Ploting
curves_inst = []
curves_cont =[]
line_styles = [':', '--', '.-', '-']
for i in range(0, len(constant_tot)):
    
    curve, = f1.plot(fixed_ste[i][0], fixed_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve, = f1.plot(fixed_neb[i][0], fixed_neb[i][1], colors[i])
    curve.set_linestyle(':')
    curve_i, = f1.plot(fixed_tot[i][0], fixed_tot[i][1], colors[i])
    curves_inst.append(curve_i)
    # the stellar enhancement
    curve, = f2.plot(fixed_ste[i][0], fixed_ste[i][1], colors[i])
    curve.set_linestyle('--')
    
    #pyplot.figure(2)
    curve, = f3.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve, = f3.plot(constant_neb[i][0], constant_neb[i][1], colors[i])
    curve.set_linestyle(':')
    curve_c, = f3.plot(constant_tot[i][0], constant_tot[i][1], colors[i])
    curves_cont.append(curve_c)
    # the stellar enhancement
    curve, = f4.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')

leg1 = pyplot.figure(1).legend(curves_inst, IMFs, bbox_to_anchor=(0.85, 0.85))  
leg2 = pyplot.figure(2).legend(curves_cont, IMFs, bbox_to_anchor=(0.87, 0.875), labelspacing=0.2)
for t in leg1.get_texts():
    t.set_fontsize(12)     
for t in leg2.get_texts():
    t.set_fontsize(12)     
pyplot.show()

'''
# Saving the plots
pyplot.figure(1)
epsfile = os.path.join(path_plots, "Inst_logAge_NEW4.eps")
pyplot.savefig(epsfile)#, dpi=(100))

pyplot.figure(2)
epsfile = os.path.join(path_plots, "Const_logAge_NEW4.eps")
pyplot.savefig(epsfile)
pyplot.ion()
'''

print('Done!')


