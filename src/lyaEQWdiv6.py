import glob
import numpy
import matplotlib
import string
import os
import table
matplotlib.use("macosx")
from matplotlib import pyplot
from random import randrange


'''
This code takes the written eqws file with title "s99LyAeqw*+IMFs[i]+.txt" and uses it to make the plots for the paper.
'''


constant_ste = []
constant_neb = []
constant_tot = []
fixed_ste = []
fixed_neb = []
fixed_tot = []
tot_lyas_files = []

# Path where the text files are
files_in_results = glob.glob('../results/s99LyAeqw*_imf*.txt')
# Path to the resulting plots
path_results = '../results/plots/'


IMFs = [#'0.5',
        '1.3',
        '1.4',
        '1.5',
        '1.6',
        '1.7',
        '1.8',
        '1.9', 
        '2.0', 
        '2.1',
        '2.2',
        '2.3']

total_chosenfiles = []

colors = ['blue', 'red', 'green', 'magenta', '#994A2B', '#2B494A', '#7E1BE0', '#1BE0D3', '#E0531B', '#94E01B', '#9868C4']
file_count = 11

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

# This loop makes sure that we are only picking the right files with the chosen IMFs
ignore_this_file = False
for text_file in files_in_results:
    if 'cutoff' in text_file:
        ignore_this_file = True
    else:
        file_name = string.split(text_file, sep='imf')
        exponent = file_name[1].replace(".txt", "")
        for imf in IMFs:
            if exponent == imf:
                if 'continous' in text_file:
                    tot_lyas_files.append(text_file)
                elif 'fixed' in text_file:
                    tot_lyas_files.append(text_file)
    
files2plot = len(tot_lyas_files)
print('Total files to obtain information for the plot: %i' % files2plot)

# lists to make the EQWs table
table_fixed = []
table_const = []
table_yrs = []
# This loop actually makes the plots
i = 1
for txt_file in tot_lyas_files:
    print('Working with %s' % txt_file)
    print('File %i of %i' % (i, files2plot))
    data = numpy.loadtxt(txt_file, dtype=float, skiprows=1, unpack=True)
    #data = ste_age, ste_eqw, neb_age, neb_eqw, log_age, tot_eqw
    #        data[0], data[1], data[2], data[3], data[4], data[5]
    # Divide continuous from instantaneous bursts
    #print(data)
    pair_ste = [data[4], data[1]]
    pair_neb = [data[4], data[3]]
    tot_age_and_eqw = [data[4], data[5]]
    word = 'fixed'
    if word in txt_file:
        # Since time is the same for both files
        for t in data[0]:
            mill_yrs = t / 1e6
            table_yrs.append(t)
        fixed_ste.append(pair_ste)
        fixed_neb.append(pair_neb)
        fixed_tot.append(tot_age_and_eqw)
        t = zip(data[1], data[3], data[5])
        table_fixed.append(t)
    else:
        constant_ste.append(pair_ste)
        constant_neb.append(pair_neb)
        constant_tot.append(tot_age_and_eqw)
        t = zip(data[1], data[3], data[5])
        table_const.append(t)
    i = i + 1
    
# The pretty table of the data
table_out = []
table_out.append(['time [10^6 yr]', 'log time', '','Instantaneous SFR', '', '', 'Constant SFR', ''])
table_out.append(['', '','stellar EQW', 'nebular EQW', 'total=nebular+stellar', 'stellar EQW', 'nebular EQW', 'total=nebular+stellar'])
table_ste_eqwf = []
table_neb_eqwf = []
table_tot_eqwf = []
table_ste_eqwc = []
table_neb_eqwc = []
table_tot_eqwc = []
for j in range(0, len(IMFs)):
    tf1, tf2, tf3 = zip(*table_fixed[j])
    tc1, tc2, tc3 = zip(*table_const[j])
    for i in range(0, len(tc1)):
        table_ste_eqwf.append(tf1[i])
        table_neb_eqwf.append(tf2[i])
        table_tot_eqwf.append(tf3[i])
        table_ste_eqwc.append(tc1[i])
        table_neb_eqwc.append(tc2[i])
        table_tot_eqwc.append(tc3[i])

print(len(table_yrs), len(table_ste_eqwf), len(table_neb_eqwf), len(table_tot_eqwf), len(table_ste_eqwc), len(table_neb_eqwc), len(table_tot_eqwc))
for i in range (0, len(table_yrs)):
    t = [repr(table_yrs[i]), repr(table_ste_eqwf[i]), repr(table_neb_eqwf[i]), repr(table_tot_eqwf[i]),
         repr(table_ste_eqwc[i]), repr(table_neb_eqwc[i]), repr(table_tot_eqwc[i])]
    table_out.append(t)
#print(table_out)
f = open('../results/SB99eqws.txt', 'w+')
#table.pprint_table(f, table_out)
for i in range (0, len(table_yrs)):
    print >> f, table_out[i]
f.close()

            
lower = 6.3
upper = 8
pyplot.figure(1, figsize=(10, 9))
pyplot.title('Instantaneous SFR')
pyplot.xlabel('log age [yr]')
pyplot.ylabel('Ly-alpha EQW [$\AA$]')
pyplot.xlim(lower, upper)
pyplot.ylim(-15, 225)
pyplot.text(6.33, -9, 'stellar component')   
pyplot.text(6.35, 100, 'total=')   
pyplot.text(6.36, 91, 'nebular')   
pyplot.text(6.41, 85, '+')  
pyplot.text(6.36, 78, 'stellar')
pyplot.text(6.6, 100, 'nebular component')  
pyplot.text(7.7, 209 , 'IMF exp')  
 
pyplot.figure(2, figsize=(10,9))
pyplot.title('Constant SFR') 
pyplot.xlabel('log age [yr]')
pyplot.ylabel('Ly-alpha EQW [$\AA$]')
pyplot.xlim(lower, upper)
#pyplot.ylim(-15, 170)
pyplot.text(6.4, 6, 'stellar component')   
pyplot.text(7.3, 72, 'total=nebular+stellar')   
pyplot.text(7.3, 165, 'nebular component')   
pyplot.text(7.75, 286, 'IMF exp')   

# Ploting
curves_inst = []
curves_cont =[]
for i in range(0, len(constant_tot)):
    pyplot.figure(1)
    curve, = pyplot.plot(fixed_ste[i][0], fixed_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve_i, = pyplot.plot(fixed_neb[i][0], fixed_neb[i][1], colors[i])
    curves_inst.append(curve_i)
    curve, = pyplot.plot(fixed_tot[i][0], fixed_tot[i][1], colors[i])
    curve.set_linestyle(':')
    pyplot.figure(2)
    curve, = pyplot.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve_c, = pyplot.plot(constant_neb[i][0], constant_neb[i][1], colors[i])
    curves_cont.append(curve_c)
    curve, = pyplot.plot(constant_tot[i][0], constant_tot[i][1], colors[i])
    curve.set_linestyle(':')

leg1 = pyplot.figure(1).legend(curves_inst, IMFs, bbox_to_anchor=(0.85, 0.85))  
leg2 = pyplot.figure(2).legend(curves_cont, IMFs, bbox_to_anchor=(0.87, 0.87), labelspacing=0.2)
for t in leg1.get_texts():
    t.set_fontsize(10)     
for t in leg2.get_texts():
    t.set_fontsize(10)     
pyplot.show()

'''
# Saving the plots
pyplot.figure(1)
epsfile = os.path.join(path_results, "Inst_log4.eps")
pyplot.savefig(epsfile)#, dpi=(100))

pyplot.figure(2)
epsfile = os.path.join(path_results, "Const_log4.eps")
pyplot.ion()
pyplot.savefig(epsfile)
'''


