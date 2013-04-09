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
#files_in_results = glob.glob('../results/SB99_round2and3/s99LyAeqw*_imf*.txt')
files_in_results = glob.glob('../results/SB99_round4/s99LyAeqw*_imf*.txt')
######
# tests of OLD and NEW versions of lyalpha EQWs:
# OLD
#files_in_results = glob.glob('../results/SB99_round4/test/OLDversionSB99/txt_files/s99LyAeqw*_imf*.txt')
# MOD
#files_in_results = glob.glob('../results/SB99_round4/test/OLDsbMODlya/txt_files/s99LyAeqw*_imf*.txt')
# NEW but only a part used
#files_in_results = glob.glob('../results/SB99_round4/test/OLDsbNEWlya/txt_files/s99LyAeqw*_imf*.txt')
# NEW
#files_in_results = glob.glob('../results/SB99_round4/test/NEWversionSB99/txt_files/s99LyAeqw*_imf*.txt')

# subtitle in the plot to know which lyman-alpha file was used
lya_used = 'NEW_lyas'

# Path to the resulting plots
path_results = '../results/plots/'


IMFs = [#'0.5',
        '1.2',
        '1.3',
        '1.4',#***
        '1.5',
        '1.6',#***
        '1.7',
        '1.8',#***
        '1.9', 
        '2.0', 
        '2.1',
        '2.2',
        '2.3',
        '2.5']

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
table_fixed5 = []
table_const5 = []
table_fixed10 = []
table_const10 = []
table_fixed15 = []
table_const15 = []
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
        fixed_ste.append(pair_ste)
        fixed_neb.append(pair_neb)
        fixed_tot.append(tot_age_and_eqw)
        # now find the eqws at 5=6.7, 10=7.0 and 15x10^6=7.2 yrs
        for j in range(0, len(data[4])):
            if (data[4][j] >= 6.7) and (data[4][j] < 6.79):
                table_fixed5.append([data[1][j], data[3][j], data[5][j]])
            #if (data[4][j] >= 7.0) and (data[4][j] < 7.05):
            #    table_fixed10.append([data[1][j], data[3][j], data[5][j]])
            if (data[4][j] >= 7.2) and (data[4][j] < 7.21):
                table_fixed15.append([data[1][j], data[3][j], data[5][j]])
    else:
        constant_ste.append(pair_ste)
        constant_neb.append(pair_neb)
        constant_tot.append(tot_age_and_eqw)
        # now find the eqws at 5=6.7, 10=7.0 and 15x10^6=7.2 yrs
        for j in range(0, len(data[4])):
            if (data[4][j] >= 6.7) and (data[4][j] < 6.79):
                table_const5.append([data[1][j], data[3][j], data[5][j]])
            #if (data[4][j] >= 7.0) and (data[4][j] < 7.05):
            #    table_const10.append([data[1][j], data[3][j], data[5][j]])
            if (data[4][j] >= 7.2) and (data[4][j] < 7.21):
                table_const15.append([data[1][j], data[3][j], data[5][j]])
    i = i + 1


# The pretty table of the data
table_out = []
table_out.append(["IMF exponent", "stellar EQW [A]", "", "nebular EQW [A]", "", "total=nebular+stellar [A]", ""])
table_out.append(["", "5x10^6 yr", "15x10^6 yr", "5x10^6 yr", "15x10^6 yr", "5x10^6 yr", "15x10^6 yr"])
table_out.append(["Instantanteous SFR", "", "", "", "", "", ""])
for i in range(0, len(IMFs)):
# for 5 million years
    t = [IMFs[i], repr(table_fixed5[i][0]), repr(table_fixed15[i][0]), repr(table_fixed5[i][1]), repr(table_fixed15[i][1]), 
         repr(table_fixed5[i][2]), repr(table_fixed15[i][2])]
    table_out.append(t)
table_out.append(["Constant SFR", "", "", "", "", "", ""])
for i in range(0, len(table_fixed5)):
    t = [IMFs[i], repr(table_const5[i][0]), repr(table_const15[i][0]), repr(table_const5[i][1]), repr(table_const15[i][1]),
         repr(table_const5[i][2]), repr(table_const15[i][2])]
    table_out.append(t)
f = open('../results/SB99eqws_NEW.txt', 'w+')
table.pprint_table(f, table_out)
f.close()

            
lower = 6.3
upper = 8
pyplot.figure(1, figsize=(10, 9))
pyplot.title('Instantaneous SFR')
pyplot.suptitle(lya_used)
pyplot.xlabel('log age [yr]')
pyplot.ylabel('Ly-alpha EQW [$\AA$]')
pyplot.xlim(lower, upper)
pyplot.ylim(-50, 225)
pyplot.text(6.33, -9, 'stellar component')   
pyplot.text(6.35, 100, 'total=')   
pyplot.text(6.36, 91, 'nebular')   
pyplot.text(6.41, 85, '+')  
pyplot.text(6.36, 78, 'stellar')
pyplot.text(6.6, 100, 'nebular component')  
pyplot.text(7.7, 209 , 'IMF exp')  
 
pyplot.figure(2, figsize=(10,9))
pyplot.title('Constant SFR') 
pyplot.suptitle(lya_used)
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
epsfile = os.path.join(path_results, "Inst_logAge_OLD.eps")
pyplot.savefig(epsfile)#, dpi=(100))

pyplot.figure(2)
epsfile = os.path.join(path_results, "Const_logAge_OLD.eps")
pyplot.ion()
pyplot.savefig(epsfile)
'''


