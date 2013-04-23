import os
import glob
import numpy
import table
import string
#matplotlib.use("macosx")
#from pprint import pprint
#from pylab import *
from random import randrange
import matplotlib
from matplotlib.ticker import MaxNLocator
from spectrum import spectrum 
from matplotlib import pyplot


'''
TESTING DIFFERENT values of SFR cor a constant case.

This code makes the plots of EQW vs log(age) using the output files of Starburst99.

In this version of the plots and simply add a zoom-in in the stellar part. 

'''

def read_S99file(txt_file):
    data = numpy.loadtxt(txt_file, dtype='float', skiprows=7, unpack=True)
    return data

# Path to the files
path_results = '../results/SB99_diffSFR/'
path_plots = '../results/plots/'
S99 = glob.glob('../Starburst99/different_SFR/*')

constant_ste = []
constant_neb = []
constant_tot = []
fixed_ste = []
fixed_neb = []
fixed_tot = []
i = 0
IMFs = ['0.5',
        '1.2',
        #'1.3',
        #'1.4',
        #'1.5',
        #'1.6',
        #'1.7',
        #'1.8',
        #'1.9', 
        #'2.0', 
        #'2.1',
        #'2.2',
        '2.3']#,
        #'2.5']

SFRs = ['0.05',
        #'0.25',
        #'0.5',
        #'0.75',
        '1.0',
        #'2.5',
        '5.0',
        '10.0']

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

# This loop only choses the selected SFRs to annalyze the files and create the plots
for each_folder in S99:
    for j in range (0, len(SFRs)):
        if SFRs[j] in each_folder:
            total_chosenfiles.append(each_folder)
print("%d chosen files" % len(total_chosenfiles))
    
print('Do you want to write the text files? (y for YES,  enter for NO)')
write_files = raw_input()
if write_files == 'y':
    for foldr in total_chosenfiles:
        print("File %i of %i" % (i+1, len(total_chosenfiles))) # for only a few SFRs
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
        #print('shape of time_and_lumLyA', shape(time_and_lumLyA))
        
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
        #print('shape of age_and_nebcontinuum', shape(age_and_nebcontinuum))
         
        
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
                #temp = [age_and_nebcontinuum[0][j], numpy.log10(age_and_nebcontinuum[0][j]), lLyA, age_and_nebcontinuum[1][j], eqw]
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

else:
    pass


files_in_results = glob.glob('../results/SB99_diffSFR/s99LyAeqw*_imf*.txt')
files2plot = len(files_in_results)
print('Total files to obtain information for the plot: %i' % files2plot)

# This loop makes sure that we are only picking the right files with the chosen SFRs
tot_lyas_files = []
for text_file in files_in_results:
    for sfr in SFRs:
        if sfr in text_file:
            tot_lyas_files.append(text_file)
            
'''
# This loop is to find the names of the selected IMFs, but in this case we only have 1 IMF
for text_file in files_in_results:
    file_name = string.split(text_file, sep='imf')
    exponent = file_name[1].replace(".txt", "")
    for imf in IMFs:
        if exponent == imf:
            if 'cont' in text_file:
                tot_lyas_files.append(text_file)
'''

# Files to compare to:
pivot23 = '../results/SB99_diffSFR/s99LyAeqwcont1.0_imf2.3.txt'
pivot_file23 = numpy.loadtxt(pivot23, dtype=float, skiprows=1, unpack=True)
diff2pivot23 = []
pivot12 = '../results/SB99_diffSFR/s99LyAeqwcont1.0_imf1.2.txt'
pivot_file12 = numpy.loadtxt(pivot12, dtype=float, skiprows=1, unpack=True)
diff2pivot12 = []

# lists to make the EQWs table
table_const5 = []
table_const10 = []
table_const15 = []
# This loop actually makes the plots
i = 1
for txt_file in tot_lyas_files:
    print('Working with %s' % txt_file)
    print('which is %i of %i' % (i, files2plot))
    data = numpy.loadtxt(txt_file, dtype=float, skiprows=1, unpack=True)
    #data = ste_age, ste_eqw, neb_age, neb_eqw, log_age, tot_eqw
    #        data[0], data[1], data[2], data[3], data[4], data[5]
    tf = string.split(txt_file, sep='imf')
    imf = tf[1].replace(".txt", "")
    f = path_results+'diff2pivot_'+'imf'+imf+'.txt'
    if imf == '1.2':
        #print >> f, 'Differences for IMF=%s to 1.2' % (imf)
        sdiff = pivot_file12[1] - data[1]
        ndiff = pivot_file12[3] - data[3]
        tdiff = pivot_file12[5] - data[5]
        diff2pivot12 = numpy.array([sdiff, ndiff, tdiff])
        output = numpy.column_stack((diff2pivot12))
        numpy.savetxt(f, output, delimiter=' ', fmt='%-10.3f')
    elif imf == '2.3':
        #print >> f, 'Differences for IMF=%s to 2.3' % (imf)
        sdiff = pivot_file23[1] - data[1]
        ndiff = pivot_file23[3] - data[3]
        tdiff = pivot_file23[5] - data[5]
        diff2pivot12 = numpy.array([sdiff, ndiff, tdiff])
        output = numpy.column_stack((diff2pivot12))
        numpy.savetxt(f, output, delimiter=' ', fmt='%-10.3f')
            
    # Divide continuous from instantaneous bursts
    #print(data)
    pair_ste = [data[4], data[1]]
    pair_neb = [data[4], data[3]]
    tot_age_and_eqw = [data[4], data[5]]
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

        
print('want to plot? (y for YES,  enter for NO)')
plt = raw_input()
if plt == 'y':
    pass
else:
    print('That was a nope...  Done!')
    exit()
    

'''
# The pretty table of the data
table_out = []
table_out.append(["IMF exponent", "stellar EQW [A]", "", "nebular EQW [A]", "", "total=nebular+stellar [A]", ""])
table_out.append(["", "5x10^6 yr", "15x10^6 yr", "5x10^6 yr", "15x10^6 yr", "5x10^6 yr", "15x10^6 yr"])
table_out.append(["Constant SFR", "", "", "", "", "", ""])
for i in range(0, len(table_fixed5)):
    t = [IMFs[i], repr(table_const5[i][0]), repr(table_const15[i][0]), repr(table_const5[i][1]), repr(table_const15[i][1]),
         repr(table_const5[i][2]), repr(table_const15[i][2])]
    table_out.append(t)
f = open('../results/SB99eqws_NEW.txt', 'w+')
table.pprint_table(f, table_out)
f.close()
'''    
    
# PLOTs of the SFRs vs log(age)            
lower = 6.3
upper = 8
font = {#'family' : 'Vera Sans',
        'weight' : 'regular',
        'size'   : 16}
# OR ONE CAN USE FOR INDIVIDUAL AXIS:
#f1.tick_params(axis='x', labelsize=17)
#f1.tick_params(axis='y', labelsize=17)
matplotlib.rc('font', **font)

figs = pyplot.figure(1, figsize=(8, 12))
figs.subplots_adjust(hspace=0.10)
#figs.suptitle()
f1 = pyplot.subplot(211)
f1.set_title('Constant SFR') 
# Set common labels
figs.text(0.5, 0.06, 'log age [yr]', ha='center', va='center')
figs.text(0.035, 0.5, 'Ly-alpha EW [$\AA$]', ha='center', va='center', rotation='vertical')
f1.set_xlim(lower, upper)
#f1.set_ylim(-15, 170)
f1.text(6.6, 6, 'stellar component')   
f1.text(7.15, 60, 'total = nebular + stellar')   
f1.text(6.8, 200, 'nebular component')   
f1.text(7.65, 275, 'SFR')  
# adding the enhancement of the stellar component  
f2 = figs.add_subplot(212)#, sharex=f3)
#f2.set_xlim(lower, upper) 
f2.set_xlim(6.5, 7.3)
#f2.set_xlim(6.5, 6.55)
f2.set_ylim(-8, -1.0)
#f2.set_ylim(-2.90, -2.85)
f2.text(6.67, -1.8, 'stellar component')   
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))

# Ploting
curves_cont =[]
line_styles = [':', '--', '.-', '-']
for i in range(0, 11):#len(constant_tot)):
    #pyplot.figure(2)
    curve, = f1.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve_c, = f1.plot(constant_neb[i][0], constant_neb[i][1], colors[i])
    curves_cont.append(curve_c)
    curve, = f1.plot(constant_tot[i][0], constant_tot[i][1], colors[i])
    curve.set_linestyle(':')
    # the stellar enhancement
    curve, = f2.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')
    print('done with file %i' % i)

leg1 = pyplot.figure(1).legend(curves_cont, IMFs, bbox_to_anchor=(0.47, 0.88), labelspacing=0.2)
leg2 = pyplot.figure(1).legend(curves_cont, SFRs, bbox_to_anchor=(0.87, 0.88), labelspacing=0.2)
for t in leg1.get_texts():
    t.set_fontsize(12)     
for t in leg2.get_texts():
    t.set_fontsize(12)     
pyplot.show()

'''
# Saving the plots
pyplot.figure(1)
epsfile = os.path.join(path_results, "Const_logAge_NEW2.eps")
pyplot.ion()
pyplot.savefig(epsfile)
'''

print('Done!   :P')
