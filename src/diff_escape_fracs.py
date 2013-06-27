import os
import glob
import numpy
#import table
import matplotlib
matplotlib.use("macosx")
#from pprint import pprint
#from pylab import *
#from random import randrange
from matplotlib.ticker import MaxNLocator
from spectrum import spectrum 
from matplotlib import pyplot

'''
This code makes the plots of EQW vs log(age) using the output files of Starburst99, and considering different escape fractions of the
nebular lyman alpha. The plots have a zoom-in in the stellar part. 

'''

def read_S99file(txt_file):
    data = numpy.loadtxt(txt_file, dtype='float', skiprows=7, unpack=True)
    return data

# Path to the files
path_results = '../results/SB99_round4/'
path_plots = '../results/plots/'
S99 = glob.glob('../Starburst99/round4/*')

constant_ste = []
constant_neb = []
constant_tot = []
fixed_ste = []
fixed_neb = []
fixed_tot = []
i = 0
IMFs = [#'0.5',
        #'0.8',
        #'1.1',
        #'1.2',
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
        #'2.5',
        #'2.6']

escape_fraction = [0.05, 0.1, 0.25, 0.5, 1.0]

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

# Once with the chosen files go over them
for foldr in total_chosenfiles:
    print("%d chosen files" % len(total_chosenfiles))
    print("File %i of %i" % (i+1, len(total_chosenfiles))) # for only a few IMFs
    print foldr
    # No need to increase i because it is being increased in the next line
    file_name = os.path.basename(foldr)
    
    # THIS IS FOR THE NEBULAR COMPONENT
    '''
    # list of the arrays of age and lyman-alpha luminosity for different escape fractions
    time_and_lumLyA_escfrac = [] 
    
    # Get luminosity of Lyman-alpha
    logLHalpha_file = '/Output/'+file_name+'.ewidth'
    txt_file = foldr+logLHalpha_file
    data = read_S99file(txt_file)
    for ef in escape_fraction:
        lumLyA = ((10**data[2]) * 8.9) * ef
        time_and_lumLyA = numpy.array([data[0], lumLyA])
        time_and_lumLyA_escfrac.append(time_and_lumLyA)
    '''
    # Get luminosity of Lyman-alpha
    logLHalpha_file = '/Output/'+file_name+'.ewidth'
    txt_file = foldr+logLHalpha_file
    data = read_S99file(txt_file)
    lumLyA = (10**data[2]) * 8.9
    time_and_lumLyA = numpy.array([data[0], lumLyA])
        
    # Nebular continuum LyA
    logneb_contLyA_file = '/Output/'+file_name+'.spectrum'
    txt_file = foldr+logneb_contLyA_file
    data = read_S99file(txt_file)
    indexes=numpy.where( data == 1195.00 )      # I chose 1195 because lyman-alpha is wide and I wanted to get it all
    age = []
    nebcontinuum =[]
    # list of the nebular continuum for the different escape fractions
    list_age_and_nebcontinuum = []
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
    pairs_neb_age_and_eqw = []
    # get the specific escape fractions of the luminosity for the selected ages
    for ef in escape_fraction:
        neb_age = []
        neb_eqw = []
        for j in range(0, len(age_and_nebcontinuum[0])):
            if age_and_nebcontinuum[0][j] in time_and_lumLyA[0]:
                lLyA = spectrum.findXinY(time_and_lumLyA[1], time_and_lumLyA[0], age_and_nebcontinuum[0][j]) * ef
                eqw = lLyA / age_and_nebcontinuum[1][j]  # dealing with anti-logs
                #eqw = lLyA - age_and_nebcont[1][j]   # dealing with logs
                neb_age.append(age_and_nebcontinuum[0][j])
                neb_eqw.append(eqw)
                temp = [repr(age_and_nebcontinuum[0][j]), repr(numpy.log10(age_and_nebcontinuum[0][j])), repr(lLyA), repr(age_and_nebcontinuum[1][j]), repr(eqw)]
                table_age_and_nebcontinuum.append(temp)
        neb_age_and_eqw = numpy.array([neb_age, neb_eqw])
        pairs_neb_age_and_eqw.append(neb_age_and_eqw)
        '''
        # Saving the table into file 
        table_out = open(path_results+'age_and_nebcontinuum_'+file_name+'.txt', 'w+')    
        table.pprint_table(table_out, table_age_and_nebcontinuum)
        table_out.close()
        print('file  %s   is now written!' % (path_results+'age_and_nebContm_'+file_name+'_EscFrac'+repr(ef)+.txt'))
        '''

    '''
    # Finding corresponding LyA luminosities and times with their nebular parts
    table_age_and_nebcontinuum = []
    hdr = ['age_nebcontinuum', 'log age', 'LyA', 'nebcontinuum', 'neb_eqw']
    table_age_and_nebcontinuum.append(hdr)
    pairs_neb_age_and_eqw = []
    for time_and_lumLyA in time_and_lumLyA_escfrac:
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
        pairs_neb_age_and_eqw.append(neb_age_and_eqw)
    '''
    
    # Stellar eqw(LyA)
    ste_eqwLyA_file = '/Output/'+file_name+'.irfeature'
    txt_file = foldr+ste_eqwLyA_file
    data = read_S99file(txt_file)
    age = []
    eqwlya = []
    # use the first neb_age_and_eqw to find the corresponding stellar part with the same age
    for j in range(0, len(data[0])):
        if data[0][j] in pairs_neb_age_and_eqw[0][0]:
            #age.append(data[0][j]/1e6)            # age in mega years    
            #age.append(numpy.log10(data[0][j]))   # age in log
            age.append(data[0][j])                 # age in scientific notation
            #eqwlya.append(data[7][j]*(-1))        # the negative sign is because EQWs are in absorption (ONLY used for initial lyas table)
            eqwlya.append(data[7][j])       # no need to use the negative sign (EQWs are in absorption), data is fixed
        
    # length check
    # use the first neb_age_and_eqw to check lengths are the same
    while len(eqwlya) > len(pairs_neb_age_and_eqw[0][1]):
        print ('Fixing dimensions of arrays so that neb_age_and_eqw = ste_age_and_eqw')
        eqwlya.pop(len(eqwlya)-1)
        age.pop(len(age)-1)
    ste_age_and_eqw = numpy.array([age, eqwlya])
    
    ni = 0
    # Total EQW = nebular + stellar
    pairs_tot_age_and_eqw = []
    # use the neb_age_and_eqw in pairs_neb_age_and_eqw of the different escape fractions
    for neb_age_and_eqw in pairs_neb_age_and_eqw:
        tot_eqw = []
        tot_age = []
        for j in range(0, len(ste_age_and_eqw[1])):
            tot = ste_age_and_eqw[1][j] + neb_age_and_eqw[1][j]
            tot_eqw.append(tot)
            tot_age.append(numpy.log10(ste_age_and_eqw[0][j]))
        tot_age_and_eqw = numpy.array([tot_age, tot_eqw])
        pairs_tot_age_and_eqw.append(tot_age_and_eqw)
        
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

        #print('fixed', len(fixed_ste), len(fixed_tot), len(fixed_neb))
        #print('constant', len(constant_ste), len(constant_tot), len(constant_neb))
        #print("Inner for loop (fill arrays) iteration: {}".format(ni))
        ni += 1
    
        '''
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
        '''
    #print("Outer for loop iteration: {}".format(i))
    i = i+1

print(len(fixed_ste), len(fixed_tot), len(fixed_neb))
print(len(constant_ste), len(constant_tot), len(constant_neb))


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

figs = pyplot.figure(1, figsize=(10, 8))
#figs.subplots_adjust(hspace=0.10)
#figs.suptitle()
# Set common labels
figs.text(0.5, 0.06, 'log age [yr]', ha='center', va='center')
figs.text(0.035, 0.5, 'Ly-alpha EW [$\AA$]', ha='center', va='center', rotation='vertical')

f1 = figs.add_subplot(111)
f1.set_title('Instantaneous SFR')
#[label.set_visible(False) for label in f1.get_xticklabels()]
#step(f1.get_xticklabels(), visible=False)    # this is the same as the previous line but for some reason it is not working
f1.set_xlim(lower, upper)
f1.set_ylim(-50, 200)
f1.text(6.35, -17, 'stellar component')   
f1.text(6.33, 114, 'total=')   
f1.text(6.48, 114, 'nebular + stellar')   
f1.text(7.2, 5, 'nebular component')  
f1.text(7.5, 173, 'Escape Fraction')
'''
# adding the enhancement of the stellar component  
f2 = figs.add_subplot(212)
f2.set_xlim(6.5, 7.3)
f2.set_ylim(-20, -2.)
f2.text(6.6, -14, 'stellar component')
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
'''

figs = pyplot.figure(2, figsize=(10, 8)) 
#figs.subplots_adjust(hspace=0.10)
f3 = pyplot.subplot(111)
f3.set_title('Constant SFR') 
# Set common labels
figs.text(0.5, 0.06, 'log age [yr]', ha='center', va='center')
figs.text(0.035, 0.5, 'Ly-alpha EW [$\AA$]', ha='center', va='center', rotation='vertical')
f3.set_xlim(lower, upper)
f3.set_ylim(-30, 250)
f3.text(6.4, -15, 'stellar component')   
f3.text(7.35, 67, 'total = nebular + stellar')   
f3.text(6.67, 150, 'nebular component')   
f3.text(7.5, 220, 'Escape Fraction')  
'''
# adding the enhancement of the stellar component  
f4 = figs.add_subplot(212)#, sharex=f3)
f4.set_xlim(6.5, 7.3)
f4.set_ylim(-10, -0.5)
f4.text(6.6, -6.4, 'stellar component')   
# remove the first tick so that they do not overlap
pyplot.gca().yaxis.set_major_locator(MaxNLocator(prune='lower'))
'''

# Ploting
curves_inst = []
curves_cont =[]
line_styles = [':', '--', '.-', '-']


for i in range(0, len(constant_tot)):
    
    curve, = f1.plot(fixed_ste[i][0], fixed_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve_i, = f1.plot(fixed_neb[i][0], fixed_neb[i][1], colors[i])
    curves_inst.append(curve_i)
    curve, = f1.plot(fixed_tot[i][0], fixed_tot[i][1], colors[i])
    curve.set_linestyle(':')
    # the stellar enhancement
    #curve, = f2.plot(fixed_ste[i][0], fixed_ste[i][1], colors[i])
    #curve.set_linestyle('--')
    
    #pyplot.figure(2)
    curve, = f3.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    curve.set_linestyle('--')
    curve_c, = f3.plot(constant_neb[i][0], constant_neb[i][1], colors[i])
    curves_cont.append(curve_c)
    curve, = f3.plot(constant_tot[i][0], constant_tot[i][1], colors[i])
    curve.set_linestyle(':')
    # the stellar enhancement
    #curve, = f4.plot(constant_ste[i][0], constant_ste[i][1], colors[i])
    #curve.set_linestyle('--')

leg1 = pyplot.figure(1).legend(curves_inst, escape_fraction, bbox_to_anchor=(0.83, 0.82))  
leg2 = pyplot.figure(2).legend(curves_cont, escape_fraction, bbox_to_anchor=(0.83, 0.82))#, labelspacing=0.2)
for t in leg1.get_texts():
    t.set_fontsize(14)     
for t in leg2.get_texts():
    t.set_fontsize(14)     
pyplot.show()

'''
# Saving the plots
pyplot.figure(1)
epsfile = os.path.join(path_plots, "Inst_diffEscFrac.eps")
pyplot.savefig(epsfile)#, dpi=(100))

pyplot.figure(2)
epsfile = os.path.join(path_plots, "Const_diffEscFrac.eps")
pyplot.ion()
pyplot.savefig(epsfile)
'''

print('Done!')


