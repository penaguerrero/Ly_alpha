import os
import glob
import numpy
import table
import matplotlib
matplotlib.use("macosx")
from spectrum import pyplot
from spectrum import spectrum
from random import randrange


def read_S99file(txt_file):
    data = numpy.loadtxt(txt_file, dtype='float', skiprows=7, unpack=True)
    return data

# Path to the files
path_results = '../results/plots'
S99 = glob.glob('../Starburst99/round2/*')

constant_ste = []
constant_neb = []
constant_tot = []
fixed_ste = []
fixed_neb = []
fixed_tot = []
i = 0
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
#colors = ['blue', 'red', 'green', 'magenta', '#994A2B', '#2B494A', '#7E1BE0', '#1BE0D3', '#E0531B', '#94E01B', '#9868C4']
file_count = 11
colors = []

def generate_colors(count):
    '''return value is list containing count number of colors '''
    colors = []
    for _ in range(file_count):
        colors.append("#%s" % "".join([hex(randrange(1, 254))[2:] for _ in range(3)]))
    return colors

colors = generate_colors(file_count)

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
    table_out = open('../results/age_and_nebcontinuum_'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_age_and_nebcontinuum)
    table_out.close()
    
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
            eqwlya.append(data[7][j]*(-1))         # the negative sign is because these EQWs are in absorption
        
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
    table_out = open('../results/time_and_lumLyA_'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_time_and_lumLyA)
    table_out.close()
        
    table_s99LyAeqw = [] 
    hdr = ['ste_age', 'ste_eqw', 'neb_age', 'neb_eqw', 'log age', 'tot_eqw']
    table_s99LyAeqw.append(hdr)
    for j in range(0, len(pair_ste[0])):
        temp = [repr(ste_age_and_eqw[0][j]), repr(ste_age_and_eqw[1][j]), repr(neb_age_and_eqw[0][j]), repr(neb_age_and_eqw[1][j]), repr(tot_age_and_eqw[0][j]), repr(tot_age_and_eqw[1][j])]
        table_s99LyAeqw.append(temp)
    #table.pprint_table(sys.stdout, table_s99LyAeqw)
    # Saving the table into file 
    table_out = open('../results/s99LyAeqw'+file_name+'.txt', 'w+')    
    table.pprint_table(table_out, table_s99LyAeqw)
    table_out.close()
        
    i = i+1
    
lower = 6.3
upper = 8
pyplot.figure(1, figsize=(10, 10))
pyplot.title('Instantaneous SFR')
pyplot.xlabel('log age [yr]')
pyplot.ylabel('Ly-alpha EQW [$\AA$]')
pyplot.xlim(lower, upper)
#pyplot.ylim(-15, 225)
pyplot.text(6.32, -15, 'stellar component')   
pyplot.text(6.33, 100, 'total=')   
pyplot.text(6.33, 89, 'nebular')   
pyplot.text(6.36, 78, '+')  
pyplot.text(6.33, 66, 'stellar')
pyplot.text(6.6, 100, 'nebular component')  
pyplot.text(8.01, 186 , 'IMF exp')  
 
pyplot.figure(2, figsize=(10,10))
pyplot.title('Constant SFR') 
pyplot.xlabel('log age [yr]')
pyplot.ylabel('Ly-alpha EQW [$\AA$]')
pyplot.xlim(lower, upper)
#pyplot.ylim(-15, 170)
pyplot.text(6.4, 6, 'stellar component')   
pyplot.text(7.3, 69, 'total=nebular+stellar')   
pyplot.text(7.3, 168, 'nebular component')   
pyplot.text(8.01, 236, 'IMF exp')   

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

leg1 = pyplot.figure(1).legend(curves_inst, IMFs, loc=7)#, bbox_to_anchor=(1.05, 1), borderaxespad=0.)    #this places the box outside
leg2 = pyplot.figure(2).legend(curves_cont, IMFs, loc=7)#, bbox_to_anchor=(1.05, 1), borderaxespad=0.) 
for t in leg1.get_texts():
    t.set_fontsize('small')     
for t in leg2.get_texts():
    t.set_fontsize('small')     
pyplot.show()
'''
# Saving the plots
pyplot.figure(1)
epsfile = os.path.join(path_results, "Inst_log2.svg")
pyplot.savefig(epsfile, dpi=(100))
epsfile = os.path.join(path_results, "Inst_log2.pdf")
pyplot.savefig(epsfile, dpi=(100))
epsfile = os.path.join(path_results, "Inst_log2.png")
pyplot.savefig(epsfile, dpi=(100))
epsfile = os.path.join(path_results, "Inst_log2.ps")
pyplot.savefig(epsfile, dpi=(100))


pyplot.figure(2)
epsfile = os.path.join(path_results, "Const_log2.eps")
pyplot.ion()
pyplot.savefig(epsfile)
'''
