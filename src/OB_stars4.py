import numpy
import string
import table
import operator
from spectrum import spectrum
from spectrum import pyplot

'''
This program produces a table of temperatures and EQWs of O and B stars from the textfiles with the selected wavelength range: 1100 - 1300 A 
from Pollux and Tlusty, respectively.

'''

def ec_VG93(Teff):
    VG_eqw = -420.0 * numpy.exp(-1*float(Teff)/6100.0) 
    return (VG_eqw)

def convert_to_tuple(ref):
    temp = []
    i = 0
    for i in range(len(ref[i])):
        temp.append(tuple([x[i] for x in ref]))
    return temp

def convert_to_list(ref):
    temp = []
    i = 0
    for i in range(len(ref[i])):
        temp.append(list([x[i] for x in ref]))
    return temp

testlist = [[1, 5, 2, 7, 2], [0.7, 1, 3.5, 2, 2.5], ['0', '1', '2', '3', '4']]
test_tuples = convert_to_tuple(testlist)
#sorted_test_tuples = sorted(test_tuples, key=lambda d: d[0])
#sorted_testlist = convert_to_list(sorted_test_tuples)
test_tuples_sorted = sorted(test_tuples, key=operator.itemgetter(0, 1))
sorted_testlist = convert_to_list(test_tuples_sorted)
print(testlist)
#print(sorted_test_tuples)
print(sorted_testlist)
#exit()

measured_temps = []
measured_loggs = []
alias = []
measured_eqws_pm5 = []
measured_eqws_pm10 = []
measured_eqws_pm20 = []
measured_eqws_pm30 = []

# Obtaining eqws, temperatures, and log g
path_OBstars = ['../results/OpolluxStars/', '../results/BTlustyStars/']
for each_folder in path_OBstars:
    print(each_folder)
    key_wd = 'Opollux'
    O_star = True
    if key_wd in each_folder:        
        txt = each_folder+'Ostar_eqws.txt'
    else:
        txt = each_folder+'eqws_COMPARISON2_updated.txt'
        O_star = False
        
    f = open(txt, 'r')
    # Ignore the first 2 lines
    hdr1 = f.readline() # sign convention
    hdr2 = f.readline() # columns titles line
    hdr3 = f.readline() # columns titles line
    for line in f:
        line = line.strip()
        columns = line.split()
        a = columns[0]
        alias.append(a)
        if O_star == True:
            eqw = float(columns[1])
            measured_eqws_pm5.append(eqw)
            measured_eqws_pm10.append(0.0)
            measured_eqws_pm20.append(0.0)
            measured_eqws_pm30.append(0.0)
        else:
            eqw1 = float(columns[1])
            measured_eqws_pm10.append(eqw1)
            eqw2 = float(columns[2])
            measured_eqws_pm20.append(eqw2)
            eqw3 = float(columns[3])
            measured_eqws_pm30.append(eqw3)

extra_lines = len(measured_eqws_pm20) - len(measured_eqws_pm5) 
for line in range(0, extra_lines):
    measured_eqws_pm5.append(0.0)
     
for i in range(0, len(alias)):
    name1 = string.split(alias[i], sep='_')
    if 'Ostar' in alias[i]:
        oname1 = string.split(name1[1], sep='g')
        Teff = float(oname1[0])
        logg = float(oname1[1])
        #print('Teff, logg', Teff, logg)    
        measured_temps.append(Teff)
        measured_loggs.append(logg)
    elif 'Norm_Bstar' in alias[i]:
        bname1 = string.split(name1[4], sep='g')
        Teff = float(bname1[0])
        bname2 = bname1[1].replace(".txt", "")
        logg = float(bname2)*0.01
        #print('Teff, logg', Teff, logg)    
        measured_temps.append(Teff)
        measured_loggs.append(logg)

OBtemps_logg = [alias, measured_temps, measured_loggs, measured_eqws_pm5, measured_eqws_pm10, measured_eqws_pm20, measured_eqws_pm30]
print(len(OBtemps_logg[0]), len(OBtemps_logg[1]), len(OBtemps_logg[2]), len(OBtemps_logg[3]), len(OBtemps_logg[4]), len(OBtemps_logg[5]), len(OBtemps_logg[6]))

# ordering by temp
sorted_by_temp = convert_to_tuple(OBtemps_logg)
sorted_test_tuples = sorted(sorted_by_temp, key=operator.itemgetter(1, 2))
sorted_OBtemps_logg = convert_to_list(sorted_test_tuples)
print(sorted_OBtemps_logg[1])
print(sorted_OBtemps_logg[2])
'''
# Printing file with desired eqws
table_data = []
#table_data.append("*** SIGN CONVENTION:   positive = emission   negative = absorption ***")
table_data.append(["Temperature", "log g", "LyA eqw +- 5A", "LyA eqw +- 20A", "Alias"])
N = len(sorted_OBtemps_logg[0])
for i in range(0, N):
    temp = [repr(sorted_OBtemps_logg[1][i]), "%.2f" % sorted_OBtemps_logg[2][i], "%.3f" % sorted_OBtemps_logg[3][i], 
            "%.3f" % sorted_OBtemps_logg[5][i], sorted_OBtemps_logg[0][i]]
    table_data.append(temp)
print(table_data[0])
f = open('../results/table_OB_desired_eqws.txt', 'w+')
table.pprint_table(f, table_data)
f.close()
exit()
'''

# Figure
fig, ax = pyplot.subplots()
axes = [ax, ax.twinx()]
bla = [ zip(sorted_OBtemps_logg[1], sorted_OBtemps_logg[2]), zip(sorted_OBtemps_logg[1], sorted_OBtemps_logg[3]), 
       zip(sorted_OBtemps_logg[1], sorted_OBtemps_logg[4]), zip(sorted_OBtemps_logg[1], sorted_OBtemps_logg[5]), 
       zip(sorted_OBtemps_logg[1], sorted_OBtemps_logg[6])]
color = ['b', 'g', 'r', 'magenta', 'k', 'cyan']
y_labels = ['log g', 'Ly-alpha EQW']
i = 0
for a, d in zip(axes, bla):
    x, y = zip(*d)
    #print("x", x)
    #print("y", y)
    a.plot(x, y, color=color[i])
    a.set_ylabel(y_labels[i])
    i = i+1 
axes[0].set_xlabel('Temperature [K]')
pyplot.show()
#exit()


# Performing the interpolation for the table 
temperatures = [15.0, 18.0, 21.0, 24.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 50.0]
loggs = [1.75, 2.0, 2.25, 2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5]
VG_eqws = []
table_data = []
table_eqws = []
#eqw0 = 
#table_eqws.append(eqw0)

# Smoothing eqw using Valls-Gabaud (1993)
for T in measured_temps:
    VGeqw = ec_VG93(T)
    #print('VGeqw', VGeqw)
    VG_eqws.append(VGeqw)

for t in temperatures:
    teff = t * 1000.0
    teff_table = spectrum.find_nearest(OBtemps_logg[1], teff)
    for g in loggs:
        logg_table = spectrum.findXinY(OBtemps_logg[2], OBtemps_logg[1], g)
        #eqw_tab = 
        line_table = [repr(teff_table), repr(logg_table), repr(eqw_table)]
        table_data.append(line_table)

        print(table_data)

table_out = open('../results/table_OBlyas.txt', 'w+')
table.pprint_table(table_out, table_data)



