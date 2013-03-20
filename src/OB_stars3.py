import numpy
import string
import table
import copy
from spectrum import spectrum
from spectrum import pyplot

'''
This program produces a table of temperatures and EQWs of O and B stars from the textfiles with the selected wavelength range: 1100 - 1300 A 
from Pollux and Tlusty, respectively.

'''

def ec_VG93(Teff):
    VG_eqw = -420.0 * numpy.exp(-1*float(Teff)/6100.0) 
    return (VG_eqw)

def sort_list(the_list, first_sort_by_col=0, second_sort_by_col=1, total_columns=6):
    '''
    This function sorts the elements of a list by comparing two elements of the first column given, list[i] and list[i+1]=list[j].
    i. e.: list = [[1, 5, 2, 7, 2], [0.7, 1, 1.5, 2, 2.5]]
           ordered_list = list
            compare i=1 with j=5
                i<j so values_in_order=1
            compare i=5 with j=2
                i>j therefore temp=5=i and then i=2=j and j=temp=5
            compare i=5 with j=7
                i<j so leave the terms as they are so values_in_order=values_in_order+1=2
            compare i=7 with j=2
                i>j therefore temp=7=i and then i=2=j and j=temp=7
            * This gives first_column=[1, 2, 5, 2, 7], which is stil not in order but the highest value is already in its right place
              hence the ordering needs to happen one less time in the next iteration.
            * Repeat the iteration until values_in_order=len(list)-1
            When comparing i=2 with j=2 use the second column of the original list for these numbers, so i_original=1.5 and j_original=2.5
                i_original < j_original so values_in_order = values_in_order+1 = len(list)-1
            which gives the first column in order but not the second in the ordered list so now:
                - find the first element of first column of the ordered list (1) in the original list
                - make its corresponding second column value (0.7) the first element of the second column of the ordered list
                * repeat for the length of the column.
            this gives ordered_list = [[1, 2, 2, 5, 7], [0.7, 1.5, 2.5, 1, 2]]
    '''
    original_indeces = []
    for idx in range(0, len(the_list[0])):
        original_indeces.append(idx)
        
    ordered_list = copy.deepcopy(the_list)
    ordered_list_indeces = copy.deepcopy(original_indeces)
    temp = 0.0
    values_in_order = 0
    i = 0
    loop = 1
    N = len(original_indeces) 
    for j in range(1, N):
        #print('loop number %i' % loop)
        if ordered_list[first_sort_by_col][i] > ordered_list[first_sort_by_col][j]:
            print('i>j, rearranging')
            temp = ordered_list[first_sort_by_col][i]
            ordered_list[first_sort_by_col][i] = ordered_list[first_sort_by_col][j]
            ordered_list[first_sort_by_col][j] = temp
            # Do the same for to move the indeces
            temp_idx = ordered_list_indeces[i]
            ordered_list_indeces[i] = ordered_list_indeces[j]
            ordered_list_indeces[j] = temp_idx
            values_in_order = values_in_order + 1
            
        elif ordered_list[first_sort_by_col][i] < ordered_list[first_sort_by_col][j]:
            values_in_order = values_in_order + 1
            
        elif ordered_list[first_sort_by_col][i] == ordered_list[first_sort_by_col][j]:
            # right indeces to compare: get the index of the current element in column 1 of original list
            right_i = spectrum.find_index_in_list(original_indeces, ordered_list_indeces[i])
            right_j = spectrum.find_index_in_list(original_indeces, ordered_list_indeces[j])
            # now compare values of second column of original list
            if the_list[second_sort_by_col][right_i] > the_list[second_sort_by_col][right_j]:
                print('i>j, rearranging')
                temp = ordered_list[first_sort_by_col][i]
                ordered_list[first_sort_by_col][i] = ordered_list[first_sort_by_col][j]
                ordered_list[first_sort_by_col][j] = temp
                # Do the same for to move the indeces
                temp_idx = ordered_list_indeces[i]
                ordered_list_indeces[i] = ordered_list_indeces[j]
                ordered_list_indeces[j] = temp_idx
                                    
            elif ordered_list[second_sort_by_col][right_i] <= ordered_list[second_sort_by_col][right_j]:
                values_in_order = values_in_order + 1
        i = i + 1
        loop = loop+1
        print('values so far', values_in_order)
        print(ordered_list[first_sort_by_col])
    
    if values_in_order != N:
        i = 0
        values_in_order = 0
        order_by_one_less_place = N-1
        for j in range(1, order_by_one_less_place):
            if ordered_list[first_sort_by_col][i] < ordered_list[first_sort_by_col][j]:
                values_in_order = values_in_order + 1
                
            elif ordered_list[first_sort_by_col][i] > ordered_list[first_sort_by_col][j]:
                print('comparing', ordered_list[first_sort_by_col][i], ordered_list[first_sort_by_col][j])
                temp = ordered_list[first_sort_by_col][i]
                ordered_list[first_sort_by_col][i] = ordered_list[first_sort_by_col][j]
                ordered_list[first_sort_by_col][j] = temp
                # Do the same for to move the indeces
                temp_idx = ordered_list_indeces[i]
                ordered_list_indeces[i] = ordered_list_indeces[j]
                ordered_list_indeces[j] = temp_idx
                values_in_order = values_in_order + 1
                print('now', ordered_list[first_sort_by_col][i], ordered_list[first_sort_by_col][j])
                
            elif ordered_list[first_sort_by_col][i] == ordered_list[first_sort_by_col][j]:
                # right indeces to compare: get the index of the current element in column 1 of original list
                right_i = spectrum.find_index_in_list(original_indeces, ordered_list_indeces[i])
                right_j = spectrum.find_index_in_list(original_indeces, ordered_list_indeces[j])
                # now compare values of second column of original list
                if the_list[second_sort_by_col][right_i] > the_list[second_sort_by_col][right_j]:
                    #print('i>j, rearranging')
                    temp = ordered_list[first_sort_by_col][i]
                    ordered_list[first_sort_by_col][i] = ordered_list[first_sort_by_col][j]
                    ordered_list[first_sort_by_col][j] = temp
                    # Do the same for to move the indeces
                    temp_idx = ordered_list_indeces[i]
                    ordered_list_indeces[i] = ordered_list_indeces[j]
                    ordered_list_indeces[j] = temp_idx
                                        
                elif ordered_list[second_sort_by_col][right_i] <= ordered_list[second_sort_by_col][right_j]:
                    values_in_order = values_in_order + 1
            i = i + 1
            loop = loop+1
            order_by_one_less_place = order_by_one_less_place - 1
            #print('values so far', values_in_order)
            #print(ordered_list[first_sort_by_col])
    
    #print('first col ordered', ordered_list[first_sort_by_col])
    # The column1 of ordered_list is ordered so now arrange the all other columns with the same order as the first
    k = 0
    for idx in ordered_list_indeces:
        for col in range(0, total_columns):
            ordered_list[col][k] = the_list[col][idx]
        k = k + 1
            
    return(ordered_list)

testlist = [[1, 5, 2, 7, 2], [0.7, 1, 1.5, 2, 2.5], ['1', '2', '3', '4', '5']]
ordertestlist = sort_list(testlist, first_sort_by_col=0, second_sort_by_col=1, total_columns=3)
print(testlist)
print(ordertestlist)
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
#print(OBtemps_logg[0])
print(OBtemps_logg[1])
#print(OBtemps_logg[2])

# ordering by temp
sorted_by_temp = sort_list(OBtemps_logg, first_sort_by_col=1, second_sort_by_col=2, total_columns=7)
print(sorted_by_temp[1])
exit()

# Figure
fig, ax = pyplot.subplots()
axes = [ax, ax.twinx()]
bla = [ zip(sorted_by_temp[1], sorted_by_temp[2]), zip(sorted_by_temp[1], sorted_by_temp[3]), zip(sorted_by_temp[1], sorted_by_temp[4]), 
       zip(sorted_by_temp[1], sorted_by_temp[5]), zip(sorted_by_temp[1], sorted_by_temp[6])]
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

#pyplot.xlabel('Temperature [K]')
#pyplot.ylabel('log g')
#pyplot.plot(sorted_temps[1], sorted_temps[2], 'b')
pyplot.show()



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

print(OBtemps_logg[3])    
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



