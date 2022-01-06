# -*- coding: utf-8 -*-
"""
Created on Fri May 28 09:28:39 2021

@author: shan0679
"""

#This script calculates mass for the whole lungs

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import glob
from pathlib import Path
#import sys 
from scipy.signal import lfilter, filtfilt, butter

#from openpyxl import Workbook
#wb = Workbook()



rootdir = Path('H:\Ronan_projekter\LungeData\pt_1')
sub_dirs = [d for d in rootdir.glob('*') if d.is_dir()]
print(sub_dirs)


# all_files = []
# for d in sub_dirs:
#     for f in d.iterdir():
#         all_files.append(f)
    
MinIndex = 7
MaxIndex = MinIndex + 140 


All_voxels = {}
PD15_voxels = {}
masses = {}
for d in sub_dirs:
    Tot_mass = list()
    mass = list()
    PD15_voxel = list()
    PD15_density = list()
    PD15_index = list()
    PD15_maxIndex = list()
    Total_lobe_volume = list()
    Tot_voxels = list()
    Tot_volume = list()
    Running_sum = list()
    test_var = list()
    test_var_sum = list()
    #test_var_sum2 = list()
    for file in sorted(d.iterdir()):
        data = pd.read_csv (file)   #read the csv file (put 'r' before the path string to address any special characters in the path, such as '\'). Don't forget to put the file name at the end of the path + ".csv"
        test = pd.DataFrame(data, columns= ['PatientName','PatientID','PatientBirthDate','PatientSize'])
        BinMin = test.iloc[MinIndex:MaxIndex, 0]
        BinMax = test.iloc[MinIndex:MaxIndex, 1]
        Voxels = test.iloc[MinIndex:MaxIndex, 2]
        volumenml = test.iloc[MinIndex:MaxIndex, 3]
    
        V_ml = volumenml.values
        BinMin = BinMin.values
        BinMax = BinMax.values
        Voxels = Voxels.values
        
        V = [float(i) for i in V_ml]
        Y = [float(i) for i in BinMin]
        Z = [float(i) for i in BinMax]
        W = [float(i) for i in Voxels]
        
        Average_HU = np.divide(np.add(Z,Y),2) 
        Average_rho = np.add(np.divide(Average_HU,1000),1)
        
        
        #plt.figure()
        #plt.plot(Average_rho,W,'k', markersize=4)
        #plt.xlabel('Density (g/ml)')
        #plt.ylabel('Voxels ()')
        #plt.title(file.stem)
        
       
        Tot_voxels.append(W)
        Tot_volume.append(V)
        #print(Tot_voxels)
        #turning the list Tot_voxels into a numpy array so that I can add the elements from the different lobes
        Tot_voxels_arr = np.array(Tot_voxels)
        Tot_volume_arr = np.array(Tot_volume)
       
        #Tot_voxels_exp = sum(Tot_voxels_arr[:,0])
        
        
        #plt.figure()
        #plt.plot(Average_rho,V,'k', markersize=4)
        #plt.xlabel('Density (g/ml)')
        #plt.ylabel('Volumen (ml)')
        #plt.title(file.stem)
        
        mass.append(np.multiply(Average_rho,V))
        Tot_mass.append(sum(mass[-1]))
        
        
    
    print(file.stem) 
    #adding data from the 5 lobes to give the number of voxels in both lungs
    Lobe_sum = np.sum(Tot_volume_arr, axis = 0)
    #finding the total numer of voxels in both lungs: Lungs_voxels
    Lungs_volume = sum(Lobe_sum)
    print("Total volume", Lungs_volume) 
    #Finding the 15th percentile point
    PD15_volume = np.multiply(Lungs_volume,0.15)
    print("15 percent =", PD15_volume)
    
    
    #Calculating the cumulative sum
    Running_sum = np.cumsum(Lobe_sum)
    #print(Running_sum)
    # finding the 15% breakpoint
    break_point = Running_sum[-1]*0.15
    #print(break_point)
    #finding the index of the entry for which the value is greater than the breakpoint
    Max_arg = np.argmax(Running_sum >= break_point) + 1
    #print(Max_arg)
    #finding the HU value that corresponds to the 15th percentile
    print("15th percentile HU = ",Average_HU[Max_arg])
    print("15th percentile rho = ",Average_rho[Max_arg])
    
    #plt.figure()
    
    #plt.xlabel('Density (g/ml)')
    #plt.ylabel('Volume (ml)')
    #plt.title(file.stem)
   
    #plt.figure()
    #plt.plot(Average_HU,Lobe_sum,'k-', markersize=4)
    #plt.xlabel('HU')
    #plt.ylabel('Volume (ml)')
    #plt.title(file.stem)
    
    
    # Filtering the data using lfilter
    #n = 20
    #b = [1.0/n]*n
    #a = 1
    #yy = lfilter(b,a,Lobe_sum)

    #plt.figure()
    #plt.plot(Average_rho,yy,'b-', markersize = 4)
    #plt.xlabel('Density (g/ml)')
    #plt.ylabel('Volume (ml)')
    #plt.title('filtered data ' + file.stem)
  
    ## doing the data analysis on the lfiltered data
    #Lungs_volume_filtered = sum(yy)
    #print("L_Filtered Total volume = ", Lungs_volume_filtered)
    ##Finding the 15th percentile point from the filtered data
    #PD15_volume_filtered = np.multiply(Lungs_volume_filtered,0.15)
    #print("L_Filtered 15 percent =", PD15_volume_filtered)
    
    #Running_sum_filtered = np.cumsum(yy)
    #break_point_filtered = Running_sum_filtered[-1]*0.15
    #Max_arg_filtered = np.argmax(Running_sum_filtered >= break_point_filtered) + 1
    
    #print("L_Filtered 15th percentile HU = ",Average_HU[Max_arg_filtered])
    #print("L_Filtered 15th percentile rho = ",Average_rho[Max_arg_filtered])
    
    
    # Filtering the data using filtfilt
    b, a = butter(3, 0.08)
    yyy = filtfilt(b, a, Lobe_sum)
    
    plt.figure()
    plt.plot(Average_rho,yyy,'b-', markersize = 4)
    plt.plot(Average_rho,Lobe_sum,'k-', markersize=4)
    plt.xlabel('Density (g/ml)')
    plt.ylabel('Volume (ml)')
    plt.title('filt_filtered data ' + file.stem)
    
     # doing the data analysis on the FILTfiltered data
    Lungs_volume_filtfiltered = sum(yyy)
    print("FiltFiltered Total volume = ", Lungs_volume_filtfiltered)
    #Finding the 15th percentile point from the filtered data
    PD15_volume_filtfiltered = np.multiply(Lungs_volume_filtfiltered,0.15)
    print("Filt_Filtered 15 percent =", PD15_volume_filtfiltered)
    
    Running_sum_filtfiltered = np.cumsum(yyy)
    break_point_filtfiltered = Running_sum_filtfiltered[-1]*0.15
    Max_arg_filtfiltered = np.argmax(Running_sum_filtfiltered >= break_point_filtfiltered) + 1
    
    print("Filt_Filtered 15th percentile HU = ",Average_HU[Max_arg_filtfiltered])
    print("Filt_Filtered 15th percentile rho = ",Average_rho[Max_arg_filtfiltered])
    
#print("Running sum", R_sum)
 
#df = pd.DataFrame(masses)
#df2 = pd.DataFrame(Voxels_both_lungs)
#df.plot()
#print(df)
#print(df2)

# Calculate the total lung mass from the 3 CT scans: Exp, Insp and LD
#Tot_Mass_exp = np.sum(masses['massEXP'])
#Tot_Mass_insp = np.sum(masses['massINSP'])
#Tot_Mass_ld = np.sum(masses['massLD'])



