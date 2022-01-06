#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 13:10:16 2021

@author: sofie
#This script  calculates the mass of the lungs for each lobe"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import glob
from pathlib import Path
#import sys 
#from scipy import integrate
#from numpy import trapz
from scipy.signal import filtfilt, butter

rootdir = Path('H:\Ronan_projekter\LungeData\pt_1')
sub_dirs = [d for d in rootdir.glob('*') if d.is_dir()]
print(sub_dirs)


# all_files = []
# for d in sub_dirs:
#     for f in d.iterdir():
#         all_files.append(f)
    
MinIndex = 7
MaxIndex = MinIndex + 140

#All_voxels = {}
Volume = {}
masses = {}
for d in sub_dirs:
    Tot_mass = list()
    mass = list()
    PD15_voxel = list()
    PD15_density = list()
    PD15_index = list()
    PD15_maxIndex = list()
    Tot_lobe_volume = list()
    Tot_voxels = list()
    Running_sum = list()
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
        
        #print(V)
        
        Average_HU = np.divide(np.add(Z,Y),2) 
        Average_rho = np.add(np.divide(Average_HU,1000),1)
        
        #plt.figure()
        #plt.plot(Average_rho,W,'k', markersize=4)
        #plt.xlabel('Density (g/ml)')
        #plt.ylabel('Voxels ()')
        #plt.title(file.stem)
        
        # Filtering the data using filtfilt
        b, a = butter(3, 0.08)
        yy = filtfilt(b, a, V)
        
        # calculating the total numbe rof voxels in the lobe
        Tot_voxels.append(sum(W))
        #Tot_lobe_volume.append(sum(V))
        Tot_lobe_volume.append(sum(yy))
        
        
        
        #mass.append(np.multiply(Average_rho,V))
        mass.append(np.multiply(Average_rho,yy))
        Tot_mass.append(sum(mass[-1]))
        print(file.stem)
    
        
        
        
        plt.figure()
        plt.plot(Average_rho,V,'k', markersize=4)
        plt.plot(Average_rho,yy,'r-', markersize = 4)
        plt.xlabel('Density (g/ml)')
        plt.ylabel('Volume (ml)')
        plt.title('filt_filtered data ' + file.stem)
    
    Volume["Lobe volume " + d.name] = Tot_lobe_volume
    #All_voxels["All voxels" + d.name]= Tot_voxels
    #PD15_voxels["PD15 voxels" + d.name] = np.multiply(Tot_voxels,0.15)
    masses["mass " + d.name] = Tot_mass
    #pd15s["PD15s" + d.name] = PD15
    #print("PD15 = ",PD15)

    

#print("Running sum", R_sum)
 
df = pd.DataFrame(masses)
df2 = pd.DataFrame(Volume)
df.plot()
df2.plot()
print(df)
print(df2)

# Calculate the total lung mass from the 3 CT scans: Exp, Insp and LD
Tot_Mass_exp = np.sum(masses['mass EXP'])
Tot_Mass_insp = np.sum(masses['mass INSP'])
Tot_Mass_ld = np.sum(masses['mass LD'])



print('Tot mass exp =',Tot_Mass_exp)
print('Tot mass insp =',Tot_Mass_insp)
print('Tot mass LD =',Tot_Mass_ld)


# Calculate the total lung volume from the 3 CT scans: Exp, Insp and LD
Tot_Vol_exp = np.sum(Volume['Lobe volume EXP'])
Tot_Vol_insp = np.sum(Volume['Lobe volume INSP'])
Tot_Vol_ld = np.sum(Volume['Lobe volume LD'])

print('Tot volume exp =',Tot_Vol_exp)
print('Tot volume insp =',Tot_Vol_insp)
print('Tot volume LD =',Tot_Vol_ld)


# finding 15% of the voxels
#PD15_exp = np.multiply(All_voxels['All voxelsEXP'],0.15)
#PD15_insp = np.multiply(All_voxels['All voxelsINSP'],0.15)
#PD15_ld = np.multiply(All_voxels['All voxelsLD'],0.15)

#finding the x-value that separates the 15% from the remaining 85%
