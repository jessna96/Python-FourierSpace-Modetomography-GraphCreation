#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 13:31:49 2022

@author: Jessica
"""

#%% ------ libraries, general styling and folder paths -----
# -- libraries --
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# -- general plot styles --
matplotlib.rcParams['font.monospace'] = "Arial"
matplotlib.rcParams['font.family'] = "monospace"
matplotlib.rc('font', size=16)
matplotlib.rcParams['axes.linewidth'] = 1

# -- data folder path --
data_path = "<path of sample data>" + "/" # enter path of sample data
filename_fourierspace = "fourier_space_data.txt" # this is the name of the sample data file
filename_realspace = "real_space_data.txt" # this is the name of the sample data file

filename = "Hyperspectral_and_Tomo_above_threshold"

# -- save option -- 
saveFigures = False

#%% functions
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
    
def get_brightest_point(data):
    brightestPoint = np.unravel_index(data.argmax(), data.shape) # y, x
    value = data[brightestPoint[0],brightestPoint[1]] 
    return brightestPoint, value

def get_intensity_k_energy_of_fourierspace_txt(path_of_txt_dispersion):
    with open(path_of_txt_dispersion) as f:
        lines = f.readlines()
        kTemp = [float(e) if isfloat(e) else e for e in lines[0].split(',')]
        kValues = kTemp[1:]
        y1 = [[float(e) if isfloat(e) else e for e in line.split(',')] for line in lines]
        y2 = y1[1:]
        energyValues = [row[0] for row in y2]
        intensityValuesTemp = [row[1:] for row in y2]
        intensityValuesTemp.reverse()
        intensityValues = np.array(intensityValuesTemp)
        
        extent = [kValues[0],kValues[len(kValues)-1],energyValues[0],energyValues[len(energyValues)-1]]
        
        data = {'IntensityValues': intensityValues, 'KValues': kValues, 'EnergyValues': energyValues, 'Extent': extent}
        
        return data
    
def get_values_of_isoEnergyCut_txt(path_of_txt_isoEnergyCut):
    with open(path_of_txt_isoEnergyCut) as f:
        lines = f.readlines()
        xTemp = [float(e) if isfloat(e) else e for e in lines[0].split(',')]
        xValues = xTemp[1:]
        y1 = [[float(e) if isfloat(e) else e for e in line.split(',')] for line in lines]
        y2 = y1[1:]
        yValues = [row[0] for row in y2]
        intensityValuesTemp = [row[1:] for row in y2]
        intensityValuesTemp.reverse()
        intensityValues = np.array(intensityValuesTemp)
        
        extent = [xValues[0],xValues[len(xValues)-1],yValues[0],yValues[len(yValues)-1]]
        
        data = {'IntensityValues': intensityValues, 'xkValues': xValues, 'ykValues': yValues, 'Extent': extent}
        
        return data
    
def slice_data(dispersionIntensityData, kData, energyData, lowerELimit, upperELimit, lowerKLimit, upperKLimit):
    for k in range(len(energyData)):
        if energyData[k] > lowerELimit:
            indexLowerE = dispersionIntensityData.shape[0] - k + 1
            break
              
    for l in range(len(energyData)):
        if energyData[l] > upperELimit:
            indexUpperE = dispersionIntensityData.shape[0] - l - 1 
            break
        
    for m in range(len(kData)):
        if kData[m] > lowerKLimit:
            indexLowerk = m - 1
            break
              
    for n in range(len(kData)):
        if kData[n] > upperKLimit:
            indexUpperk = n + 1
            break
        
    intensityValuesSliced = dispersionIntensityData[indexUpperE:indexLowerE, indexLowerk:indexUpperk]
    intensityValuesSliced_normalized = intensityValuesSliced/get_brightest_point(intensityValuesSliced)[1]
    
    extent=[kData[indexLowerk],kData[indexUpperk],energyData[dispersionIntensityData.shape[0]-indexLowerE],energyData[dispersionIntensityData.shape[0]-indexUpperE]]
    
    data = {'Data': intensityValuesSliced_normalized, 'Extent': extent}
    
    return data

    
def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    

#%% ------ create custom colormaps from colorlist ------
colorlist_fourierspace = [(1,1,1), (62/255, 14/255, 56/255), (103/255, 24/255, 101/255), (248/255, 101/255, 46/255), (255/255, 255/255, 220/255)]  # R -> G -> B
n_bins_fourierspace = 200  # discretizes the interpolation into bins
cmap_name_fourierspace = 'cmap_fourierspace'
colormap_fourierspace = LinearSegmentedColormap.from_list(cmap_name_fourierspace, colorlist_fourierspace, N=n_bins_fourierspace)

colorlist_realspace = [(255/255,255/255,255/255),(43/255, 0/255, 0/255),(155/255,0/255,0/255),(255/255,35/255,0/255),(255/255,139/255,0/255),(255/255,255/255,34/255),(255/255,255/255,219/255)]
n_bins_realspace = 500
cmap_name_realspace = 'cmap_realspace'
colormap_realspace = LinearSegmentedColormap.from_list(cmap_name_realspace, colorlist_realspace, N=n_bins_realspace)


#%% ------ limit the energy range and realspace ranges ------
selectedlowerEnergyLimit = 1.606
selectedupperEnergyLimit = 1.614
selectedlowerkLimit = -3.0
selectedupperkLimit = 3.0
# ----------------------------------------
selectedYmax = 30
selectedYmin = 10
selectedXmax = 30
selectedXmin = 10

#%% ------ plot data ------
# -- format fourier space --
fourierspace_data = get_intensity_k_energy_of_fourierspace_txt(data_path + filename_fourierspace)
sliced_fourierspace_data = slice_data(fourierspace_data['IntensityValues'], fourierspace_data['KValues'], fourierspace_data['EnergyValues'],selectedlowerEnergyLimit, selectedupperEnergyLimit, selectedlowerkLimit, selectedupperkLimit)

# -- format real space --
realspace_data = get_values_of_isoEnergyCut_txt(data_path + filename_realspace)
sliced_frealspace_data = slice_data(realspace_data['IntensityValues'], realspace_data['xkValues'], realspace_data['ykValues'], selectedYmin, selectedYmax, selectedXmin, selectedXmax)

# -- set aspect ratios --
kDiff = selectedupperkLimit-selectedlowerkLimit
energyDiff = selectedupperEnergyLimit-selectedlowerEnergyLimit
aspect_fourierspace = kDiff/energyDiff 
aspect_tomo = 1
    
# -- fourier space plot --
plt.imshow(sliced_fourierspace_data['Data'], cmap=colormap_fourierspace, aspect=aspect_fourierspace, extent=sliced_fourierspace_data['Extent'])
ax = plt.gca();
forceAspect(ax,aspect=1.0)
ax.tick_params(which="major", direction='in')
ax.tick_params(which="minor", direction='in')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xlabel('$k_\mathregular{y}$ ($\mu \mathregular{m}^{-1}$)')
plt.ylabel('Energy (eV)')
v1 = np.linspace(0, 1, 6, endpoint=True)
cbar = plt.colorbar(pad = 0.02, ticks=v1)
cbar.set_label('Intensity (arb. u.)', fontsize=12, rotation=90, labelpad=-38, color="white")
if saveFigures:    
    plt.savefig(data_path + 'fourier_space_image' + '.png', format='png', bbox_inches='tight', dpi=1000)
plt.show()

# -- real space plot --
plt.imshow(sliced_frealspace_data['Data'], cmap=colormap_realspace, aspect=aspect_tomo, extent=[0,20,0,20])
ax = plt.gca();
forceAspect(ax,aspect=1.0)
ax.tick_params(which="major", direction='in')
ax.tick_params(which="minor", direction='in')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xlabel('y ($\mu \mathregular{m}$)')
plt.ylabel('x ($\mu \mathregular{m}$)')
v1 = np.linspace(0, 1, 6, endpoint=True)
cbar = plt.colorbar(pad = 0.02, ticks=v1)
cbar.set_label('Intensity (arb. u.)', fontsize=12, rotation=90, labelpad=-38, color="white")
if saveFigures:    
    plt.savefig(data_path + 'real_space_image' + '.png', format='png', bbox_inches='tight', dpi=3000)
plt.show()

