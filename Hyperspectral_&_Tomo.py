#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 13:31:49 2022

@author: Jessica
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import sif_reader
from palettable.cmocean.sequential import Tempo_9
from palettable.cubehelix import cubehelix1_16
from copy import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import ticker

matplotlib.rcParams['font.monospace'] = "Arial"
matplotlib.rcParams['font.family'] = "monospace"
matplotlib.rc('font', size=16)
matplotlib.rcParams['axes.linewidth'] = 1

fontsize = 16

saveFigures = True

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

def get_intensity_k_energy_of_dispersion_txt(path_of_txt_dispersion):
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
    
def slice_dispersionData(dispersionIntensityData, kData, energyData, lowerELimit, upperELimit, lowerKLimit, upperKLimit):
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

def show_dispersion_data(data, extent):
    plt.imshow(data, cmap='turbo', aspect='auto', extent=extent)
    plt.colorbar()
    plt.show()
    
def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    
def get_char(i, alphabet = "abcdefghijklmnopqrstuvwxyz"):
    return(alphabet[i % len(alphabet)])

#%% data 
path_disperionAbovePth = "/Users/Jessica/Library/CloudStorage/OneDrive-UniversitätWürzburg/Meine Dateien/0_Messdaten/2022-03-22_C4892_Kagome_Tomos_Spotsize_Barrieren/A41_b3_f40_20erObj_Tomo_ueberSchwelle/Hyperspectral/A41_b3_1p9mW_1.txt"
# path_tomoAbovePth = "/Users/Jessica/Library/CloudStorage/OneDrive-UniversitätWürzburg/Meine Dateien/0_Messdaten/2022-03-22_C4892_Kagome_Tomos_Spotsize_Barrieren/A41_b3_f40_20erObj_Tomo_ueberSchwelle/energy1p60999.txt"
path_tomoAbovePth = "/Users/Jessica/Library/CloudStorage/OneDrive-UniversitätWürzburg/Meine Dateien/1_Auswertung/4_2__PolaritonCondensation/4_2_3_Variation_Excitonic_Fraction/Modetomographies/b3/energy1p60944-1p61046.txt"

saveFolderPath_daten = "/Users/Jessica/Library/CloudStorage/OneDrive-UniversitätWürzburg/Meine Dateien/0_Messdaten/2022-03-23_C4892_Kagome_Powerseries/A41_b3_f40_Powerserie/"
saveFolderPath_auswertung = "/Users/Jessica/Library/CloudStorage/OneDrive-UniversitätWürzburg/Meine Dateien/1_Auswertung/4_2__PolaritonCondensation/4_2_2_Ground_state_condensation/"
filename = "Hyperspectral_and_Tomo_above_threshold"

#%%create custom colormap
palette = copy(plt.get_cmap('viridis_r'))
palette.set_under('white', 1.0)  # 1.0 represents not transparent

# create custom colormap from colorlist
colorlist = [(1,1,1), (62/255, 14/255, 56/255), (103/255, 24/255, 101/255), (248/255, 101/255, 46/255), (255/255, 255/255, 220/255)]  # R -> G -> B
n_bins = 200  # Discretizes the interpolation into bins
cmap_name = 'my_list'
colormapCustom = LinearSegmentedColormap.from_list(cmap_name, colorlist, N=n_bins)

colorlist_interferogram = [(225/255,237/255,247/255),(33/255, 62/255, 106/255),(103/255,28/255,80/255),(155/255,43/255,121/255),(248/255,228/255,242/255)]
n_bins_interferogram = 500
cmap_name_interferogram = 'cmap_interferogram'
colormap_interferogram = LinearSegmentedColormap.from_list(cmap_name_interferogram, colorlist_interferogram, N=n_bins_interferogram)

colorlist_tomo = [(255/255,255/255,255/255),(43/255, 0/255, 0/255),(155/255,0/255,0/255),(255/255,35/255,0/255),(255/255,139/255,0/255),(255/255,255/255,34/255),(255/255,255/255,219/255)]
n_bins_tomo = 500
cmap_name_tomo = 'cmap_tomo'
colormap_tomo = LinearSegmentedColormap.from_list(cmap_name_tomo, colorlist_tomo, N=n_bins_tomo)


#%% --- option to limit the energy range ---
selectedlowerEnergyLimit = 1.606
selectedupperEnergyLimit = 1.614
selectedlowerkLimit = -3.0
selectedupperkLimit = 3.0
# ----------------------------------------
selectedYmax = 30
selectedYmin = 10
selectedXmax = 30
selectedXmin = 10

#%% plot
# dispersion
dispersionData_abovePth = get_intensity_k_energy_of_dispersion_txt(path_disperionAbovePth)
data_abovePth_dispersion = slice_dispersionData(dispersionData_abovePth['IntensityValues'], dispersionData_abovePth['KValues'], dispersionData_abovePth['EnergyValues'],selectedlowerEnergyLimit, selectedupperEnergyLimit, selectedlowerkLimit, selectedupperkLimit)
# show_dispersion_data(data_abovePth_dispersion['Data'], data_abovePth_dispersion['Extent'])

# tomo
tomoData_abovePth = get_values_of_isoEnergyCut_txt(path_tomoAbovePth)
# show_dispersion_data(tomoData_abovePth['IntensityValues'], tomoData_abovePth['Extent'])
data_abovePth_tomo = slice_dispersionData(tomoData_abovePth['IntensityValues'], tomoData_abovePth['xkValues'], tomoData_abovePth['ykValues'], selectedYmin, selectedYmax, selectedXmin, selectedXmax)
# show_dispersion_data(data_abovePth_tomo['Data'], data_abovePth_tomo['Extent'])

kDiff = selectedupperkLimit-selectedlowerkLimit
energyDiff = selectedupperEnergyLimit-selectedlowerEnergyLimit
aspect_dispersion = kDiff/energyDiff 

aspect_tomo = 1
    

# dispersion plot
plt.imshow(data_abovePth_dispersion['Data'], cmap=colormapCustom, aspect=aspect_dispersion, extent=data_abovePth_dispersion['Extent'])
ax = plt.gca();
forceAspect(ax,aspect=1.0)
ax.tick_params(which="major", direction='in')#, width=1.3)
ax.tick_params(which="minor", direction='in')#, width=1.3)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xlabel('$k_\mathregular{y}$ ($\mu \mathregular{m}^{-1}$)')
plt.ylabel('Energy (eV)')
v1 = np.linspace(0, 1, 6, endpoint=True)
cbar = plt.colorbar(pad = 0.02, ticks=v1)
cbar.set_label('Intensity (arb. u.)', fontsize=12, rotation=90, labelpad=-40, color="white")
if saveFigures:    
    plt.savefig(saveFolderPath_auswertung + '4_2_2_hyperspectral_above_threshold' + '.png', format='png', bbox_inches='tight', dpi=1000)
plt.show()

# tomo plot
plt.imshow(data_abovePth_tomo['Data'], cmap=colormap_tomo, aspect=aspect_tomo, extent=[0,20,0,20])
ax = plt.gca();
forceAspect(ax,aspect=1.0)
ax.tick_params(which="major", direction='in')#, width=1.3)
ax.tick_params(which="minor", direction='in')#, width=1.3)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xlabel('y ($\mu \mathregular{m}$)')
plt.ylabel('x ($\mu \mathregular{m}$)')
v1 = np.linspace(0, 1, 6, endpoint=True)
cbar = plt.colorbar(pad = 0.02, ticks=v1)
cbar.set_label('Intensity (arb. u.)', fontsize=12, rotation=90, labelpad=-40, color="white")
if saveFigures:    
    plt.savefig(saveFolderPath_auswertung + '4_2_2_modetomography_above_threshold' + '.png', format='png', bbox_inches='tight', dpi=3000)
plt.show()

