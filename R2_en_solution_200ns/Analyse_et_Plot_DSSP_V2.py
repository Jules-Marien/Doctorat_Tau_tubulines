#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 15:43:27 2022

@author: marien
"""


import numpy as np
import matplotlib.pyplot as plt

#Load the dssp data
data_dssp = np.genfromtxt('/ibpc/tethys/marien/Documents/Analyse_Tau_seule/DSSP_traj.txt',dtype='str')

"""The idea is to create an array with specific numbers for each type of structure

The structure code is the following :
    
   1: ’H’ : Alpha helix

   2: ‘B’ : Residue in isolated beta-bridge

   3: ‘E’ : Extended strand, participates in beta ladder

   4: ‘G’ : 3-helix (3/10 helix)

   5: ‘I’ : 5 helix (pi helix)

   6: ‘T’ : hydrogen bonded turn

   7: ‘S’ : bend

   0: ‘0‘ : Loops and irregular elements
    
"""


def treatment_dssp(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ’H’ : Alpha helix  and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_treated = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'H':
                array_treated[i,j] = 1
                
                
            if data[i,j] == 'B':
                array_treated[i,j] = 2
                
            if data[i,j] == 'E':
                array_treated[i,j] = 3
                
                
            if data[i,j] == 'G':
                array_treated[i,j] = 4
                
            if data[i,j] == 'I':
                array_treated[i,j] = 5
                
                
            if data[i,j] == 'T':
                array_treated[i,j] = 6
                
                
            if data[i,j] == 'S':
                array_treated[i,j] = 7
                
                
            if data[i,j] == '0':
                array_treated[i,j] = 0
            
            
                
    return array_treated






array_dssp_treated = treatment_dssp(data_dssp)










nbre_of_frame = len(data_dssp[:,0])
nbre_of_residues = len(data_dssp[0,:])

time = np.linspace(0,200,1000)
array_residus = np.linspace(1,27,27)


"""Creation of the corresponding colormap"""
#Inspired by : https://riptutorial.com/matplotlib/example/20692/custom-discrete-colormap

from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm

#List of the colors
map_perso = ListedColormap(['white','blue','black','red','grey','purple','yellow','green'])

#List of the boundaries of values 
boundaries = [-0.5, 0.5, 1.5 , 2.5 , 3.5 , 4.5 , 5.5, 6.5 , 7.5]

norm = BoundaryNorm(boundaries, map_perso.N, clip =True)

#the map will take color i between boundaries i and i+1


plt.figure()
plt.pcolormesh(time,array_residus,np.transpose(array_dssp_treated),cmap=map_perso,norm=norm,shading='auto')


"""Legend"""


#inpired by : https://www.geeksforgeeks.org/how-to-manually-add-a-legend-with-a-color-box-on-a-matplotlib-figure/
import matplotlib.patches as mpatches

patch_coil = mpatches.Patch(color='white',label = 'Coil')
patch_b_sheet = mpatches.Patch(color='red',label = 'B-sheet')
patch_b_bridge = mpatches.Patch(color='black',label = 'B-bridge')
patch_bend = mpatches.Patch(color='green',label='Bend')
patch_turn = mpatches.Patch(color='yellow',label='Turn')
patch_alpha_helix = mpatches.Patch(color='blue',label='A-Helix')
patch_3_helix = mpatches.Patch(color='grey',label='3-Helix')
patch_5_helix = mpatches.Patch(color='purple',label='5-Helix')


#bbox_to_anchor = allows to place the legend where we wish
#ncol = specifiy the number of columns we want in our legend

plt.legend(edgecolor='black',handles = [patch_coil,patch_b_sheet,patch_b_bridge,patch_bend,patch_turn,patch_alpha_helix,patch_3_helix,patch_5_helix],bbox_to_anchor =(1,-0.2),ncol=4)

#Force the edges of the patches to become black
ax = plt.gca()
leg = ax.get_legend()

for i in range(8):
    leg.legendHandles[i].set_edgecolor('black')



plt.xlabel('Temps (ns)')
plt.ylabel('Residu')


#bbox_inches : allows to keep the legend in the saved figure
plt.savefig('Plot_DSSP_tau_seule.png',bbox_inches='tight',dpi=1200)



"""The idea is to create an array with specific numbers for each type of structure

The structure code is the following :
    
   1:  blue   : ’H’ : Alpha helix

   2:  black  : ‘B’ : Residue in isolated beta-bridge

   3:  red    : ‘E’ : Extended strand, participates in beta ladder

   4:  grey   : ‘G’ : 3-helix (3/10 helix)

   5:  purple : ‘I’ : 5 helix (pi helix)

   6:  yellow : ‘T’ : hydrogen bonded turn

   7:  green  : ‘S’ : bend

   0:  white  : ‘0‘ : Loops and irregular elements
    
"""
