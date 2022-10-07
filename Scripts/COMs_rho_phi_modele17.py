#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:11:04 2022

@author: marien
"""


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 
from mpl_toolkits import mplot3d



import pandas as pd

from numba import jit


#largely inspired from https://userguide.mdanalysis.org/1.1.1/examples/analysis/distances_and_contacts/contacts_within_cutoff.html



modele = '17'



u = mda.Universe( '../../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')


#Requires to select the segid to not select water and ions for some reason

#beta1_chainA = "segid PROA and resid 429:444"


selection_cores = "protein and ((resid 1:426) or (resid 445:882) or (resid 896:1321))"

selection_chainA = "segid PROA and resid 427:444"

selection_chainB= "segid PROB and resid 883:895"

selection_chainC = "segid PROC and resid 1322:1339"




def Centers_of_mass(u,selection):
    """Return the centers of mass for each frame"""
    
    centers_of_mass = np.empty((len(u.trajectory),3))
    
    for i in range(len(u.trajectory)):
        
        #Update trajectory to each frame
        u.trajectory[i]
        
        #Calculate centers of mass
        centers_of_mass[i,:] = u.select_atoms(selection).center_of_mass()
        
    return centers_of_mass
        
    
COM_chainA = Centers_of_mass(u,selection_chainA)

COM_chainB = Centers_of_mass(u,selection_chainB)

COM_chainC = Centers_of_mass(u,selection_chainC)


#%%

COM_cores = Centers_of_mass(u, selection_cores)

#np.savetxt('COM_modele17.txt', centers_of_mass)




#%%


"""
R = Distance between the two COM

Angle Phi : in the xy plane, with regards to the x axis

Angle Theta : vertical displacement in the z direction 


R = sqrt(x**2 + y**2 + z**2)

Theta = arccos(z/R)

Phi = arctan(y/x)

The arctan function is tricky, it has to be defined correctly to avoid quadrant errors

Fortunately, in the "math" module, there is the atan2(y,x) function which returns arctan(y/x) correctly



"""

import math 

def radius(x,y,z):
    return math.sqrt(x**2 + y**2 + z**2)


def Theta(z,r):
    return math.acos(z/r)

def Phi(x,y):
    return math.atan2(y, x)


def sphericals(array_COM_cores,array_COM_chain):
    """Return the values R, Phi and Theta for a single frame
    
    Phi and theta are converted in degrees"""
    
    coord = array_COM_chain - array_COM_cores
    
    R = radius(coord[0],coord[1],coord[2])
    
    theta = Theta(coord[2],R) *(180/math.pi)
    
    phi = Phi(coord[1],coord[0])  *(180/math.pi)
    
    return R,phi,theta





def pos_spherical(array_cores, array_chain):
    
    len_traj = len(array_cores[:,0])
    
    array_pos_sphericals = np.empty((3,len_traj))
    
    for i in range(len_traj):
        array_pos_sphericals[:,i] = sphericals(array_cores[i], array_chain[i])
        
    return array_pos_sphericals




spheri_COM_chainA = pos_spherical(COM_cores, COM_chainA)

spheri_COM_chainB = pos_spherical(COM_cores, COM_chainB)

spheri_COM_chainC = pos_spherical(COM_cores, COM_chainC)




np.savetxt('COM_chainA_modele'+modele+'.txt',COM_chainA)
np.savetxt('COM_chainB_modele'+modele+'.txt',COM_chainB)
np.savetxt('COM_chainC_modele'+modele+'.txt',COM_chainC)
np.savetxt('COM_cores_modele'+modele+'.txt',COM_cores)

np.savetxt('COM_spherical_chainA_modele'+modele+'.txt',spheri_COM_chainA)
np.savetxt('COM_spherical_chainB_modele'+modele+'.txt',spheri_COM_chainB)
np.savetxt('COM_spherical_chainC_modele'+modele+'.txt',spheri_COM_chainC)





#%%
#Plot with raw colors

markersize = 100
markercolor = 'orange'

plt.figure()
plt.title('Mobilité des CTT pour le modele '+modele)

plt.xlabel('φ',weight='bold')
plt.ylabel('θ',weight='bold')



#ChainA
plt.scatter(spheri_COM_chainA[1],spheri_COM_chainA[2],color='blue',label='Chaine A')

#ChainB
plt.scatter(spheri_COM_chainB[1],spheri_COM_chainB[2],color='red',label='Chaine B')
plt.scatter(spheri_COM_chainB[1,0],spheri_COM_chainB[2,0], color =markercolor,marker = 'v', s= markersize)

#ChainC
plt.scatter(spheri_COM_chainC[1],spheri_COM_chainC[2],color='black',label='Chaine C')
plt.scatter(spheri_COM_chainC[1,0],spheri_COM_chainC[2,0], color =markercolor,marker = 'v', s= markersize)

plt.scatter(spheri_COM_chainA[1,0],spheri_COM_chainA[2,0], color =markercolor,label ='Positions initiales',marker = 'v', s= markersize)

plt.legend()

plt.savefig('Plot_COM_modele'+modele+'.png')

plt.show()






#%%

#Plot with time gradient of colors

markersize = 100
markercolor = 'orange'

plt.figure()
plt.title('Mobilité des CTT pour le modele '+modele)

plt.xlabel('φ',weight='bold')
plt.ylabel('θ',weight='bold')

#Get the color maps
cm_A = plt.get_cmap('Blues')
cm_B = plt.get_cmap('Reds')
cm_C = plt.get_cmap('Greys')

#Create the corresponding colormaps
len_traj = len(spheri_COM_chainA[0])

col_A = [cm_A(float(i)/len_traj) for i in range(len_traj)]
col_B = [cm_B(float(i)/len_traj) for i in range(len_traj)]
col_C = [cm_C(float(i)/len_traj) for i in range(len_traj)]


#ChainA
plt.scatter(spheri_COM_chainA[1],spheri_COM_chainA[2],c=col_A)
plt.scatter(spheri_COM_chainA[1,0],spheri_COM_chainA[2,0],color='blue',label='Chaine A',s=25)

#ChainB
plt.scatter(spheri_COM_chainB[1],spheri_COM_chainB[2],c=col_B)
plt.scatter(spheri_COM_chainB[1,0],spheri_COM_chainB[2,0],color='red',label='Chaine B',s=25)
plt.scatter(spheri_COM_chainB[1,0],spheri_COM_chainB[2,0], color =markercolor,marker = 'v', s= markersize)


#ChainC
plt.scatter(spheri_COM_chainC[1],spheri_COM_chainC[2],c=col_C)
plt.scatter(spheri_COM_chainC[1,0],spheri_COM_chainC[2,0],color='black',label='Chaine C',s=25)
plt.scatter(spheri_COM_chainC[1,0],spheri_COM_chainC[2,0], color =markercolor,marker = 'v', s= markersize)


plt.scatter(spheri_COM_chainA[1,0],spheri_COM_chainA[2,0], color =markercolor,label ='Positions initiales',marker = 'v', s= markersize)

plt.legend()

plt.savefig('Plot_COM_gradient_modele'+modele+'.png')

plt.show()







#%%


u.trajectory[0]
structure_ref = u.select_atoms('protein')

atoms_ref =  structure_ref.atoms.positions

u.trajectory[0]


#Coloring by segnames 
#Getting rid of hydrogens for simplicity
structure_ref_chainA = u.select_atoms('protein and segid PROA and not (name H*)')

atoms_ref_chainA =  structure_ref_chainA.atoms.positions


structure_ref_chainB = u.select_atoms('protein and segid PROB and not (name H*)')

atoms_ref_chainB =  structure_ref_chainB.atoms.positions


structure_ref_chainC = u.select_atoms('protein and segid PROC and not (name H*)')

atoms_ref_chainC =  structure_ref_chainC.atoms.positions



structure_ref_chainD = u.select_atoms('protein and segid PROD and not (name H*)')

atoms_ref_chainD =  structure_ref_chainD.atoms.positions









fig = plt.figure()
ax = plt.axes(projection = '3d')

#ax.scatter3D(COM_chainA[:,0], COM_chainA[:,1], COM_chainA[:,2], color ='purple', s = 50)
ax.scatter3D(COM_chainA[:,0], COM_chainA[:,1], COM_chainA[:,2], c=col_A, s = 50)

#ax.scatter3D(COM_chainB[:,0], COM_chainB[:,1], COM_chainB[:,2], color ='green', s = 50)
ax.scatter3D(COM_chainB[:,0], COM_chainB[:,1], COM_chainB[:,2], c=col_B, s = 50)

#ax.scatter3D(COM_chainC[:,0], COM_chainC[:,1], COM_chainC[:,2], color ='grey', s = 50)
ax.scatter3D(COM_chainC[:,0], COM_chainC[:,1], COM_chainC[:,2], c=col_C, s = 50)

ax.scatter3D(COM_cores[:,0], COM_cores[:,1], COM_cores[:,2], color ='yellow', s = 50)


ax.scatter3D(atoms_ref_chainA[:,0], atoms_ref_chainA[:,1], atoms_ref_chainA[:,2], color = 'blue', s= 1)

ax.scatter3D(atoms_ref_chainB[:,0], atoms_ref_chainB[:,1], atoms_ref_chainB[:,2], color = 'red', s= 1)

ax.scatter3D(atoms_ref_chainC[:,0], atoms_ref_chainC[:,1], atoms_ref_chainC[:,2], color = 'black', s= 1)

ax.scatter3D(atoms_ref_chainD[:,0], atoms_ref_chainD[:,1], atoms_ref_chainD[:,2], color = 'orange', s= 1)


#plt.savefig("Plot_center_of_mass_test.png")                
plt.show()


