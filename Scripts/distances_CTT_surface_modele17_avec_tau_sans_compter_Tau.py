#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 13:25:04 2022

@author: marien
"""


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 

import MDAnalysis.analysis.distances as mdaDist


import pandas as pd



modele = '17'



u = mda.Universe( '../../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')




#Requires to select the segid to not select water and ions for some reason

#beta1_chainA = "segid PROA and resid 429:444"


selection_cores = "protein and ((resid 1:426) or (resid 445:882) or (resid 896:1321)) and name CA"

selection_end_chainA = "segid PROA and resid 444 and name CA"

selection_end_chainB= "segid PROB and resid 895 and name CA"

selection_end_chainC = "segid PROC and resid 1339 and name CA"






def dist(array_1, array_2):
    return np.sqrt( (array_1[0]-array_2[0])**2 + (array_1[1]-array_2[1])**2 + (array_1[2]-array_2[2])**2 )




def distances_CTT_surface(u,selection_surface, selection_CTT_end):
    """Return the centers of mass for each frame"""
    
    all_dist_CTT_surface = np.empty((len(u.trajectory), len(u.select_atoms(selection_surface).positions[:,0]) ))  
    
    
    
    dist_CTT_surface = np.empty((len(u.trajectory), 1))
    
    
    
    

	
    
    for t in range(len(u.trajectory)):
        
        #Update trajectory to each frame
        u.trajectory[t]
    
        print(t)
        
        coordinates_end_CTT = u.select_atoms(selection_CTT_end).positions[0,:]
        
        
        coordinates_surface = u.select_atoms(selection_surface).positions
    
        
        
        
        for atom in range(len(u.select_atoms(selection_surface).positions[:,0])):
            
        
            
            all_dist_CTT_surface[t,atom] = dist(coordinates_surface[atom,:], coordinates_end_CTT ) 
        
        
        #Calculate distances between 
        
        dist_CTT_surface[t,0] = np.min(all_dist_CTT_surface[t,:])
	
        
    return dist_CTT_surface





dist_CTT_surface_chainA = distances_CTT_surface(u, selection_cores, selection_end_chainA)

dist_CTT_surface_chainB = distances_CTT_surface(u, selection_cores, selection_end_chainB)

dist_CTT_surface_chainC = distances_CTT_surface(u, selection_cores, selection_end_chainC)






np.savetxt('dist_CTT_surface_chainA_sans_compter_tau_modele'+modele+'.txt',dist_CTT_surface_chainA)
np.savetxt('dist_CTT_surface_chainB_sans_compter_tau_modele'+modele+'.txt',dist_CTT_surface_chainB)
np.savetxt('dist_CTT_surface_chainC_sans_compter_tau_modele'+modele+'.txt',dist_CTT_surface_chainC)



#%%




"""Plots"""




#ChainA

plt.figure()
plt.title('Distance entre la surface et le résidu de fin de la CTT \n de la chaine A du modele '+modele+' en fonction du temps')

plt.xlabel('Temps (ns)')
plt.ylabel('Distance (Å)')

plt.ylim(0,60)


mean_value = np.mean(dist_CTT_surface_chainA)
len_traj = len(dist_CTT_surface_chainA)

time = np.arange(0,len_traj,1)/10

plt.plot(time, np.ones(len_traj)*mean_value,color='red',label = "Valeur moyenne")
plt.plot(time, dist_CTT_surface_chainA, label ="Distances")

plt.legend()

plt.savefig('dist_CTT_surface_chainA_sans_compter_tau_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_CTT_surface_chainA, columns = ['Distances'])

plt.figure()

plt.xlim(0,60)

Distances_pandas.Distances.plot.density()
plt.axvline(x=mean_value, color='r', linestyle='-',label ="Valeur moyenne")
plt.title('Graphique de densité de la distance entre la surface \n et le résidu de fin de la CTT de la chaine A du modele '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Densité')
plt.legend()

plt.savefig('Density_plot_dist_CTT_surface_chainA_sans_compter_tau_modele'+modele+'.png')

plt.show()














#ChainB

plt.figure()
plt.title('Distance entre la surface et le résidu de fin de la CTT \n de la chaine B du modele '+modele+' en fonction du temps')

plt.xlabel('Temps (ns)')
plt.ylabel('Distance (Å)')


plt.ylim(0,60)


mean_value = np.mean(dist_CTT_surface_chainB)
len_traj = len(dist_CTT_surface_chainB)

time = np.arange(0,len_traj,1)/10

plt.plot(time, np.ones(len_traj)*mean_value,color='red',label = "Valeur moyenne")
plt.plot(time, dist_CTT_surface_chainB, label ="Distances")

plt.legend()

plt.savefig('dist_CTT_surface_chainB_sans_compter_tau_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_CTT_surface_chainB, columns = ['Distances'])

plt.figure()

plt.xlim(0,60)

Distances_pandas.Distances.plot.density()
plt.axvline(x=mean_value, color='r', linestyle='-',label ="Valeur moyenne")
plt.title('Graphique de densité de la distance entre la surface \n et le résidu de fin de la CTT de la chaine B du modele '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Densité')
plt.legend()

plt.savefig('Density_plot_dist_CTT_surface_chainB_sans_compter_tau_modele'+modele+'.png')

plt.show()











#ChainC

plt.figure()
plt.title('Distance entre la surface et le résidu de fin de la CTT \n de la chaine C du modele '+modele+' en fonction du temps')

plt.xlabel('Temps (ns)')
plt.ylabel('Distance (Å)')

plt.ylim(0,60)


mean_value = np.mean(dist_CTT_surface_chainC)
len_traj = len(dist_CTT_surface_chainC)

time = np.arange(0,len_traj,1)/10

plt.plot(time, np.ones(len_traj)*mean_value,color='red',label = "Valeur moyenne")
plt.plot(time, dist_CTT_surface_chainC, label ="Distances")

plt.legend()

plt.savefig('dist_CTT_surface_chainC_sans_compter_tau_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_CTT_surface_chainC, columns = ['Distances'])

plt.figure()

plt.xlim(0,60)

Distances_pandas.Distances.plot.density()
plt.axvline(x=mean_value, color='r', linestyle='-',label ="Valeur moyenne")
plt.title('Graphique de densité de la distance entre la surface \n et le résidu de fin de la CTT de la chaine C du modele '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Densité')
plt.legend()

plt.savefig('Density_plot_dist_CTT_surface_chainC_sans_compter_tau_modele'+modele+'.png')

plt.show()




