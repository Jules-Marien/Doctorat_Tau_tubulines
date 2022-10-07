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



u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')




#Requires to select the segid to not select water and ions for some reason

#beta1_chainA = "segid PROA and resid 429:444"


selection_CTT_chainA = "segid PROA and resid 427:444 and name CA"

selection_CTT_chainB= "segid PROB and resid 883:895 and name CA"

selection_CTT_chainC = "segid PROC and resid 1322:1339 and name CA"






def dist(array_1, array_2):
    return np.sqrt( (array_1[0]-array_2[0])**2 + (array_1[1]-array_2[1])**2 + (array_1[2]-array_2[2])**2 )




def distances_min_CTT_CTT(u,selection_CTT_1, selection_CTT_2):
    """Return the centers of mass for each frame"""
    
    all_dist_CTT_CTT = np.empty((len(u.trajectory), len(u.select_atoms(selection_CTT_1).positions[:,0]) , len(u.select_atoms(selection_CTT_2).positions[:,0])))  
    
    
    
    dist_min_CTT_CTT = np.empty((len(u.trajectory), 1))
    
    
    
    

	
    
    for t in range(len(u.trajectory)):
        
        #Update trajectory to each frame
        u.trajectory[t]
    
        print(t)
        
        coordinates_CTT_1 = u.select_atoms(selection_CTT_1).positions

        coordinates_CTT_2 = u.select_atoms(selection_CTT_2).positions 
        
        
        for resid_1 in range(len(coordinates_CTT_1)):

            for resid_2 in range(len(coordinates_CTT_2)):
          
                all_dist_CTT_CTT[t,resid_1,resid_2] = dist(coordinates_CTT_1[resid_1,:], coordinates_CTT_2[resid_2,:])


        
        
        #Calculate distances between 
        
        dist_min_CTT_CTT[t,0] = np.min(all_dist_CTT_CTT[t,:,:])
	
        
    return dist_min_CTT_CTT





dist_min_CTT_CTT_chainA_chain_B = distances_min_CTT_CTT(u, selection_CTT_chainA, selection_CTT_chainB)

dist_min_CTT_CTT_chainB_chain_C = distances_min_CTT_CTT(u, selection_CTT_chainB, selection_CTT_chainC)

dist_min_CTT_CTT_chainA_chain_C = distances_min_CTT_CTT(u, selection_CTT_chainA, selection_CTT_chainC)




np.savetxt('dist_min_CTT_CTT_chainA_chainB_modele'+modele+'.txt',dist_min_CTT_CTT_chainA_chain_B)
np.savetxt('dist_min_CTT_CTT_chainB_chainC_modele'+modele+'.txt',dist_min_CTT_CTT_chainB_chain_C)
np.savetxt('dist_min_CTT_CTT_chainA_chainC_modele'+modele+'.txt',dist_min_CTT_CTT_chainA_chain_C)


#%%




"""Plots"""

distance_bases = 50



#ChainA vs chainB

plt.figure()
plt.title('Minimum distance between the Cα of the CTTs  \n of the chains A (βI) and B (αI) in presence of R2 for model '+modele)

plt.xlabel('Time (ns)')
plt.ylabel('Minimum distance (Å)')

plt.ylim(0,80)
plt.xlim(0,200)


len_traj = len(dist_min_CTT_CTT_chainA_chain_B)

time = np.arange(0,len_traj,1)/10

plt.plot(time, np.ones(len_traj)*distance_bases,color='red', label="Distance at the bases")
plt.plot(time, dist_min_CTT_CTT_chainA_chain_B, label ="Distances")

plt.legend()

plt.savefig('dist_min_CTT_CTT_chainA_chain_B_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_min_CTT_CTT_chainA_chain_B, columns = ['Distances'])

plt.figure()

plt.xlim(0,80)

Distances_pandas.Distances.plot.density()
plt.axvline(x=distance_bases, color='r', linestyle='-',label ="Distance at the bases")
plt.title('Density plot of the minimum distance between the Cα of the CTTs  \n of the chains A (βI) and B (αI) in presence of R2 for model '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_dist_min_CTT_CTT_chainA_chain_B_modele'+modele+'.png')

plt.show()

















#ChainB vs chainC

plt.figure()
plt.title('Minimum distance between the Cα of the CTTs  \n of the chains B (αI) and C (βI) in presence of R2 for model '+modele)

plt.xlabel('Time (ns)')
plt.ylabel('Minimum distance (Å)')

plt.ylim(0,80)
plt.xlim(0,200)


len_traj = len(dist_min_CTT_CTT_chainB_chain_C)

time = np.arange(0,len_traj,1)/10

plt.plot(time, np.ones(len_traj)*distance_bases,color='red', label="Distance at the bases")
plt.plot(time, dist_min_CTT_CTT_chainB_chain_C, label ="Distances")

plt.legend()

plt.savefig('dist_min_CTT_CTT_chainB_chain_C_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_min_CTT_CTT_chainB_chain_C, columns = ['Distances'])

plt.figure()

plt.xlim(0,80)

Distances_pandas.Distances.plot.density()
plt.axvline(x=distance_bases, color='r', linestyle='-',label ="Distance at the bases")
plt.title('Density plot of the minimum distance between the Cα of the CTTs  \n of the chains B (αI) and C (βI) the presence of R2 for model '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_dist_min_CTT_CTT_chainB_chain_C_modele'+modele+'.png')

plt.show()















#ChainA vs chainB

plt.figure()
plt.title('Minimum distance between the Cα of the CTTs  \n of the chains A (βI) and C (βI) in presence of R2 for model '+modele)

plt.xlabel('Time (ns)')
plt.ylabel('Minimum distance (Å)')

plt.ylim(0,80)
plt.xlim(0,200)


len_traj = len(dist_min_CTT_CTT_chainA_chain_C)

time = np.arange(0,len_traj,1)/10

#plt.plot(time, np.ones(len_traj)*distance_bases,color='red', label="Distance at the bases")
plt.plot(time, dist_min_CTT_CTT_chainA_chain_C, label ="Distances")

plt.legend()

plt.savefig('dist_min_CTT_CTT_chainA_chain_C_modele'+modele+'.png')

plt.show()





Distances_pandas = pd.DataFrame(dist_min_CTT_CTT_chainA_chain_C, columns = ['Distances'])

plt.figure()

plt.xlim(0,80)

Distances_pandas.Distances.plot.density()
#plt.axvline(x=distance_bases, color='r', linestyle='-',label ="Distance at the bases")
plt.title('Density plot of the minimum distance between the Cα of the CTTs  \n of the chains A (βI) and C (βI) in presence of R2 for model '+modele)
plt.xlabel('Distances (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_dist_min_CTT_CTT_chainA_chain_C_modele'+modele+'.png')

plt.show()




