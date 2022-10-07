#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 17:48:57 2022

@author: marien
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 


from MDAnalysis.analysis import contacts

import pandas as pd

from numba import jit


#largely inspired from https://userguide.mdanalysis.org/1.1.1/examples/analysis/distances_and_contacts/contacts_within_cutoff.html



modele = '17'

u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')




"""User part"""

radius_cutoff = 5.



selection_CTT_1 = "segid PROA and resid 427:444 and not (name H*)"

selection_CTT_2 = "segid PROB and resid 883:895 and not (name H*)"








"""Program"""


"""Groupe 2 corresponds to the tubuline complex !"""

#Set the trajectory to the first frame in order to select the native contacts
u.trajectory[0]
groupe_2 = u.select_atoms(selection_CTT_2)


#Set back the trajectory to the entirety of it
u.trajectory[:]
groupe_1 = u.select_atoms(selection_CTT_1)






#%%


groupe_1_resid_numbers = groupe_1.resnums

groupe_1_resnames = groupe_1.resnames

groupe_1_indices_atoms = groupe_1.indices




#get rid of duplicates
groupe_1_filtered_resid_numbers = np.unique(groupe_1_resid_numbers)


#Table of the indices at which the residues change
index = 0
groupe_1_indices_change = np.array([0])
for i in range(1,len(groupe_1_resnames)):
    if groupe_1_resid_numbers[i-1] != groupe_1_resid_numbers[i]:
        groupe_1_indices_change = np.append(groupe_1_indices_change, i)
        
#Add the last index by hand
groupe_1_indices_change = np.append(groupe_1_indices_change, len(groupe_1_resid_numbers)-1)


#Table of the sequence of residues in the groupe
groupe_1_filtered_resnames = np.array([])
for i in range(len(groupe_1_indices_change)):
    groupe_1_filtered_resnames = np.append(groupe_1_filtered_resnames, groupe_1_resnames[groupe_1_indices_change[i]])









groupe_2_resid_numbers = groupe_2.resnums

groupe_2_resnames = groupe_2.resnames

groupe_2_indices_atoms = groupe_2.indices





#get rid of duplicates
groupe_2_filtered_resid_numbers = np.unique(groupe_2_resid_numbers)


#Table of the indices at which the residues change
index = 0
groupe_2_indices_change = np.array([0])
for i in range(1,len(groupe_2_resnames)):
    if groupe_2_resid_numbers[i-1] != groupe_2_resid_numbers[i]:
        groupe_2_indices_change = np.append(groupe_2_indices_change, i)
        
#Add the last index by hand
groupe_2_indices_change = np.append(groupe_2_indices_change, len(groupe_2_resid_numbers)-1)
        

#Table of the sequence of residues in the groupe
groupe_2_filtered_resnames = np.array([])
for i in range(len(groupe_2_indices_change)):
    groupe_2_filtered_resnames = np.append(groupe_2_filtered_resnames, groupe_2_resnames[groupe_2_indices_change[i]])
    


"""We can use the indices in groupe_1_indices_change and groupe_1_indices_change to target the proper elements in the contact matrix !"""
    







#%%


#print(u.trajectory)


def contacts_within_cutoff(u, group_a, group_b, radius):
    timeseries_ncontacts = []
    timeseries_matrix_contacts = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        
        # determine which distances <= radius
        matrix_contacts = contacts.contact_matrix(dist, radius)
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        
        timeseries_ncontacts.append([ts.frame, n_contacts])
        timeseries_matrix_contacts.append([ts.frame, matrix_contacts])
    return np.array(timeseries_ncontacts),np.array(timeseries_matrix_contacts)




"""Try to see if we can accelerate the calculations using numba


#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def contacts_within_cutoff(u, group_a, group_b, radius):
    timeseries_ncontacts = []
    timeseries_matrix_contacts = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        
        # determine which distances <= radius
        matrix_contacts = contacts.contact_matrix(dist, radius)
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        
        timeseries_ncontacts.append([ts.frame, n_contacts])
        timeseries_matrix_contacts.append([ts.frame, matrix_contacts])
    return np.array(timeseries_ncontacts),np.array(timeseries_matrix_contacts)
"""


#Select all heavy atoms by removing the hydrogens with "and (not H*)"
#Requires to select the segid to not select water and ions for some reason





number_contacts, matrices_contacts = contacts_within_cutoff(u, groupe_1 , groupe_2 , radius=radius_cutoff)



number_contacts_df = pd.DataFrame(number_contacts, columns=['Frame',
                                  '# Contacts'])
number_contacts_df.head()



matrices_contacts_df = pd.DataFrame(matrices_contacts, columns=['Frame',
                                  'Matrices des contacts'])




#%%


#Take only the matrices with True or False
matrices_contacts_brut = matrices_contacts[:,1]


#Array which contains the number of contacts per frame
nombre_contacts = np.empty((len(matrices_contacts[:,1])))


#Loop over time
for t in range(len(matrices_contacts[:,1])):
    
    contact_count = 0
    
    #Loop over Tubulines residues indices (groupe 2)
    for i in range(1,len(groupe_2_indices_change)):
        
        #boundary i-1
        boundary_left_2 = groupe_2_indices_change[i-1]
        
        #boundary i
        boundary_right_2 = groupe_2_indices_change[i]
        
        
        #Loop over Tau residues indices (groupe 1)
        for j in range(1,len(groupe_1_indices_change)):
            
            #boundary j-1
            boundary_left_1 = groupe_1_indices_change[j-1]
            
            #boundary j
            boundary_right_1 = groupe_1_indices_change[j]
            
            #Tell if there's a contact happenning in the submatrix defined as (residue of Tau) x (residue of groupe_2, aka tubuline)
            
            """Answers the question : is there a contact between residue j in Tau and residue i in tubuline at time t?"""
            
            contact_per_residue_of_CTT_1 = np.any( matrices_contacts_brut[t][boundary_left_1:boundary_right_1,boundary_left_2:boundary_right_2]  )
            
            if contact_per_residue_of_CTT_1 == True:
                
                contact_count = contact_count + 1
            
    nombre_contacts[t] = contact_count




#%%

np.savetxt('nombre_total_de_contacts_CTT_CTT_chainA_chainB_modele'+modele+'_avec_tau.txt',nombre_contacts)



#%%


time = matrices_contacts[:,0] /10


"""PLOT"""


plt.figure()

plt.xlim(0,200)
plt.ylim(0,10)

plt.title('Number of contacts residue to residue \n between the CTTs of chain A and B of model '+modele+' with R2')

plt.xlabel('Time (ns)')
plt.ylabel('Number of contacts')

plt.plot(time, nombre_contacts)

plt.savefig('Plot_nbre_total_contacts_CTT_CTT_chainA_chainB_modele'+modele+'_avec_tau.png')



plt.show()

