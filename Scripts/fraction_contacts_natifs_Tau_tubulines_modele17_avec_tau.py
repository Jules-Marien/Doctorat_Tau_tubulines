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



selection_fragment_tau = "segid PROD and resid 1340:1366 and not (name H*)"

#selection_beta1_chainA = "segid PROA and resid 427:444 and not (name H*)"

#selection_chain = "(around "+ str(radius_cutoff)+ " (segid PROD and resid 1340:1366 and not (name H*)) ) and not (segid PROD)"

#selection_chain = "(around "+ str(radius_cutoff)+ " (segid PROD and resid 1340:1366 and not (name H*))) and not segid PROD"

selection_chain = "(around "+ str(radius_cutoff)+ " (segid PROD)) and not (name H*) and not segid PROD"


residues_Tau = range(1340,1367)













"""Program"""


"""Groupe 2 corresponds to the tubuline complex !"""

#Set the trajectory to the first frame in order to select the native contacts
u.trajectory[0]
groupe_2 = u.select_atoms(selection_chain)


#Set back the trajectory to the entirety of it
u.trajectory[:]
groupe_1 = u.select_atoms(selection_fragment_tau)






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



"""Idea : Comparing matrices with True value when the residues are in contact"""


#Take only the matrices with True or False
matrices_contacts_brut = matrices_contacts[:,1]



def matrice_contacts_residues(matrix_of_contacts, indices_change_groupe_1, indices_change_groupe_2):
    """
    matrix_of_contacts = matrices_contacts_brut[t]
       indices_change_groupe_1 = groupe_1_indices_change
       indices_change_groupe_2 = groupe_2_indices_change
    """
    
    matrix_contact_True_or_False = np.full((len(indices_change_groupe_1)-1,len(indices_change_groupe_2)-1), False, dtype=bool)

    
    
    #Loop over Tubulines residues indices (groupe 2)
    for i in range(1,len(indices_change_groupe_2)):
        
        #boundary i-1
        boundary_left_2 = indices_change_groupe_2[i-1]
        
        #boundary i
        boundary_right_2 = indices_change_groupe_2[i]
        
        
        #Loop over Tau residues indices (groupe 1)
        for j in range(1,len(indices_change_groupe_1)):
            
            #boundary j-1
            boundary_left_1 = indices_change_groupe_1[j-1]
            
            #boundary j
            boundary_right_1 = indices_change_groupe_1[j]
            
            #Tell if there's a contact happenning in the submatrix defined as (residue of Tau) x (residue of groupe_2, aka tubuline)
            
            """Answers the question : is there a contact between residue j in Tau and residue i in tubuline at time t?"""
            
            contact_per_residue_of_Tau = np.any( matrix_of_contacts[boundary_left_1:boundary_right_1,boundary_left_2:boundary_right_2]  )
            
            if contact_per_residue_of_Tau == True:
                
                matrix_contact_True_or_False[j-1,i-1] = True
                
    return matrix_contact_True_or_False


matrix_contacts_natifs = matrice_contacts_residues(matrices_contacts_brut[0], groupe_1_indices_change, groupe_2_indices_change)
       
    


#%%


matrices_contacts_reduites = np.full((len(matrices_contacts[:,1]), len(matrix_contacts_natifs[:,0]), len(matrix_contacts_natifs[0,:])), False, dtype=bool)



for t in range(len(matrices_contacts[:,1])):
    
    matrices_contacts_reduites[t,:,:] = matrice_contacts_residues(matrices_contacts_brut[t], groupe_1_indices_change, groupe_2_indices_change)






#%%

"""Comparison of native state and the trajectory"""


matrices_comparaisons = np.full((len(matrices_contacts[:,1]), len(matrix_contacts_natifs[:,0]), len(matrix_contacts_natifs[0,:])), False, dtype=bool)

matrices_comparaisons[0,:,:] = matrix_contacts_natifs
    
for t in range(1,len(matrices_contacts_reduites[:,0,0])):
    
    for i in range(len(matrix_contacts_natifs[:,0])):
        
        for j in range(len(matrix_contacts_natifs[0,:])):
            
            if matrices_contacts_reduites[0,i,j] == matrices_contacts_reduites[t,i,j] and matrices_contacts_reduites[0,i,j] == True:
            
                matrices_comparaisons[t,i,j] = True
        
        



#%%



"""Calculation of fractions of native contacts"""

"""
#The number of preserved native contacts between pairs of residues is the number of True occurence in one matrix_contacts_reduites !

fractions_native_contacts = np.empty((len(matrices_contacts[:,1])))

preserved_native_contacts = np.empty((len(matrices_contacts[:,1])))

preserved_native_contacts[0] = np.count_nonzero(matrices_contacts_reduites[0,:,:])

fractions_native_contacts[0] = 1.

for t in range(1,len(matrices_contacts[:,1])):
    
    preserved_native_contacts[t] = np.count_nonzero(matrices_contacts_reduites[t,:,:])
    
    fractions_native_contacts[t] = np.count_nonzero(matrices_contacts_reduites[t,:,:]) / np.count_nonzero(matrices_contacts_reduites[0,:,:]) 

"""




#The number of preserved native contacts between pairs of residues is the number of True occurence in one matrix_comparaison !

fractions_native_contacts = np.empty((len(matrices_contacts[:,1])))

preserved_native_contacts = np.empty((len(matrices_contacts[:,1])))

preserved_native_contacts[0] = np.count_nonzero(matrices_comparaisons[0,:,:])

fractions_native_contacts[0] = 1.

for t in range(1,len(matrices_contacts[:,1])):
    
    preserved_native_contacts[t] = np.count_nonzero(matrices_comparaisons[t,:,:])
    
    fractions_native_contacts[t] = np.count_nonzero(matrices_comparaisons[t,:,:]) / np.count_nonzero(matrices_comparaisons[0,:,:]) 







#%%

np.savetxt('fraction_contacts_natifs_modele'+modele+'.txt',fractions_native_contacts)



#%%


time = matrices_contacts[:,0] /10


"""PLOT"""


plt.figure()

plt.title('Fraction de contacts natifs préservés \n résidu à résidu entre le fragment de Tau et les tubulines \n du modele '+modele+' au cours de la trajectoire')

plt.xlabel('Temps (ns)')
plt.ylabel('Fraction de contacts natifs')

plt.plot(time, fractions_native_contacts)

plt.savefig('fraction_contacts_natifs_Tau_tubulines_modele'+modele+'_avec_tau.png')

plt.show()













#%%





def dataframe(matrice,groupe_1,groupe_2):
    list_names_groupe_1 = []
    list_names_groupe_2 = []
    
    for i in range(len(groupe_1)):
        list_names_groupe_1.append(str(groupe_1[i]))
        
    for j in range(len(groupe_2)):
        list_names_groupe_2.append(str(groupe_2[j]))
        
    matrice_contact_complete_df = pd.DataFrame(matrice[0,1], columns = list_names_groupe_2 , index = list_names_groupe_1)
    
    return matrice_contact_complete_df

         





def contact_pairs(matrix,number_contacts):
    """Return the infos about the pairs of contacts in matrix produced as matrices_contacts""" 
    
    #Create the array of names
    array_names_groupe_1 = np.array([],dtype=object)
    array_names_groupe_2 = np.array([],dtype=object)
    
    for i in range(len(groupe_1)):
        array_names_groupe_1 = np.append(array_names_groupe_1, str(groupe_1[i]))
        
    for j in range(len(groupe_2)):
        array_names_groupe_2 = np.append(array_names_groupe_2, str(groupe_2[j]))



    nbre_frames = len(matrix[:,0])
    
    array_contacts_name = np.empty((nbre_frames,2),dtype=object)
    
    #Run on all frames
    for i in range(nbre_frames):
        
        #put the number of the frame:
        array_contacts_name[i,0] = i
        
        
        #Sets the array which will contain the names of the contacts in the i frame
        nbre_contacts = number_contacts[i,1]
        
        array_contacts = np.empty((nbre_contacts,2),dtype=object)
        
        contact = 0        
        #Run on all elements of the matix
        for j in range(len(matrix[i][1][:,0])):
            for k in range(len(matrix[i][1][0,:])):
                
                if matrix[i][1][j,k] == True:
                    array_contacts[contact,0] = array_names_groupe_1[j]
                    array_contacts[contact,1] = array_names_groupe_2[k]
                    
                    #allows to move in the array_contacts as the contacts are detected
                    contact = contact + 1
                    
                    
        array_contacts_name[i,1] = array_contacts

        
    
                    
    return array_contacts_name



                    



#Library containing the findall() function allowing to split by words
import re


def array_names_residues_contact(array_contact_pairs):
    """return an array with all the residues numbers involved in a contact"""

    nbre_frames = len(array_contact_pairs[:,0])

    array_contact_numbers = np.empty((nbre_frames,2),dtype=object)

    for i in range(nbre_frames):
        array_contact_numbers[i,0] = i
        
        #how many contacts in the i frame
        nbre_contacts = len(array_contact_pairs[i,1])
        
        array_numbers = np.empty((nbre_contacts,2),dtype=object)
        
        #Looping on all the contacts
        for j in range(nbre_contacts):
            
            #Splitting in words with findall, the resid number is in 11th position
            resid_number_1 = re.findall(r'\w+' , array_contact_pairs[i,1][j][0])[10]
            resid_number_2 = re.findall(r'\w+' , array_contact_pairs[i,1][j][1])[10]
    
            
            array_numbers[j,0] = resid_number_1
            array_numbers[j,1] = resid_number_2
            
            
        array_contact_numbers[i,1] = array_numbers
        
    return array_contact_numbers









def array_involved_residues(array_contact_numbers):
    """return an array with the non repeating residue numbers involved in a contact for each frame"""
    
    nbre_frames = len(array_contact_numbers[:,0])
    
    array_involved_residues = np.empty((nbre_frames,2), dtype=object)
    
    for i in range(nbre_frames):
        array_involved_residues[i,0] = i 
        
        
        #Sort the array 
        array_numbers = np.sort(array_contact_numbers[i,1])
        
        #Keeps only one occurance of each mentionned residue
        array_numbers = np.unique(array_numbers)
    
        
        array_involved_residues[i,1] = array_numbers
        
        
    return array_involved_residues




            
            
    
  
list_names_groupe_1 = []
list_names_groupe_2 = []

for i in range(len(groupe_1)):
    list_names_groupe_1.append(str(groupe_1[i]))
    
for j in range(len(groupe_2)):
    list_names_groupe_2.append(str(groupe_2[j]))
        
    

def array_contacts_per_residue(matrix_contacts,res_interest_number,selection,list_names_1,list_names_2):
    """Return an array of size "time" with True if the residue of interest at res_interest_number is in contact with the selection at time t and False otherwise
    
       res_interest_number is in list_names_2
       selection is in list_names_1"""
    
    frames = len(matrix_contacts[:,0])
    
    #all the occurences of the residue of interest in list_names_2 == all its heavy atoms
    range_resid = np.array([])
    
    for i in range(len(list_names_2)):
        residue_number = int(re.findall(r'\w+' ,list_names_2[i])[10])
        
        if residue_number == res_interest_number:
            range_resid = np.append(range_resid,i)
            
    range_resid = np.array(range_resid , dtype=int)    
            
    
    
    #all the occurences of the selection in list_names_1 == all heavy atoms of the selection
    
    range_selec = np.array([])
    
    #loop over all the selection
    for resid_number_selection in selection:
        
        for i in range(len(list_names_1)):
            residue_number = int(re.findall(r'\w+' ,list_names_1[i])[10])
            
            if residue_number == resid_number_selection:
                range_selec = np.append(range_selec,i)
                
                
    range_selec = np.array(range_selec ,  dtype=int)
    

    
    
    #All contacts are set to False to begin
    array_contacts_over_time = np.zeros((frames), dtype=bool)
    
    
    #loop over all heavy atoms of the residue of interest
    for resid in range_resid:
        
        #loop over time 
        for t in range(frames):
            
            #loop over all elements of the selection:
            for selec in range_selec:
                
                if matrix_contacts[t,1][selec,resid] == True:
                    array_contacts_over_time[t] = True
                
    return array_contacts_over_time


#res_interest_number corresponds to the CTT residue
 
#selection are the residues in the Tau fragment

  
        
    
    
