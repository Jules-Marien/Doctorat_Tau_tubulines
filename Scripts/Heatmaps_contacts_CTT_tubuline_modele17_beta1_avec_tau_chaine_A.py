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

chaine ='A'

u = mda.Universe( '../../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')




"""User part"""

radius_cutoff = 5.



selection_chaine_CTT = "segid PROA and resid 427:444 and not (name H*)"

#selection_beta1_chainA = "segid PROA and resid 427:444 and not (name H*)"

selection_contacts_natifs_de_Tau = "around "+ str(radius_cutoff)+ " (segid PROD and resid 1340:1366 and not (name H*))"















"""Program"""

#Set the trajectory to the first frame in order to select the native contacts
u.trajectory[0]
groupe_2 = u.select_atoms(selection_contacts_natifs_de_Tau)


#Set back the trajectory to the entirety of it
u.trajectory[:]
groupe_1 = u.select_atoms(selection_chaine_CTT)




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
    
    array_contacts_per_residue first takes the contact matrix, 
    then an integer (residue number of the residue in groupe 1), 
    then a LIST of integer(s) (residue number of the residue in groupe 2),
    then the lists of the caracterisitics of the selections
    
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

  
        
    
    
#%%  
    
    
"PLOTS"
    
    
    
    
time = np.linspace(0, 0.1*len(matrices_contacts[:,0]), len(matrices_contacts[:,0]))     
    


#Create the lists of residues names and numbers for groupe1
residues_names_groupe1_init = []

#loop over all the selection
for i in range(len(list_names_groupe_1)):
    
    residue_name_groupe1 = str(re.findall(r'\w+' ,list_names_groupe_1[i])[10]+ "_"+  re.findall(r'\w+' ,list_names_groupe_1[i])[8])

    residues_names_groupe1_init.append(residue_name_groupe1)
    

#get rid of duplicates
residues_names_groupe1 =[] 
for i in residues_names_groupe1_init:
    if i not in residues_names_groupe1:
        residues_names_groupe1.append(i)



#Create the list of residue numbers as string for the ticks

ticks_groupe1_init = []

#loop over all the selection
for i in range(len(list_names_groupe_1)):

    residue_groupe1 = str(re.findall(r'\w+' ,list_names_groupe_1[i])[10])
        
    ticks_groupe1_init.append(residue_groupe1)
   
    
   
#get rid of duplicates
ticks_groupe1 =[] 
for i in ticks_groupe1_init:
    if i not in ticks_groupe1:
        ticks_groupe1.append(i)








#Create the lists of residues names and numbers for groupe2
residues_names_groupe2_init = []

#loop over all the selection
for i in range(len(list_names_groupe_2)):
    
    residue_name_groupe2 = str(re.findall(r'\w+' ,list_names_groupe_2[i])[10]+ "_"+  re.findall(r'\w+' ,list_names_groupe_2[i])[8])

    residues_names_groupe2_init.append(residue_name_groupe2)
    

#get rid of duplicates
residues_names_groupe2 =[] 
for i in residues_names_groupe2_init:
    if i not in residues_names_groupe2:
        residues_names_groupe2.append(i)



#Create the list of residue numbers as string for the ticks

ticks_groupe2_init = []

#loop over all the selection
for i in range(len(list_names_groupe_2)):

    residue_groupe2 = str(re.findall(r'\w+' ,list_names_groupe_2[i])[10])
        
    ticks_groupe2_init.append(residue_groupe2)
   
    
   
#get rid of duplicates
ticks_groupe2 =[] 
for i in ticks_groupe2_init:
    if i not in ticks_groupe2:
        ticks_groupe2.append(i)












#%%            
           
len_groupe1 = len(residues_names_groupe1)
len_groupe2 = len(residues_names_groupe2)





matrix_heatmap = np.empty((len_groupe1,len_groupe2))

for i in range(len_groupe1):
    
    for j in range(len_groupe2):
        
        #array_contacts_per_residue first takes the contact matrix, then an integer (residue number of the residue in groupe 1), then a LIST of integer(s) (residue number of the residue in groupe 2), 
        #then the lists of the caracterisitics of the selections
        
        array_of_contact = array_contacts_per_residue(matrices_contacts, int(ticks_groupe2[j]), [int(ticks_groupe1[i])], list_names_groupe_1, list_names_groupe_2)
        
        percentage = np.sum(array_of_contact)*100/np.size(array_of_contact)
        
        matrix_heatmap[i,j] = percentage
        
        
        

 





#%%


len_Tau = 27

import matplotlib.ticker as mticker

plt.figure(figsize =(20,12))
ax=plt.axes()

plt.grid(True)

#Careful to the dimensions here
plt.imshow(matrix_heatmap[:,:-len_Tau], cmap='Purples',vmin = 0, vmax =100)

plt.title('Pourcentage de contacts observés entre résidus de la CTT de la chaine'+chaine+' \n et ceux des coeurs de tubuline en contact natif avec le fragment Tau \n du modele '+modele+' sur la durée de la trajectoire', fontsize=13,weight='bold')

plt.colorbar().set_label(label='% de contacts',fontsize=13,weight='bold')


#Setting ticks 

#imshow automatically place the biggest length as the y axis, so we need to adapt for the ticks

#if 

ax.xaxis.set_major_locator(mticker.FixedLocator(range(len(ticks_groupe2))))
ax.xaxis.set_major_formatter(mticker.FixedFormatter(residues_names_groupe2))

ax.yaxis.set_major_locator(mticker.FixedLocator(range(len(ticks_groupe1))))
ax.yaxis.set_major_formatter(mticker.FixedFormatter(residues_names_groupe1))
    

plt.xlabel('Résidus des coeurs de tubuline en contact natif avec Tau')
plt.ylabel('Résidus de la CTT de la chaine'+chaine)
    
plt.xticks(rotation=90)

plt.savefig('Heatmap_contacts_CTT_tubulines_modele'+modele+'_chaine'+chaine+'_test2.png')

plt.show()






