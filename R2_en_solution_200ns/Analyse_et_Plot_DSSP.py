#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 11:11:50 2022

@author: marien
"""

import numpy as np
import matplotlib.pyplot as plt

#Load the dssp data
data_dssp = np.genfromtxt('/ibpc/tethys/marien/Documents/Analyse_Tau_seule/DSSP_traj.txt',dtype='str')

"""The idea is to create arrays for each type of structure, of format residue x frame, with value 1 if the residue is of the considered structure at 
this frame, and 0 if it isn't

The structure code is the following :
    
    ’H’ : Alpha helix

    ‘B’ : Residue in isolated beta-bridge

    ‘E’ : Extended strand, participates in beta ladder

    ‘G’ : 3-helix (3/10 helix)

    ‘I’ : 5 helix (pi helix)

    ‘T’ : hydrogen bonded turn

    ‘S’ : bend

    ‘0‘ : Loops and irregular elements
    
"""


def alpha_helix(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ’H’ : Alpha helix  and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_alpha_helix = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'H':
                array_alpha_helix[i,j] = 1
            
            else :
                array_alpha_helix[i,j] = 0
                
    return array_alpha_helix


def beta_bridge(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘B’ : Residue in isolated beta-bridge  and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_beta_bridge = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'B':
                array_beta_bridge[i,j] = 1
            
            else :
                array_beta_bridge[i,j] = 0
                
    return array_beta_bridge



def extended_strand(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘E’ : Extended strand, participates in beta ladder and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_extended_strand = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'E':
                array_extended_strand[i,j] = 1
            
            else :
                array_extended_strand[i,j] = 0
                
    return array_extended_strand



def G_3_helix(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘G’ : 3-helix (3/10 helix) and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_G_3_helix = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'G':
                array_G_3_helix[i,j] = 1
            
            else :
                array_G_3_helix[i,j] = 0
                
    return array_G_3_helix




def I_5_helix(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘I’ : 5 helix (pi helix) and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_I_5_helix = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'I':
                array_I_5_helix[i,j] = 1
            
            else :
                array_I_5_helix[i,j] = 0
                
    return array_I_5_helix





def hydro_bonded_turn(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘T’ : hydrogen bonded turn and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_hydro_bonded_turn = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'T':
                array_hydro_bonded_turn[i,j] = 1
            
            else :
                array_hydro_bonded_turn[i,j] = 0
                
    return array_hydro_bonded_turn



def bend(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘S’ : bend and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_bend = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == 'S':
                array_bend[i,j] = 1
            
            else :
                array_bend[i,j] = 0
                
    return array_bend





def loops_and_irregular(data):
    """Take the dssp data numpy array and returns an array of the same shape,
       with 1s if the code is ‘0‘ : Loops and irregular elements and 0s otherwise"""
       
    nbre_of_frame = len(data[:,0])
    nbre_of_residues = len(data[0,:])
    
    array_loops_and_irregular = np.empty((nbre_of_frame,nbre_of_residues))
    
    for i in range(nbre_of_frame):
        for j in range(nbre_of_residues):
            
            if data[i,j] == '0':
                array_loops_and_irregular[i,j] = 1
            
            else :
                array_loops_and_irregular[i,j] = 0
                
    return array_loops_and_irregular


array_alpha_helix = alpha_helix(data_dssp)
array_beta_bridge = beta_bridge(data_dssp)
array_extended_strand = extended_strand(data_dssp)
array_3_helix = G_3_helix(data_dssp)
array_I_5_helix = I_5_helix(data_dssp)
array_hydro_bonded_turn = hydro_bonded_turn(data_dssp)
array_bend = bend(data_dssp)
array_loops_and_irregular = loops_and_irregular(data_dssp)









nbre_of_frame = len(data_dssp[:,0])
nbre_of_residues = len(data_dssp[0,:])

time = np.linspace(0,200,1000)
array_residus = np.linspace(1,27,27)


"""Creation of the corresponding colormaps"""
from matplotlib.colors import ListedColormap

#white if =0, colored if =1

cm_alpha_helix = ListedColormap(['white','blue']) 
cm_beta_bridge = ListedColormap(['white','black'])
cm_extended_strand = ListedColormap(['white','red'])
cm_3_helix = ListedColormap(['white','grey'])
cm_I_5_helix = ListedColormap(['white','purple'])
cm_hydro_bonded_turn = ListedColormap(['white','yellow'])
cm_bend = ListedColormap(['white','green'])

#Irregular structures left white
cm_loops_and_irregular = ListedColormap(['white','white'])



plt.figure()

#plt.pcolormesh(time,array_residus,np.transpose(array_alpha_helix),cmap=cm_alpha_helix)


#plt.pcolormesh(time,array_residus,np.transpose(array_beta_bridge),cmap=cm_beta_bridge)


#plt.pcolormesh(time,array_residus,np.transpose(array_extended_strand),cmap=cm_extended_strand)


#plt.pcolormesh(time,array_residus,np.transpose(array_3_helix),cmap=cm_3_helix)
#plt.pcolormesh(time,array_residus,np.transpose(array_I_5_helix),cmap=cm_I_5_helix)
#plt.pcolormesh(time,array_residus,np.transpose(array_hydro_bonded_turn),cmap=cm_hydro_bonded_turn)
#plt.pcolormesh(time,array_residus,np.transpose(array_bend),cmap=cm_bend)

plt.pcolormesh(time,array_residus,np.transpose(array_alpha_helix),cmap=cm_alpha_helix)


plt.pcolormesh(time,array_residus,np.transpose(array_beta_bridge),cmap=cm_beta_bridge)


plt.pcolormesh(time,array_residus,np.transpose(array_extended_strand),cmap=cm_extended_strand)


plt.pcolormesh(time,array_residus,np.transpose(array_3_helix),cmap=cm_3_helix)
plt.pcolormesh(time,array_residus,np.transpose(array_I_5_helix),cmap=cm_I_5_helix)
plt.pcolormesh(time,array_residus,np.transpose(array_hydro_bonded_turn),cmap=cm_hydro_bonded_turn)
plt.pcolormesh(time,array_residus,np.transpose(array_bend),cmap=cm_bend)



#plt.pcolormesh(time,array_residus,np.transpose(array_loops_and_irregular),cmap=cm_loops_and_irregular)





plt.xlabel('Temps (ns)')
plt.ylabel('Residu')



