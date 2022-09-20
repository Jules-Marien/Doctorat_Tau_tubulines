#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 18:36:55 2022

@author: marien
"""

import mdtraj as md
import numpy as np
import time



traj = md.load('/ibpc/tethys/marien/Documents/Analyse_Tau_seule/namd_analyse/traj_complete.dcd',top = '/ibpc/tethys/marien/Documents/Analyse_Tau_seule/namd_analyse/step3_input.psf')


nbre_of_frames = 1000


#We generate the dssp data with all structuration possible
dssp_brut = md.compute_dssp(traj[0:nbre_of_frames], simplified=False)


#We only keep the data concerning the protein, meaning the 27 first results (bc 27 amino acids)

len_protein = 27


dssp = np.array(dssp_brut[:,:len_protein])

dssp_treated = np.copy(dssp)

#We complete the table in order to obtain an easier output to handle with numpy
for i in range(len_protein):
    for j in range(nbre_of_frames):
        if dssp_treated[j,i] == ' ':
            dssp_treated[j,i] = 0

#We obtain an array with 0s instead of empty spaces -> 0s represent a lack of structure of any kind


np.savetxt('DSSP_traj.txt',dssp_treated,fmt='%s')
