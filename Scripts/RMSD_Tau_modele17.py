#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:17:57 2022

@author: marien
"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 


from MDAnalysis.analysis import rms

import pandas as pd

from numba import jit


#largely inspired from https://userguide.mdanalysis.org/1.1.1/examples/analysis/distances_and_contacts/contacts_within_cutoff.html



modele = '17'

data = np.genfromtxt('RMSD_Tau_modele'+modele+'.dat')


plt.figure()
plt.title('RMSD du squelette du fragment de Tau pour le modele '+modele+' \n par rapport à la structure de départ ')

plt.xlabel('Temps (ns)')
plt.ylabel('RMSD (Å)')

plt.plot(data[:,0]/10, data[:,1])

plt.savefig('RMSD_Tau_modele'+modele+'.png')

plt.show()