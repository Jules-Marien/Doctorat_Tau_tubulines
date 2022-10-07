#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:08:38 2022

@author: marien
"""

import numpy as np

modele = '17'


Comment = 'COMs of modele '+modele +': chain A is C, chain B is N, chain C is O'



atom_type_chainA ='C'
atom_type_chainB ='N'
atom_type_chainC= 'O'


COM_chainA = np.genfromtxt('../COM_rho_phi/COM_chainA_modele'+modele+'.txt')
COM_chainB = np.genfromtxt('../COM_rho_phi/COM_chainB_modele'+modele+'.txt')
COM_chainC = np.genfromtxt('../COM_rho_phi/COM_chainC_modele'+modele+'.txt')

count_A = len(COM_chainA[:,0])
count_B = len(COM_chainB[:,0])
count_C = len(COM_chainC[:,0])

total_count = count_A + count_B + count_C




xyz_file = np.empty((total_count +2), dtype=object)

xyz_file[0] = str(total_count) + "\n"

xyz_file[1] = Comment + "\n"



for i in range(count_A):
    xyz_file[2+i] = atom_type_chainA +" "+ str(COM_chainA[i,0]) +" " + str(COM_chainA[i,1]) +" " + str(COM_chainA[i,2]) + "\n"
    
for j in range(count_B):
    xyz_file[2+count_A+j] = atom_type_chainB +" "+ str(COM_chainB[j,0]) +" " + str(COM_chainB[j,1]) +" " + str(COM_chainB[j,2]) + "\n"

for k in range(count_C):
    xyz_file[2+count_A+count_B+k] = atom_type_chainC +" "+ str(COM_chainC[k,0]) +" " + str(COM_chainC[k,1]) +" " + str(COM_chainC[k,2]) + "\n"


with open('COM_modele'+modele+'.xyz','w') as outputfile:    
    for i in range(total_count +2):
            
        outputfile.write(xyz_file[i])
        
        print(i)

"""
count = len(data[:,0])

Comment = "COM of modele " +modele

new_array = np.empty((len(data[:,0])+2), dtype=object)

new_array[0] = str(count) +"\n"
new_array[1] = Comment +"\n"

for i in range(count):
    new_array[2+i] = atom_type + " "+ str(data[i,0]) +" " + str(data[i,1]) + " " + str(data[i,2]) + "\n"
    
    
with open('COM_modele17.xyz','w') as outputfile:   
    
    for i in range(len(new_array)):
            
        outputfile.write(new_array[i])
        
        print(i)
    
    """
