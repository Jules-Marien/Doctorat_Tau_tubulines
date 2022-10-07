# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:27:13 2022

@author: yacaj
"""

import numpy as np
import matplotlib.pyplot as plt


modele = '17'
serine = '1351'


data = np.genfromtxt('Hbonds_Ser'+serine+'_modele'+modele+'.dat')

frames = data[:2000,0]

Hbonds = data[:2000,1]



time = frames /10





plt.figure(figsize=(12,10))

police = 20

plt.title('Nombre de liaisons hydrogène en fonction du temps \n de la serine '+serine+ ' du modele '+modele+ ' avec tous les autres résidus' , fontsize = 15)

plt.xlabel('Temps (ns)', fontsize = police)
plt.ylabel('Nombre de liaisons hydrogène', fontsize = police)

plt.scatter(time, Hbonds, s=1)

plt.savefig('Hbonds_Ser'+serine+'_modele'+modele+'.png')

plt.show()
