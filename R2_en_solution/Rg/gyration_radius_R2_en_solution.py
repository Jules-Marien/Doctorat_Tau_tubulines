# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 22:15:14 2022

@author: yacaj
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd




#Import data


Rg_R2_simple_modele17 = np.loadtxt('../R2_simple/modele17/Rg/Rg_R2_simple_modele17.txt')

Rg_R2_simple_modele22 = np.loadtxt('../R2_simple/modele22/Rg/Rg_R2_simple_modele22.txt')

Rg_R2_simple_modele84 = np.loadtxt('../R2_simple/modele84/Rg/Rg_R2_simple_modele84.txt')





Rg_R2_phospho_Ser1351_modele17 = np.loadtxt('../R2_phospho/modele17/Ser1351/Rg/Rg_R2_modele17_phospho_Ser1351.txt')

Rg_R2_phospho_Ser1351_modele22 = np.loadtxt('../R2_phospho/modele22/Ser1351/Rg/Rg_R2_modele22_phospho_Ser1351.txt')

Rg_R2_phospho_Ser1351_modele84 = np.loadtxt('../R2_phospho/modele84/Ser1351/Rg/Rg_R2_modele84_phospho_Ser1351.txt')





Rg_R2_phospho_Ser1355_modele17 = np.loadtxt('../R2_phospho/modele17/Ser1355/Rg/Rg_R2_modele17_phospho_Ser1355.txt')

Rg_R2_phospho_Ser1355_modele22 = np.loadtxt('../R2_phospho/modele22/Ser1355/Rg/Rg_R2_modele22_phospho_Ser1355.txt')

Rg_R2_phospho_Ser1355_modele84 = np.loadtxt('../R2_phospho/modele84/Ser1355/Rg/Rg_R2_modele84_phospho_Ser1355.txt')





Rg_R2_phospho_Ser1359_modele17 = np.loadtxt('../R2_phospho/modele17/Ser1359/Rg/Rg_R2_modele17_phospho_Ser1359.txt')

Rg_R2_phospho_Ser1359_modele22 = np.loadtxt('../R2_phospho/modele22/Ser1359/Rg/Rg_R2_modele22_phospho_Ser1359.txt')

Rg_R2_phospho_Ser1359_modele84 = np.loadtxt('../R2_phospho/modele84/Ser1359/Rg/Rg_R2_modele84_phospho_Ser1359.txt')





Rg_R2_phospho_toutes_Ser_modele17 = np.loadtxt('../R2_phospho/modele17/toutes_Ser/Rg/Rg_R2_modele17_phospho_toutes_Ser.txt')

Rg_R2_phospho_toutes_Ser_modele22 = np.loadtxt('../R2_phospho/modele22/toutes_Ser/Rg/Rg_R2_modele22_phospho_toutes_Ser.txt')

Rg_R2_phospho_toutes_Ser_modele84 = np.loadtxt('../R2_phospho/modele84/toutes_Ser/Rg/Rg_R2_modele84_phospho_toutes_Ser.txt')





#Concatenate data


data_Rg_R2_simple = np.concatenate((Rg_R2_simple_modele17, Rg_R2_simple_modele22, Rg_R2_simple_modele84))

data_Rg_R2_phospho_Ser1351 = np.concatenate((Rg_R2_phospho_Ser1351_modele17, Rg_R2_phospho_Ser1351_modele22, Rg_R2_phospho_Ser1351_modele84))

data_Rg_R2_phospho_Ser1355 = np.concatenate((Rg_R2_phospho_Ser1355_modele17, Rg_R2_phospho_Ser1355_modele22, Rg_R2_phospho_Ser1355_modele84))

data_Rg_R2_phospho_Ser1359 = np.concatenate((Rg_R2_phospho_Ser1359_modele17, Rg_R2_phospho_Ser1359_modele22, Rg_R2_phospho_Ser1359_modele84))

data_Rg_R2_phospho_toutes_Ser = np.concatenate((Rg_R2_phospho_toutes_Ser_modele17, Rg_R2_phospho_toutes_Ser_modele22, Rg_R2_phospho_toutes_Ser_modele84))




"""Plots"""

x_lim_min = 5
x_lim_max = 25

y_lim_min = 0
y_lim_max = 0.2


Rg_R2_simple_pandas = pd.DataFrame(data_Rg_R2_simple, columns = ['Rg'])

plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)

Rg_R2_simple_pandas.Rg.plot.density()


plt.title('Density plot of the gyration radius Rg of R2 in solution (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_Rg_R2_in_solution.png')

plt.show()








Rg_R2_phospho_Ser1351_pandas = pd.DataFrame(data_Rg_R2_phospho_Ser1351, columns = ['Rg'])

plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)

Rg_R2_phospho_Ser1351_pandas.Rg.plot.density()


plt.title('Density plot of the gyration radius Rg of R2 \n phosphorylated on Ser285 in solution (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_Rg_R2_phospho_Ser1351_in_solution.png')

plt.show()









Rg_R2_phospho_Ser1355_pandas = pd.DataFrame(data_Rg_R2_phospho_Ser1355, columns = ['Rg'])

plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)

Rg_R2_phospho_Ser1355_pandas.Rg.plot.density()


plt.title('Density plot of the gyration radius Rg of R2 \n phosphorylated on Ser289 in solution (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_Rg_R2_phospho_Ser1355_in_solution.png')

plt.show()








Rg_R2_phospho_Ser1359_pandas = pd.DataFrame(data_Rg_R2_phospho_Ser1359, columns = ['Rg'])

plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)

Rg_R2_phospho_Ser1359_pandas.Rg.plot.density()


plt.title('Density plot of the gyration radius Rg of R2 \n phosphorylated on Ser293 in solution (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_Rg_R2_phospho_Ser1359_in_solution.png')

plt.show()







Rg_R2_phospho_toutes_Ser_pandas = pd.DataFrame(data_Rg_R2_phospho_toutes_Ser, columns = ['Rg'])

plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)

Rg_R2_phospho_toutes_Ser_pandas.Rg.plot.density()


plt.title('Density plot of the gyration radius Rg of R2 \n phosphorylated on all Serines in solution (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Density_plot_Rg_R2_phospho_toutes_Ser_in_solution.png')

plt.show()











#Comparison




plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)



Rg_R2_phospho_Ser1351_pandas.Rg.plot.density(label="p-Ser285", color="lightgreen")

Rg_R2_phospho_Ser1355_pandas.Rg.plot.density(label="p-Ser289", color="limegreen")

Rg_R2_phospho_Ser1359_pandas.Rg.plot.density(label="p-Ser293", color="forestgreen")

Rg_R2_phospho_toutes_Ser_pandas.Rg.plot.density(label="all phospho", color="red")

Rg_R2_simple_pandas.Rg.plot.density(label="no phospho", color="black")


plt.title('Comparison of density plot of the gyration radius Rg of R2 \n in solution with and without phosphorylation(s) (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Comparison_density_plot_Rg_R2_in_solution.png')

plt.show()






plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)



Rg_R2_phospho_Ser1351_pandas.Rg.plot.density(label="p-Ser285", color="lightgreen")

Rg_R2_phospho_Ser1355_pandas.Rg.plot.density(label="p-Ser289", color="limegreen")

Rg_R2_simple_pandas.Rg.plot.density(label="no phospho", color="black")


plt.title('Comparison of density plot of the gyration radius Rg of R2 \n in solution with and without phosphorylation(s) (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Comparison_beta_folders_density_plot_Rg_R2_in_solution.png')

plt.show()







plt.figure()

plt.xlim(x_lim_min,x_lim_max)
plt.ylim(y_lim_min,y_lim_max)


Rg_R2_phospho_Ser1359_pandas.Rg.plot.density(label="p-Ser293", color="forestgreen")

Rg_R2_phospho_toutes_Ser_pandas.Rg.plot.density(label="all phospho", color="red")

Rg_R2_simple_pandas.Rg.plot.density(label="no phospho", color="black")


plt.title('Comparison of density plot of the gyration radius Rg of R2 \n in solution with and without phosphorylation(s) (600ns)')
plt.xlabel('Gyration radius (Å)')
plt.ylabel('Density')
plt.legend()

plt.savefig('Comparison_alpha_folders_density_plot_Rg_R2_in_solution.png')

plt.show()





