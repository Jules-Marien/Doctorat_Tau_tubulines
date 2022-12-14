import matplotlib.pyplot as plt 
import numpy as np


"""Import data"""



#alpha1_beta1_avec_tau

RMSD_alpha1_beta1_avec_tau_modele17 = np.genfromtxt('../data/alpha1_beta1_avec_tau/RSMDs_alpha1_beta1_avec_tau_modele17.txt')[:,6]
RMSD_alpha1_beta1_avec_tau_modele22 = np.genfromtxt('../data/alpha1_beta1_avec_tau/RSMDs_alpha1_beta1_avec_tau_modele22.txt')[:,6]
RMSD_alpha1_beta1_avec_tau_modele84 = np.genfromtxt('../data/alpha1_beta1_avec_tau/RSMDs_alpha1_beta1_avec_tau_modele84.txt')[:,6]



#alpha1_beta1_sans_tau

RMSD_alpha1_beta1_sans_tau_modele17 = np.genfromtxt('../data/alpha1_beta1_sans_tau/RSMDs_alpha1_beta1_sans_tau_modele17.txt')[:,6]
RMSD_alpha1_beta1_sans_tau_modele22 = np.genfromtxt('../data/alpha1_beta1_sans_tau/RSMDs_alpha1_beta1_sans_tau_modele22.txt')[:,6]
RMSD_alpha1_beta1_sans_tau_modele84 = np.genfromtxt('../data/alpha1_beta1_sans_tau/RSMDs_alpha1_beta1_sans_tau_modele84.txt')[:,6]



#a1_b3_avec_tau

RMSD_a1_b3_avec_tau_modele16 = np.genfromtxt('../data/a1_b3_avec_tau/RSMDs_a1_b3_avec_tau_modele16.txt')[:,6]
RMSD_a1_b3_avec_tau_modele38 = np.genfromtxt('../data/a1_b3_avec_tau/RSMDs_a1_b3_avec_tau_modele38.txt')[:,6]
RMSD_a1_b3_avec_tau_modele65 = np.genfromtxt('../data/a1_b3_avec_tau/RSMDs_a1_b3_avec_tau_modele65.txt')[:,6]



#a1_b3_sans_tau

RMSD_a1_b3_sans_tau_modele16 = np.genfromtxt('../data/a1_b3_sans_tau/RSMDs_a1_b3_sans_tau_modele16.txt')[:,6]
RMSD_a1_b3_sans_tau_modele38 = np.genfromtxt('../data/a1_b3_sans_tau/RSMDs_a1_b3_sans_tau_modele38.txt')[:,6]
RMSD_a1_b3_sans_tau_modele65 = np.genfromtxt('../data/a1_b3_sans_tau/RSMDs_a1_b3_sans_tau_modele65.txt')[:,6]



#sans_CTT

RMSD_sans_CTT_modele17 = np.genfromtxt('../data/modeles_sans_CTT/RSMDs_modele17_sans_CTT.txt')[:,6]
RMSD_sans_CTT_modele22 = np.genfromtxt('../data/modeles_sans_CTT/RSMDs_modele22_sans_CTT.txt')[:,6]
RMSD_sans_CTT_modele84 = np.genfromtxt('../data/modeles_sans_CTT/RSMDs_modele84_sans_CTT.txt')[:,6]



time = np.arange(0, len(RMSD_alpha1_beta1_avec_tau_modele17)) /10


y_min = 0
y_max = 35









plt.figure(figsize=(8,6))

plt.title('RMSD of the CTTs \n for trajectories with $\\bf{??I}$ $\\bf{CTTs}$ and $\\bf{with}$ R2')

plt.xlim(0,200)
plt.ylim(y_min,y_max)

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (??)')


plt.plot(time, RMSD_alpha1_beta1_avec_tau_modele17, color='skyblue', label = 'model 17' )

plt.plot(time, RMSD_alpha1_beta1_avec_tau_modele22, color='dodgerblue', label = 'model 22' ) 

plt.plot(time, RMSD_alpha1_beta1_avec_tau_modele84, color='blue', label = 'model 84' ) 


plt.legend()

plt.savefig('RMSD_CTTs_beta1_avec_tau.png')








plt.figure(figsize=(8,6))

plt.title('RMSD of the CTTs \n for trajectories with $\\bf{??I}$ $\\bf{CTTs}$ and $\\bf{without}$ R2')

plt.xlim(0,200)
plt.ylim(y_min,y_max)

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (??)')


plt.plot(time, RMSD_alpha1_beta1_sans_tau_modele17, color='skyblue', label = 'model 17' )

plt.plot(time, RMSD_alpha1_beta1_sans_tau_modele22, color='dodgerblue', label = 'model 22' ) 

plt.plot(time, RMSD_alpha1_beta1_sans_tau_modele84, color='blue', label = 'model 84' ) 


plt.legend()

plt.savefig('RMSD_CTTs_beta1_sans_tau.png')













plt.figure(figsize=(8,6))

plt.title('RMSD of the CTTs \n for trajectories with $\\bf{??III}$ $\\bf{CTTs}$ and $\\bf{with}$ R2')

plt.xlim(0,200)
plt.ylim(y_min,y_max)

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (??)')


plt.plot(time, RMSD_a1_b3_avec_tau_modele16, color='orange', label = 'model 16' )

plt.plot(time, RMSD_a1_b3_avec_tau_modele38, color='red', label = 'model 38' ) 

plt.plot(time, RMSD_a1_b3_avec_tau_modele65, color='firebrick', label = 'model 65' ) 


plt.legend()

plt.savefig('RMSD_CTTs_b3_avec_tau.png')










plt.figure(figsize=(8,6))

plt.title('RMSD of the CTTs \n for trajectories with $\\bf{??III}$ $\\bf{CTTs}$ and $\\bf{without}$ R2')

plt.xlim(0,200)
plt.ylim(y_min,y_max)

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (??)')


plt.plot(time, RMSD_a1_b3_sans_tau_modele16, color='orange', label = 'model 16' )

plt.plot(time, RMSD_a1_b3_sans_tau_modele38, color='red', label = 'model 38' ) 

plt.plot(time, RMSD_a1_b3_sans_tau_modele65, color='firebrick', label = 'model 65' ) 


plt.legend()

plt.savefig('RMSD_CTTs_b3_sans_tau.png')










"""


plt.figure(figsize=(8,6))

plt.title('RMSD of the CTTs \n for trajectories $\\bf{without}$ $\\bf{CTTs}$ and $\\bf{with}$ R2')

plt.xlim(0,200)
plt.ylim(y_min,y_max)

plt.xlabel('Time (ns)')
plt.ylabel('RMSD (??)')


plt.plot(time, RMSD_sans_CTT_modele17, color='lightgrey', label = 'model 17' )

plt.plot(time, RMSD_sans_CTT_modele22, color='grey', label = 'model 22' ) 

plt.plot(time, RMSD_sans_CTT_modele84, color='black', label = 'model 84' ) 


plt.legend()

plt.savefig('RMSD_CTTs_sans_CTT.png')



"""




"""Mean values"""


with open('Valeurs_moyennes_RMSD_CTTs.txt','w') as outputfile:
     outputfile.write("Mean values of RMSD on backbone for the CTTs" + "\n")

     outputfile.write("\n")
    
     outputfile.write("alpha1/beta1 with R2" + "\n")

     outputfile.write("model 17 = " + str(np.mean(RMSD_alpha1_beta1_avec_tau_modele17)) + "\n")
     outputfile.write("model 22 = " + str(np.mean(RMSD_alpha1_beta1_avec_tau_modele22)) + "\n")
     outputfile.write("model 84 = " + str(np.mean(RMSD_alpha1_beta1_avec_tau_modele84)) + "\n")
     outputfile.write("all 3 combined = " + str(np.mean(np.concatenate((RMSD_alpha1_beta1_avec_tau_modele17, RMSD_alpha1_beta1_avec_tau_modele22, RMSD_alpha1_beta1_avec_tau_modele84)))) + "\n")
     


     outputfile.write("\n")
    
     outputfile.write("alpha1/beta1 without R2" + "\n")

     outputfile.write("model 17 = " + str(np.mean(RMSD_alpha1_beta1_sans_tau_modele17)) + "\n")
     outputfile.write("model 22 = " + str(np.mean(RMSD_alpha1_beta1_sans_tau_modele22)) + "\n")
     outputfile.write("model 84 = " + str(np.mean(RMSD_alpha1_beta1_sans_tau_modele84)) + "\n")
     outputfile.write("all 3 combined = " + str(np.mean(np.concatenate((RMSD_alpha1_beta1_sans_tau_modele17, RMSD_alpha1_beta1_sans_tau_modele22, RMSD_alpha1_beta1_sans_tau_modele84)))) + "\n")





     outputfile.write("\n")
    
     outputfile.write("a1/b3 with R2" + "\n")

     outputfile.write("model 16 = " + str(np.mean(RMSD_a1_b3_avec_tau_modele16)) + "\n")
     outputfile.write("model 38 = " + str(np.mean(RMSD_a1_b3_avec_tau_modele38)) + "\n")
     outputfile.write("model 65 = " + str(np.mean(RMSD_a1_b3_avec_tau_modele65)) + "\n")
     outputfile.write("all 3 combined = " + str(np.mean(np.concatenate((RMSD_a1_b3_avec_tau_modele16, RMSD_a1_b3_avec_tau_modele38, RMSD_a1_b3_avec_tau_modele65)))) + "\n")
     
     

     outputfile.write("\n")
    
     outputfile.write("a1/b3 without R2" + "\n")

     outputfile.write("model 16 = " + str(np.mean(RMSD_a1_b3_sans_tau_modele16)) + "\n")
     outputfile.write("model 38 = " + str(np.mean(RMSD_a1_b3_sans_tau_modele38)) + "\n")
     outputfile.write("model 65 = " + str(np.mean(RMSD_a1_b3_sans_tau_modele65)) + "\n")
     outputfile.write("all 3 combined = " + str(np.mean(np.concatenate((RMSD_a1_b3_sans_tau_modele16, RMSD_a1_b3_sans_tau_modele38, RMSD_a1_b3_sans_tau_modele65)))) + "\n")
     


"""


     outputfile.write("\n")
    
     outputfile.write("no CTTs" + "\n")

     outputfile.write("model 17 = " + str(np.mean(RMSD_sans_CTT_modele17)) + "\n")
     outputfile.write("model 22 = " + str(np.mean(RMSD_sans_CTT_modele22)) + "\n")
     outputfile.write("model 84 = " + str(np.mean(RMSD_sans_CTT_modele84)) + "\n")
     outputfile.write("all 3 combined = " + str(np.mean(np.concatenate((RMSD_sans_CTT_modele17, RMSD_sans_CTT_modele22, RMSD_sans_CTT_modele84)))) + "\n")


"""


