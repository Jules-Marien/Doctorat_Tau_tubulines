import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

import numpy as np

modele = '84'








"""RMSD of the whole system (alignment on the whole system)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "protein"

rmsd_selection = "protein"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

#aligner = align.AlignTraj(u, ref, select=align_selection,
 #                          filename='aligned_to_first_frame.dcd').run()
#u = mda.Universe('../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', 'aligned_to_first_frame.dcd')






RMSD_whole_system = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_whole_system[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1
    













"""RMSD of all cores (alignment on all cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_all_cores = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_all_cores[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1















"""RMSD of CTT A (alignment on CTT A)"""

"""
u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 427:444"

rmsd_selection = "resid 427:444"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')



"""


RMSD_CTT_A = np.ones((len_traj))*(-1)












"""RMSD of CTT B (alignment on CTT B)"""

"""
u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 883:895"

rmsd_selection = "resid 883:895"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')


"""



RMSD_CTT_B = np.ones((len_traj))*(-1)












"""RMSD of CTT C (alignment on CTT C)"""

"""
u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1322:1339"

rmsd_selection = "resid 1322:1339"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')

"""




RMSD_CTT_C = np.ones((len_traj))*(-1)


















"""RMSD of all CTTs  (alignment on all CTTs)"""

"""
u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 427:444 or resid 883:895 or resid 1322:1339"

rmsd_selection = "resid 427:444 or resid 883:895 or resid 1322:1339"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')


"""



RMSD_all_CTTs = np.ones((len_traj))*(-1)


















"""RMSD of R2 (alignment on R2)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_sans_CTT.pdb', '../aligned_traj_reduite_modele'+modele+'_sans_CTT.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1340:1366"

rmsd_selection = "resid 1340:1366"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True).run()



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_R2 = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_R2[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1







#%%


"""Creation of the final .txt file"""




RMSD = np.array([np.arange(0,len_traj),
                  RMSD_whole_system,
                  RMSD_all_cores,
                  RMSD_CTT_A,
                  RMSD_CTT_B,
                  RMSD_CTT_C,
                  RMSD_all_CTTs,
                  RMSD_R2])


#%%





import pandas as pd


column_names = ["Frame", 
                "RMSD of the whole system (alignment on whole system)", 
                "RMSD of all cores (alignement on cores)",
                "RMSD of CTT A (beta chain, alignment on CTT A)",
                "RMSD of CTT B (beta chain, alignment on CTT B)",
                "RMSD of CTT C (beta chain, alignment on CTT C)",
                "RMSD of all CTTs (alignment on all CTTs)",
                "RMSD of R2 (alignment on R2)"]

df = pd.DataFrame(data= RMSD.T,
                  index = np.arange(0,len_traj),
                  columns=column_names)



RMSD_txt = df.to_numpy()

np.savetxt("RSMDs_modele"+modele+"_sans_CTT.txt", RMSD_txt)



"""Write data to files"""





with open("Header_RMSD_file.txt", "w") as header_RMSD:
    for column in column_names:
        header_RMSD.write(column + '\n')


df.to_csv(r'RMSDs_modele'+modele+'_sans_CTT.csv', index=False)






