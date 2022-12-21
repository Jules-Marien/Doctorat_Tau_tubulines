import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

import numpy as np

modele = '84'








"""RMSD of the whole system (alignment on the cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "protein"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_whole_system = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_whole_system[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1
    






"""RMSD of core A (alignment on core A)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426"

rmsd_selection = "resid 1:426"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_core_A = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_core_A[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1










"""RMSD of core B (alignment on core B)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 445:882"

rmsd_selection = "resid 445:882"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_core_B = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_core_B[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1














"""RMSD of core C (alignment on core C)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')


len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 896:1321"

rmsd_selection = "resid 896:1321"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_core_C = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_core_C[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1











"""RMSD of all cores (alignment on all cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



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















"""RMSD of CTT A (alignment on all cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "resid 427:444"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_CTT_A = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_CTT_A[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1












"""RMSD of CTT B (alignment on all cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "resid 883:895"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_CTT_B = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_CTT_B[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1












"""RMSD of CTT C (alignment on all cores)"""


u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_sans_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_sans_tau.dcd')


len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1:426 or resid 445:882 or resid 896:1321"

rmsd_selection = "resid 1322:1339"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')






RMSD_CTT_C = np.zeros((len_traj))


ref_backbone = ref.select_atoms(rmsd_selection)


frame = 0
for ts in u.trajectory:
    
    u_backbone = u.select_atoms(rmsd_selection)

    RMSD_CTT_C[frame] = rms.rmsd(u_backbone.positions, ref_backbone.positions, superposition=False)

    frame = frame +1













"""RMSD of R2 (alignment on all cores)"""

"""
u = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')

ref = mda.Universe( '../PDB_protein_modele'+modele+'_beta1_avec_tau.pdb', '../aligned_traj_reduite_modele'+modele+'_beta1_avec_tau.dcd')



len_traj = len(u.trajectory)



ref.trajectory[0] #set ref to first frame 



align_selection = "resid 1340:1366"

rmsd_selection = "resid 1340:1366"

align_selection = "backbone and "  + align_selection

rmsd_selection = "backbone and " + rmsd_selection



#Align based on align_selection

aligner = align.AlignTraj(u,ref, select=align_selection, in_memory=True)



#If you don’t have enough memory to do that, write the trajectory out to a file and reload it into MDAnalysis (uncomment the cell below). Otherwise, you don’t have to run it.

# aligner = align.AlignTraj(u, ref, select=align_selection,
#                           filename='aligned_to_first_frame.dcd').run()
# u = mda.Universe(PSF, 'aligned_to_first_frame.dcd')


"""



RMSD_R2 = np.ones((len_traj))*(-1)










#%%


"""Creation of the final .txt file"""




RMSD = np.array([np.arange(0,len_traj),
                  RMSD_whole_system,
                  RMSD_core_A,
                  RMSD_core_B,
                  RMSD_core_C,
                  RMSD_all_cores,
                  RMSD_CTT_A,
                  RMSD_CTT_B,
                  RMSD_CTT_C,
                  RMSD_R2])


#%%





import pandas as pd


column_names = ["Frame", 
                "RMSD of the whole system (alignment on cores)", 
                "RMSD of tubulin core A (alignement on tubulin core A)", 
                "RMSD of tubulin core B (alignement on tubulin core B)",
                "RMSD of tubulin core C (alignement on tubulin core C)",
                "RMSD of all cores (alignement on cores)",
                "RMSD of CTT A (beta chain, alignment on core A)",
                "RMSD of CTT B (beta chain, alignment on core B)",
                "RMSD of CTT C (beta chain, alignment on core C)",
                "RMSD of R2 (alignment on R2)"]

df = pd.DataFrame(data= RMSD.T,
                  index = np.arange(0,len_traj),
                  columns=column_names)



RMSD_txt = df.to_numpy()

np.savetxt("RSMDs_alpha1_beta1_sans_tau_modele"+modele+".txt", RMSD_txt)



"""Write data to files"""





with open("Header_RMSD_file.txt", "w") as header_RMSD:
    for column in column_names:
        header_RMSD.write(column + '\n')


df.to_csv(r'RMSDs_alpha1_beta1_sans_tau_modele'+modele+'.csv', index=False)






