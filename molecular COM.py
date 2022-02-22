import numpy as np
from utils import  split_xyz, calc_center_of_mass

# Loading the drug trajectory

drug_atoms =  np.loadtxt(r"C:\Users\mhanowar\Downloads\rfl1_200k_0.txt")
drug_atoms = drug_atoms.reshape(drug_atoms.shape[0], drug_atoms.shape[1] // 3 , 3)

atomic_mass_list = {"O" : 15.9994, "C" : 12.0107, "H" : 1.00794}

drug_atoms_number = {"O" : 2, "C" : 15, "H" : 12}

# As coordinates of the drug atom are written sequentially, we can save O, C and H atomic masses in an array

atomic_mass = np.concatenate((np.repeat(atomic_mass_list.get('O'), drug_atoms_number.get('O')), 
                              np.repeat (atomic_mass_list.get('C'), drug_atoms_number.get('C')), 
                              np.repeat (atomic_mass_list.get('H'), drug_atoms_number.get('H'))), axis = 0)


# Total time frame = 200,000

init_frame = 0
final_frame = 200000
total_frame = np.arange(init_frame,final_frame, dtype=int)

center_of_mass_traj = np.empty((0,3), dtype = float, order='C')

for frame_num in total_frame:
    
    drug_atoms_xyz = drug_atoms[frame_num]
    xcoords, ycoords, zcoords = split_xyz (drug_atoms_xyz)
    com_x, com_y, com_z = calc_center_of_mass (atomic_mass, xcoords, ycoords, zcoords)
    center_of_mass_traj = np.append(center_of_mass_traj, [[com_x, com_y, com_z]] , axis = 0)
    
np.savetxt(r"C:\Users\mhanowar\Desktop\RFL1_Center of Mass trajectory_0-200K.txt", center_of_mass_traj)
