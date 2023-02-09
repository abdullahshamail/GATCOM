import numpy as np
from utils import  perpendicular_distance, donor_acceptor_dist

# Load the representative atom
ether_o =  np.loadtxt(r"C:\Users\mhanowar\Downloads\os atoms 0-800k\os_200k_0.txt")
ether_o = ether_o.reshape(ether_o.shape[0], ether_o.shape[1] // 3 , 3)

# Load the donor atom
donor_nh =  np.loadtxt(r"C:\Users\mhanowar\Downloads\nitrogen 0-800k\nitrogen200k-0.txt")
donor_nh = donor_nh.reshape(donor_nh.shape[0], donor_nh.shape[1] // 3 , 3)

# Load COM trajectory
center_of_mass_traj =  np.loadtxt(r"C:/Users/mhanowar/Downloads/RFL1_center of Mass trajectory_0-200K.txt")

init_frame = 0
final_frame = 200000

total_frame = np.arange(init_frame,final_frame, dtype=int)    


perpendicular_dist = perpendicular_distance(center_of_mass_traj = center_of_mass_traj, total_frame = total_frame)

for perpendicular_threshold in np.arange (-5, 6, 1):
    
    perpendicular_frame = total_frame [perpendicular_dist < perpendicular_threshold]
    perpendicular_frame = np.array(perpendicular_frame,int)
    
    for hbond_threshold in np.arange(3.5,5.5,0.5):
        
        donor_acceptor_timestep, _ = donor_acceptor_dist (atom1_traj = ether_o [:, 0, :], atom2_traj = donor_nh, total_frame = perpendicular_frame, hbond_threshold = hbond_threshold)
        
        np.savetxt(r"C:\Users\mhanowar\Desktop\Varation of Distance Threshold\hbond threshold"+str(hbond_threshold)+ "A with perpendicular threshold"+ str(perpendicular_threshold)+"A.csv", donor_acceptor_timestep, delimiter=",")
    

