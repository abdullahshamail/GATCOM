import numpy as np
import time
from utils import donor_acceptor_dist

# Load the representative atom
ether_o =  np.loadtxt(r"C:\Users\mhanowar\Downloads\os atoms 0-800k\os_200k_0.txt")
ether_o = ether_o.reshape(ether_o.shape[0], ether_o.shape[1] // 3 , 3)

# Load the donor atom
donor_nh =  np.loadtxt(r"C:\Users\mhanowar\Downloads\nitrogen 0-800k\nitrogen200k-0.txt")
donor_nh = donor_nh.reshape(donor_nh.shape[0], donor_nh.shape[1] // 3 , 3)

# Load the time frame after classical compression
spm3D_compressed_traj = np.loadtxt(r"C:\Users\mhanowar\Box\Iowa State Research\Shared Materials_Abdullah_Hasan\rdp and spm compressed ether oxygen and true com\etherOspm0.5.csv",delimiter=",")
spm3D_compressed_frame = spm3D_compressed_traj[:,3].astype(int)
epsilon = 0.5

start = time.time()

# compute the HB instances
donor_acceptor_timestep, _ = donor_acceptor_dist (atom1_traj = ether_o [:, 0, :], atom2_traj = donor_nh, total_frame = spm3D_compressed_frame, hbond_threshold = 3.5 + 2* epsilon)

print (time.time() - start)    
