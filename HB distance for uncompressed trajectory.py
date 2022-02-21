import numpy as np

def split_xyz (coords):
    """
    Parameters
    ----------
    coords : 3D Coordinates

    Returns
    -------
    Coordinates on each axis
 
    """
    # Calculate the x coordinates
    X1 = coords[:,0]

    # Calculate the y coordinates
    Y1 = coords[:,1]

    # Calculate the z coordinates
    Z1 = coords[:,2]
    
    return X1, Y1, Z1

def distance_periodicity (X1,Y1,Z1,X2,Y2,Z2):
    
    """
    Apply periodic boundaries
    
    Periodic Boundaries in X and Y direction is 72 A.
    Z axis does not have any periodic boundary
    
    So, the periodicity conditions

    Distance in X axis = if |x2 − x1|> 36:
                                72 −|x2 − x1|
                            
                         else:
                                x2 − x1
    
    Distance in Y axis = same as X
    
    Euclidean Shortest Distance = √(〖"(Distance in X axis)" 〗^2+〖"(Distance in Y axis)" 〗^2+〖"(Distance in Z axis)" 〗^2 )

    Parameters
    ----------
    X1 : Numpy 1D array
        x axis coordinates of atom 1 trajectory
    Y1 : Numpy 1D array
        y axis coordinates of atom 1 trajectory 
    Z1 : Numpy 1D array
        z axis coordinates of atom 1 trajectory 
    X2 : Numpy 1D array
        x axis coordinates of atom 2 trajectory 
    Y2 : Numpy 1D array
        y axis coordinates of atom 2 trajectory 
    Z2 : Numpy 1D array
        z axis coordinates of atom 2 trajectory 

    Returns
    -------
    distance: Numpy 1D array
            Euclidean distance of (X1,Y1,Z1) and (X2,Y2,Z2)

    """

    x_coord_diff = np.array([(72-item) if abs(item)>36 else item for item in (X1-X2)])
    
    y_coord_diff = np.array([(72-item) if abs(item)>36 else item for item in (Y1-Y2)])
        
    z_coord_diff = np.array((Z1 - Z2))

    euclidean_distance = np.sqrt(np.square(x_coord_diff) + np.square(y_coord_diff) + np.square(z_coord_diff))
    
    return euclidean_distance, x_coord_diff, y_coord_diff, z_coord_diff

def donor_acceptor_dist (atom1_traj, atom2_traj, total_frame ):
    
    """
        Calculate the time frames and donor atoms where donor acceptor distance is within 3.5A
        
    Parameters
    ----------
    atom1_traj : 3D Coordinates of drug acceptor
    
    atom2_traj : 3D Coordinates of Polymer Donor
    
    total_frame : Timesteps from either perpendicular distance or interpolation

    Returns
    -------
    
    donor_acceptor_timestep : Timesteps where exists an instance of donor-acceptance which are within 3.5 A
    
    donor_atom_index : Index of donor (nh) atoms in every timestep that is within 3.5 A of acceptor
                        If there is multiple donor atoms that met the criteria, the closest donor atom is returned.
    """
    
    donor_acceptor_timestep = np.empty((0,1), dtype=int, order='C')
    
    donor_atom_index = np.empty((0,1), dtype=int, order='C')
            
    for frame_num in total_frame:
        
        atom1_traj_xyz = atom1_traj[frame_num]
        
        atom2_traj_xyz = atom2_traj[frame_num]
    
        atom2_traj_x, atom2_traj_y, atom2_traj_z = split_xyz (atom2_traj_xyz)
    
        donor_acceptor_dist,x_coord_diff, _, _ = distance_periodicity (atom1_traj_xyz[0], atom1_traj_xyz[1],atom1_traj_xyz[2], atom2_traj_x, atom2_traj_y, atom2_traj_z)

    
        if min(donor_acceptor_dist) <= 3.5 :
            
            index = np.argmin(donor_acceptor_dist) 
        
            donor_acceptor_timestep = np.append(donor_acceptor_timestep, frame_num)
            
            donor_atom_index = np.append(donor_atom_index, index)
        
    return donor_acceptor_timestep, donor_atom_index

def loadData(link):
    l = np.loadtxt(link)
    data = l.reshape(l.shape[0], l.shape[1] // 3 , 3)
    
    return data

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\oxygenTraj.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogenTraj.txt"

    oxygen = loadData(oxygen_location)
    nitrogen = loadData(nitrogen_location)

    a, b = donor_acceptor_dist(oxygen[:,0,:], nitrogen, np.arange(200000))

    np.savetxt(r"original trajectory without compression.csv", a, delimiter=",")