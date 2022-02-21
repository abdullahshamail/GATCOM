import numpy as np
import plotly.graph_objects as go

def doInterpolation(data):
    for i in range(data.shape[0] - 1):
        data = interpolation(data[i,:], data[i+1,:], data)
    return data[data[:,3].argsort()]

def interpolation(p0, p1, data):
    x0, y0, z0, t0 = p0
    x1, y1, z1, t1 = p1
    unitVector = np.array([x1-x0, y1-y0, z1-z0])
    magnitudeUnit = np.linalg.norm(unitVector)
#     x = [x0]
#     y = [y0]
#     z = [z0]
    howMany = int(t1 - t0)
#     print(howMany)
    if (howMany > 1):
        for i in range(1, howMany):
            tempx = x0 + (magnitudeUnit*i/howMany)*(x1-x0)/(magnitudeUnit)
            tempy = y0 + (magnitudeUnit*i/howMany)*(y1-y0)/(magnitudeUnit)
            tempz = z0 + (magnitudeUnit*i/howMany)*(z1-z0)/(magnitudeUnit)
            tempt = t0 + i
            data = np.append(data, np.array([[tempx, tempy, tempz, tempt]]), axis=0)
#             print("Hello")
#     x.append(x1)
#     y.append(y1)
#     z.append(z1)
    return data

def plotCoord(x, y, z):
    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z,
                                   mode='markers')])
    fig.show()

def interpole(p0, p1, i):
    x0, y0, z0, t0 = p0
    x1, y1, z1, t1 = p1
    unitVector = np.array([x1-x0, y1-y0, z1-z0])
    magnitudeUnit = np.linalg.norm(unitVector)

    howMany = int(t1 - t0)
    tempx = x0 + (magnitudeUnit*(i - t0)/howMany)*(x1-x0)/(magnitudeUnit)
    tempy = y0 + (magnitudeUnit*(i - t0)/howMany)*(y1-y0)/(magnitudeUnit)
    tempz = z0 + (magnitudeUnit*(i - t0)/howMany)*(z1-z0)/(magnitudeUnit)
    tempt = i

    return np.array([[tempx, tempy, tempz, tempt]])

def doInter(data, timesteps):
    j = 0
    for i, t in enumerate(timesteps):
        if (t != data[j][3]):
                temp = interpole(data[j-1], data[j], t)
                data = np.append(data, temp, axis = 0)
        if(t == data[j][3]):
            j = j+1
    return data

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

    x_coord_diff = np.array([(72-abs(item)) if abs(item)>36 else abs(item) for item in (X1-X2)])
    
    y_coord_diff = np.array([(72-abs(item)) if abs(item)>36 else abs(item) for item in (Y1-Y2)])
        
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

def magnitude(vector):
    return np.linalg.norm(vector)

def numerator(currentPoint, basePointA, basePointB):
    xp, yp, zp = currentPoint
    xa, ya, za = basePointA
    xb, yb, zb = basePointB
    
    return magnitude(
            np.cross(
                    [xp - xa, yp - ya, zp - za], 
                    [xb - xa, yb - ya, zb - za]
                )
            )

def denominator(basePointA, basePointB):
    xa, ya, za = basePointA
    xb, yb, zb = basePointB
    
    return magnitude([xb - xa, yb - ya, zb - za])

def giveDistance(point, a, b):
    return numerator(point, a, b) / denominator(a, b)


def spm3d(data, epsilon, location):
    drugTraj = data[:,location,:]

    a = drugTraj[0]
    b = drugTraj[-1]

    newTrajectory = np.array([[a[0], a[1], a[2], 0]])

    for i in range(1, len(drugTraj) - 1):
        p = drugTraj[i]
        distance = giveDistance(p, a, b)
        if distance > epsilon:
            newTrajectory = np.append(newTrajectory, np.array([[p[0], p[1], p[2], i]]), axis=0)
            a = p
    newTrajectory = np.append(newTrajectory, np.array([[b[0], b[1], b[2], len(drugTraj) - 1]]), axis=0)

    return newTrajectory