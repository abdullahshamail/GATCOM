import numpy as np
import pytraj as pt
import os

if __name__ == "__main__":
    currDir = '/home/abdullah/Downloads'
    filepathDCD = os.path.join(currDir, 'traj_1.dcd')

    filepathPRM = os.path.join(currDir, '4csp_no_sol.prmtop')

    atoms = {
        'nitrogen': ':1-4@%n',
        'oxygen': ':5@%os'
    }

    atomCount = {
        'nitrogen': 216,
        'oxygen': 1
    }

    count = 0
    saveAtom = 'oxygen'
    data = np.empty((0, atomCount[saveAtom], 3), dtype=float)
    while count < 5000:
        traj = pt.iterload(filepathDCD, filepathPRM, frame_slice=(count,count + 500,1))
        test = traj[atoms[saveAtom]].xyz
        data = np.append(data, test, axis=0)
        count += 500
        # print(count)
    # reshape and save
    print(data.shape)
    data_r = data.reshape(data.shape[0], -1)
    np.savetxt(saveAtom+"Traj.txt", data_r)

    # load and reshape to confirm
    l = np.loadtxt(saveAtom+"Traj.txt")
    data_w = l.reshape(l.shape[0], l.shape[1] // 3 , 3)
    print(data_w.shape)