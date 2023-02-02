from os import listdir
from os.path import isfile, join
import numpy as np
from utils import loadData, donor_acceptor_dist

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\os_200_1.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogen200k.txt"

    oxygen = loadData(oxygen_location)
    nitrogen = loadData(nitrogen_location)

    whichCOM = ["OS atom as COM", "True COM"]

    timeStepsFolder = r"C:\Users\ashamail\OneDrive - Iowa State University\Molecular Simulation Data\Text Data Files\Timesteps using Generalization\RFL1 Timesteps  with Generalization" #these are the timesteps where d_delta < 3.5 
    saveFolder = r"GATCOM"
    for typ in whichCOM:
        path = join(timeStepsFolder, typ)
        files = [f for f in listdir(path) if isfile(join(path, f))]
        for f in files:
            if f.endswith('.csv'):
                timesteps = np.loadtxt(join(path, f)).astype(int)
                a, b = donor_acceptor_dist(oxygen[:,0,:], nitrogen, timesteps)
                # np.savetxt(join(saveFolder, typ, f), a, delimiter=",")
                print(typ, f, a.shape)