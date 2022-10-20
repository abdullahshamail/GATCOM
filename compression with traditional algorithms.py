import numpy as np
import pandas as pd
import time
from functools import reduce
import os
from utils import doInterpolation, loadData, spm3d, doInter, performRDP, rdpOxygen, giveTimesteps, donor_acceptor_dist
from os.path import join, isfile
import sys

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\os_200_1.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogen200k.txt"
    hydrogen_location = r"C:\Users\ashamail\Box\Shared Materials_Abdullah_Hasan\ADBIS PAPER_Everything\EXTENSION\Data\hn_200k_0.txt"
    savePathDistances = r"C:\Users\ashamail\Desktop\GATCOM\data\conventional algorithms\distances"
    savePath = r"C:\Users\ashamail\Desktop\GATCOM\data\conventional algorithms"
    savePathInterpolatedTimes = join(savePath, "times")
    savePathInterpolatedData = join(savePath, "interpolated")
    savePathCompressed = join(savePath, "compressed")
    epsilons = [0.5, 0.75, 1.0]
    currentDrug = 0 #information for RFL1 trajcetory, the drug we used
    start = time.time()
    print("Loading data...")
    oxygenData = loadData(oxygen_location)
    print("Loaded oxygen...")
    nitrogenData = loadData(nitrogen_location)
    print("Loaded nitrogen...")
    hydrogenData = loadData(hydrogen_location)
    print(f"Loaded hydrogen and all data in {time.time() - start} seconds")

    method = "spm"

    for e in epsilons:

        savePathInterpolatedTimes = join(savePath, method, str(e), "times")
        savePathInterpolatedData = join(savePath, method, str(e), "interpolated")
        savePathCompressed = join(savePath, method, str(e), "compressed")
        print(f"Starting epsilon {e}")
        startTime = time.time()
        oxygenTraj = spm3d(oxygenData, e, currentDrug)
        np.savetxt(join(savePathCompressed,"oxygenTrajspm3dEpsilon-"+str(e)+".csv"), oxygenTraj, delimiter=",")

        oxygenTime = oxygenTraj[:,3]
        oxygenInterpolated = doInterpolation(oxygenTraj)
        np.savetxt(join(savePathInterpolatedTimes,"oxygenTrajspm3dinterpolated-"+str(e)+".csv"), oxygenTime, delimiter=",")
        np.savetxt(join(savePathInterpolatedData,"oxygenTrajspm3dinterpolated-"+str(e)+".csv"), oxygenInterpolated, delimiter=",")

        print("Starting N and H")

        for i in range(216): # for the 216 nitrogen trajectories
            nitrogenTrajectory = spm3d(nitrogenData, e, i)
            np.savetxt(join(savePathCompressed,"nitrogenTrajspm3dEpsilon-"+str(e)+"--"+str(i)+".csv"), nitrogenTrajectory, delimiter=",")

            hydrogenTrajectory = spm3d(nitrogenData, e, i)
            np.savetxt(join(savePathCompressed,"hydrogenTrajspm3dEpsilon-"+str(e)+"--"+str(i)+".csv"), hydrogenTrajectory, delimiter=",")

            nitrogenTime = nitrogenTrajectory[:,3]
            hydrogenTime = hydrogenTrajectory[:,3]
            union = reduce(np.union1d, (oxygenTime, nitrogenTime, hydrogenTime))
            nitrogenInterpolated = doInter(nitrogenTrajectory, union)
            hydrogenInterpolated = doInter(hydrogenTrajectory, union)

            nitrogenInterpolated = nitrogenInterpolated[nitrogenInterpolated[:,3].argsort()]
            hydrogenInterpolated = hydrogenInterpolated[hydrogenInterpolated[:,3].argsort()]


            np.savetxt(join(savePathInterpolatedTimes,"nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenTime, delimiter=",")
            np.savetxt(join(savePathInterpolatedData,"nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenInterpolated, delimiter=",")

            np.savetxt(join(savePathInterpolatedTimes,"hydrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), hydrogenTime, delimiter=",")
            np.savetxt(join(savePathInterpolatedData,"hydrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), hydrogenInterpolated, delimiter=",")
            print(f"Done {i} for epsilon {e}")
        print(f"Done epsilon {e} in {time.time() - startTime} seconds")
    
    # print("Compressed and interpolated for SPM.")

    # for e in epsilons:
    #     oxygenTraj = rdpOxygen(oxygenData, e)
    #     oxygenTime = oxygenTraj[:,3]

    #     oxygenInterpolated = doInterpolation(oxygenTraj)
    #     np.savetxt(join(savePathInterpolatedTimes,"oxygenTrajrdpinterpolated-"+str(e)+".csv"), oxygenTime, delimiter=",")
    #     np.savetxt(join(savePathInterpolatedData,"oxygenTrajrdpinterpolated-"+str(e)+".csv"), oxygenInterpolated, delimiter=",")


    #     for i in range(216): # for the 216 nitrogen trajectories
    #         nitrogenTrajectory = performRDP(nitrogenData, e, i)
    #         np.savetxt(join(savePathCompressed,"nitrogenTrajrdpEpsilon-"+str(e)+"--"+str(i)+".csv"), nitrogenTrajectory, delimiter=",")

    #         nitrogenTime = nitrogenTrajectory[:,3]
    #         union = np.union1d(oxygenTime, nitrogenTime)
    #         nitrogenInterpolated = doInter(nitrogenTrajectory, union)
    #         nitrogenInterpolated = nitrogenInterpolated[nitrogenInterpolated[:,3].argsort()]

    #         np.savetxt(join(savePathInterpolatedTimes,"nitrogenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenTime, delimiter=",")
    #         np.savetxt(join(savePathInterpolatedData,"nitrogenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenInterpolated, delimiter=",")

    for e in epsilons:
        savePathInterpolatedTimes = join(savePath, method, str(e), "times")
        savePathInterpolatedData = join(savePath, method, str(e), "interpolated")
        savePathCompressed = join(savePath, method, str(e), "compressed")
        oxygenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajspm3dinterpolated-"+str(e)+".csv"), delimiter=",")

        union = np.array([])
        combinedNitrogens = np.empty((0,4))
        for i in range(216):
            print(i)
            nitrogenData = np.loadtxt(join(savePathInterpolatedData, "nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), delimiter=",")
            combinedNitrogens = np.append(combinedNitrogens, nitrogenData, axis=0)
        
        df = pd.DataFrame(combinedNitrogens)
        times = np.unique(combinedNitrogens[:,3])
        nitrogenDict = {}

        print(len(times))
        for i in times:
            # print("-", i)
            if i%10000 == 0: print("-", i)
            nitrogenDict[i] = np.array(df[df[3] == i][[0,1,2]])
        timesteps = oxygenData[:,3].astype(int)
        oxygenData = pd.DataFrame(oxygenData)
        tempOxygen = np.zeros((200000, 3))
        uniquetimesteps = oxygenData[3].values

        for i in uniquetimesteps:
                xyz = oxygenData[oxygenData[3] == i].values
                i = int(i)
                tempOxygen[i][0] = xyz[0][0]
                tempOxygen[i][1] = xyz[0][1]
                tempOxygen[i][2] = xyz[0][2]
        output, _ = donor_acceptor_dist(tempOxygen, nitrogenDict, list(nitrogenDict.keys()))
        np.savetxt(join(savePathDistances, "distances for spm epsilon-"+str(e), output ,delimiter=","))


    # for e in epsilons:
    #     oxygenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajrdpinterpolated-"+str(e)+".csv"), delimiter=",")

    #     union = np.array([])
    #     combinedNitrogens = np.empty((0,4))
    #     for i in range(216):
    #         nitrogenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), delimiter=",")
    #         combinedNitrogens = np.append(combinedNitrogens, nitrogenData, axis=0)  
    #     df = pd.DataFrame(combinedNitrogens)
    #     times = np.unique(combinedNitrogens[:,3])
    #     nitrogenDict = {}

    #     for i in times:
    #         nitrogenDict[i] = np.array(df[df[3] == i][[0,1,2]])
    #     timesteps = oxygenData[:,3].astype(int)
    #     oxygenData = pd.DataFrame(oxygenData)
    #     tempOxygen = np.zeros((200000, 3))
    #     uniquetimesteps = oxygenData[3].values

    #     for i in uniquetimesteps:
    #             xyz = oxygenData[oxygenData[3] == i].values
    #             i = int(i)
    #             tempOxygen[i][0] = xyz[0][0]
    #             tempOxygen[i][1] = xyz[0][1]
    #             tempOxygen[i][2] = xyz[0][2]
    #     output, _ = donor_acceptor_dist(tempOxygen, nitrogenDict, list(nitrogenDict.keys()))
    #     np.savetxt("distances for rdp epsilon-"+str(e) , output, delimiter=',')

