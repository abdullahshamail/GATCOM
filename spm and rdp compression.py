import numpy as np
import os
from utils import doInterpolation, loadData, spm3d, doInter, performRDP, rdpOxygen, giveTimesteps
from os.path import join, isfile

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\oxygenTraj.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogenTraj.txt"
    savePath = r"C:\Users\ashamail\Desktop\drug-polymers\data"
    savePathInterpolatedTimes = join(savePath, "times")
    savePathInterpolatedData = join(savePath, "interpolated")
    savePathCompressed = join(savePath, "compressed")
    epsilons = [0.5, 0.75, 1.0]
    currentDrug = 0 #information for RFL1 trajcetory, the drug we used

    oxygenData = loadData(oxygen_location)
    nitrogenData = loadData(nitrogen_location)

    for e in epsilons:
        oxygenTraj = spm3d(oxygenData, e, currentDrug)
        np.savetxt(join(savePathCompressed,"oxygenTrajspm3dEpsilon-"+str(e)+".csv"), oxygenTraj, delimiter=",")


        oxygenTime = oxygenTraj[:,3]
        oxygenInterpolated = doInterpolation(oxygenTraj)
        np.savetxt(join(savePathInterpolatedTimes,"oxygenTrajspm3dinterpolated-"+str(e)+".csv"), oxygenTime, delimiter=",")
        np.savetxt(join(savePathInterpolatedData,"oxygenTrajspm3dinterpolated-"+str(e)+".csv"), oxygenInterpolated, delimiter=",")


        for i in range(216): # for the 216 nitrogen trajectories
            nitrogenTrajectory = spm3d(nitrogenData, e, i)
            np.savetxt(join(savePathCompressed,"nitrogenTrajspm3dEpsilon-"+str(e)+"--"+str(i)+".csv"), nitrogenTrajectory, delimiter=",")

            nitrogenTime = nitrogenTrajectory[:,3]
            union = np.union1d(oxygenTime, nitrogenTime)
            nitrogenInterpolated = doInter(nitrogenTrajectory, union)
            nitrogenInterpolated = nitrogenInterpolated[nitrogenInterpolated[:,3].argsort()]

            np.savetxt(join(savePathInterpolatedTimes,"nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenTime, delimiter=",")
            np.savetxt(join(savePathInterpolatedData,"nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenInterpolated, delimiter=",")
    
    print("Compressed and interpolated for SPM.")

    for e in epsilons:
        oxygenTraj = rdpOxygen(oxygenData, e)
        oxygenTime = oxygenTraj[:,3]

        oxygenInterpolated = doInterpolation(oxygenTraj)
        np.savetxt(join(savePathInterpolatedTimes,"oxygenTrajrdpinterpolated-"+str(e)+".csv"), oxygenTime, delimiter=",")
        np.savetxt(join(savePathInterpolatedData,"oxygenTrajrdpinterpolated-"+str(e)+".csv"), oxygenInterpolated, delimiter=",")


        for i in range(216): # for the 216 nitrogen trajectories
            nitrogenTrajectory = performRDP(nitrogenData, e, i)
            np.savetxt(join(savePathCompressed,"nitrogenTrajrdpEpsilon-"+str(e)+"--"+str(i)+".csv"), nitrogenTrajectory, delimiter=",")

            nitrogenTime = nitrogenTrajectory[:,3]
            union = np.union1d(oxygenTime, nitrogenTime)
            nitrogenInterpolated = doInter(nitrogenTrajectory, union)
            nitrogenInterpolated = nitrogenInterpolated[nitrogenInterpolated[:,3].argsort()]

            np.savetxt(join(savePathInterpolatedTimes,"nitrogenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenTime, delimiter=",")
            np.savetxt(join(savePathInterpolatedData,"nitrogenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), nitrogenInterpolated, delimiter=",")

    for e in epsilons:
        oxygenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajspm3dinterpolated-"+str(e)+".csv"), delimiter=",")

        union = np.array([])
        for i in range(216):
            nitrogenData = np.loadtxt(join(savePathInterpolatedData, "nitrogenTrajspm3dinterpolated-"+str(e)+"--"+str(i)+".csv"), delimiter=",")

            oD = oxygenData[nitrogenData[:,3].astype(int)][:,:3]
            timesteps = nitrogenData[:,3].astype(int)
            nit = nitrogenData[:,:3]

            timeSteps = giveTimesteps(oD, nit, timesteps, 0)
            union = np.union1d(union, timeSteps)
            np.savetxt("distances for spm epsilon-"+str(e) , union, delimiter=',')

    for e in epsilons:
        oxygenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajrdpinterpolated-"+str(e)+".csv"), delimiter=",")

        union = np.array([])
        for i in range(216):
            nitrogenData = np.loadtxt(join(savePathInterpolatedData, "oxygenTrajrdpinterpolated-"+str(e)+"--"+str(i)+".csv"), delimiter=",")

            oD = oxygenData[nitrogenData[:,3].astype(int)][:,:3]
            timesteps = nitrogenData[:,3].astype(int)
            nit = nitrogenData[:,:3]

            timeSteps = giveTimesteps(oD, nit, timesteps, 0)
            union = np.union1d(union, timeSteps)
            np.savetxt("distances for rdp epsilon-"+str(e) , union, delimiter=',')

