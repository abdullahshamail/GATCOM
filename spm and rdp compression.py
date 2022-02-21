import numpy as np
from utils import doInterpolation, loadData, spm3d, doInter
from os.path import join

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

        # oxygendata = np.loadtxt(join(savePathOriginal,"oxygenTrajspm3dOriginal.csv"), delimiter=",")
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
