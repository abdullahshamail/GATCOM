import numpy as np
import pandas as pd
import time
from functools import reduce
from utils import doInterpolation, loadData, spm3d, doInter
from os.path import join
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        epsilon = sys.argv[1]
    else:
        print("Please add epsilon.")
        sys.exit()
    oxygen_location = r"C:\Users\ashamail\Desktop\data\os_200_1.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogen200k.txt"
    hydrogen_location = r"C:\Users\ashamail\Box\Shared Materials_Abdullah_Hasan\ADBIS PAPER_Everything\EXTENSION\Data\hn_200k_0.txt"
    savePathDistances = r"C:\Users\ashamail\Desktop\GATCOM\data\conventional algorithms\distances"
    savePath = r"C:\Users\ashamail\Desktop\GATCOM\data\conventional algorithms"
    savePathInterpolatedTimes = join(savePath, "times")
    savePathInterpolatedData = join(savePath, "interpolated")
    savePathCompressed = join(savePath, "compressed")
    # epsilons = [0.5, 0.75, 1.0]
    
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


    savePathInterpolatedTimes = join(savePath, method, str(epsilon), "times")
    savePathInterpolatedData = join(savePath, method, str(epsilon), "interpolated")
    savePathCompressed = join(savePath, method, str(epsilon), "compressed")
    print(f"Starting epsilon {epsilon}")
    startTime = time.time()
    oxygenTraj = spm3d(oxygenData, epsilon, currentDrug)
    np.savetxt(join(savePathCompressed,"oxygenTrajspm3dEpsilon-"+str(epsilon)+".csv"), oxygenTraj, delimiter=",")

    oxygenTime = oxygenTraj[:,3]
    oxygenInterpolated = doInterpolation(oxygenTraj)
    np.savetxt(join(savePathInterpolatedTimes,"oxygenTrajspm3dinterpolated-"+str(epsilon)+".csv"), oxygenTime, delimiter=",")
    np.savetxt(join(savePathInterpolatedData,"oxygenTrajspm3dinterpolated-"+str(epsilon)+".csv"), oxygenInterpolated, delimiter=",")

    print("Starting N and H")

    for i in range(216): # for the 216 nitrogen trajectories
        nitrogenTrajectory = spm3d(nitrogenData, epsilon, i)
        np.savetxt(join(savePathCompressed,"nitrogenTrajspm3dEpsilon-"+str(epsilon)+"--"+str(i)+".csv"), nitrogenTrajectory, delimiter=",")

        hydrogenTrajectory = spm3d(nitrogenData, epsilon, i)
        np.savetxt(join(savePathCompressed,"hydrogenTrajspm3dEpsilon-"+str(epsilon)+"--"+str(i)+".csv"), hydrogenTrajectory, delimiter=",")

        nitrogenTime = nitrogenTrajectory[:,3]
        hydrogenTime = hydrogenTrajectory[:,3]
        union = reduce(np.union1d, (oxygenTime, nitrogenTime, hydrogenTime))
        nitrogenInterpolated = doInter(nitrogenTrajectory, union)
        hydrogenInterpolated = doInter(hydrogenTrajectory, union)

        nitrogenInterpolated = nitrogenInterpolated[nitrogenInterpolated[:,3].argsort()]
        hydrogenInterpolated = hydrogenInterpolated[hydrogenInterpolated[:,3].argsort()]


        np.savetxt(join(savePathInterpolatedTimes,"nitrogenTrajspm3dinterpolated-"+str(epsilon)+"--"+str(i)+".csv"), nitrogenTime, delimiter=",")
        np.savetxt(join(savePathInterpolatedData,"nitrogenTrajspm3dinterpolated-"+str(epsilon)+"--"+str(i)+".csv"), nitrogenInterpolated, delimiter=",")

        np.savetxt(join(savePathInterpolatedTimes,"hydrogenTrajspm3dinterpolated-"+str(epsilon)+"--"+str(i)+".csv"), hydrogenTime, delimiter=",")
        np.savetxt(join(savePathInterpolatedData,"hydrogenTrajspm3dinterpolated-"+str(epsilon)+"--"+str(i)+".csv"), hydrogenInterpolated, delimiter=",")
        print(f"Done {i} for epsilon {epsilon}")
        print(f"Done epsilon {epsilon} in {time.time() - startTime} seconds")
    
    print(f"Compressed and interpolated for SPM epsilon {epsilon}")