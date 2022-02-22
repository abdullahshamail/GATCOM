from os import listdir
from os.path import isfile, join
import numpy as np
from utils import loadData, spm3d, rdpOxygen, doInter, doInterpolation, giveTimesteps


def compression(timeStepsFolder, methods, comTypes, epsilons):
    for method in methods:
        for epsilon in epsilons:
            for comType in comTypes:
                print("Working on", method, epsilon, comType)
                tempFolder = join(timeStepsFolder, comType)
                timeStepFiles = [f for f in listdir(tempFolder) if isfile(join(tempFolder, f))]
                for f in timeStepFiles:
                    if f.endswith('.csv'):
                        originalTimesteps = np.loadtxt(join(tempFolder, f))
                    
                        rflOxygen = oxygenData[:,0,:]
                        originalTimesteps = originalTimesteps.astype(int)
                        toappend = np.arange(200000)
                        toappend = toappend.reshape(toappend.shape[0], 1)
                        rflOxygen = np.append(rflOxygen, toappend, axis = 1)
                        rflOxygen = rflOxygen[originalTimesteps,:]
                        if (method == "spm"):
                            oxygenTraj = spm3d(rflOxygen, float(epsilon), 0)
                        else: 
                            oxygenTraj = rdpOxygen(rflOxygen, float(epsilon))

                        np.savetxt(join(r"compressed\oxygens","oxygen"+method+epsilon+comType+f), oxygenTraj, delimiter=",")

                        for i in range(216):
                            
                            nitrogeni = nitrogenData[:,i,:]
                            toappend = np.arange(200000)
                            toappend = toappend.reshape(toappend.shape[0], 1)
                            nitrogeni = np.append(nitrogeni, toappend, axis = 1)
                            nit = nitrogeni[originalTimesteps,:]

                            if (method == "spm"):
                                nitrogenTrajectory = spm3d(nit, float(epsilon), i)
                            else:
                                nitrogenTrajectory = rdpOxygen(nit, float(epsilon)) #note - rdp oxygen used here is deliberate, since that is the function we want to use here
                            np.savetxt(join(r"compressed\nitrogens","nitrogen"+method+epsilon+comType+"-"+str(i)+f), nitrogenTrajectory, delimiter=",")

def interpolationLocal(timeStepsFolder, methods, comTypes, epsilons):
    for method in methods:
        for epsilon in epsilons:
            for comType in comTypes:
                print("Working on", method, epsilon, comType)
                tempFolder = join(timeStepsFolder, comType)
                timeStepFiles = [f for f in listdir(tempFolder) if isfile(join(tempFolder, f))]
                for f in timeStepFiles:


                    oxygendata = np.loadtxt(join(r"compressed\oxygens","oxygen"+method+epsilon+comType+f), delimiter=",")
                    oxygenTime = oxygendata[:,3]
                    oxygenInterpolated = doInterpolation(oxygendata)


                    np.savetxt(join(r"interpolated\oxygens","oxygen"+method+epsilon+comType+f),oxygenInterpolated, delimiter=",")

                    for i in range(216):
                        nitrogenTrajectory = np.loadtxt(join(r"compressed\nitrogens","nitrogen"+method+epsilon+comType+"-"+str(i)+f), delimiter=",")
                        nitrogenTime = nitrogenTrajectory[:,3]
                        union = np.union1d(oxygenTime, nitrogenTime)

                        nitrogenInterpolated = doInter(nitrogenTrajectory, union)
                        print(nitrogenInterpolated.shape)

                        np.savetxt(join(r"interpolated\nitrogens","nitrogen"+method+epsilon+comType+"-"+str(i)+f), nitrogenInterpolated, delimiter=",")

def calcDistances(timeStepsFolder, methods, comTypes, epsilons):
    for method in methods:
        for epsilon in epsilons:
            for comType in comTypes:
                print("Working on", method, epsilon, comType)
                tempFolder = join(timeStepsFolder, comType)
                timeStepFiles = [f for f in listdir(tempFolder) if isfile(join(tempFolder, f))]
                for f in timeStepFiles:
                    if f.endswith('.csv'):
                        oxygenData = np.loadtxt(join(r"interpolated\oxygens","oxygen"+method+epsilon+comType+f), delimiter=",")
                        # nitrogenFolder = join(mainFolder, method, epsilon)
                        # files = [f for f in listdir(nitrogenFolder) if isfile(join(nitrogenFolder, f))]
                        union = np.array([])
                        for adjustment in [0, float(epsilon)]:
                            for i in range(216):
                                nitrogenData = np.loadtxt(join(r"interpolated\nitrogens","nitrogen"+method+epsilon+comType+"-"+str(i)+f), delimiter=",")

                                filt = np.asarray(nitrogenData[:,3])
                                oD = oxygenData[np.in1d(oxygenData[:,3], filt)]

                                timesteps = nitrogenData[:,3].astype(int)

                                timeSteps = giveTimesteps(oD, nitrogenData, timesteps, adjustment)
                                union = np.union1d(union, timeSteps)
                            np.savetxt(join(r"data\hybrid gatcom", method, epsilon, "distances with adjustment"+ str(adjustment)+" for"+f), union ,delimiter=",")

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\oxygenTraj.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogenTraj.txt"

    oxygenData = loadData(oxygen_location)
    nitrogenData = loadData(nitrogen_location)

    timeStepsFolder = r"Timesteps  with Generalization" #these are the timesteps where d_delta < 3.5 
    comTypes = ["OS atom as COM"]
    epsilons = ["0.5", "0.75", "1.0"]
    methods = ["rdp", "spm"]

    compression(timeStepsFolder, methods, comTypes, epsilons)
    interpolationLocal(timeStepsFolder, methods, comTypes, epsilons)
    calcDistances(timeStepsFolder, methods, comTypes, epsilons)
    