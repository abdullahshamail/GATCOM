import numpy as np
from utils import donor_acceptor_dist, loadData

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\oxygenTraj.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogenTraj.txt"

    oxygen = loadData(oxygen_location)
    nitrogen = loadData(nitrogen_location)

    a, b = donor_acceptor_dist(oxygen[:,0,:], nitrogen, np.arange(200000))

    np.savetxt(r"original trajectory without compression.csv", a, delimiter=",")