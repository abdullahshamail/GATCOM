import numpy as np
import time
from utils import donor_acceptor_dist, loadData

if __name__ == "__main__":
    oxygen_location = r"C:\Users\ashamail\Desktop\data\os_200_1.txt"
    nitrogen_location = r"C:\Users\ashamail\Desktop\data\nitrogen200k.txt"

    start = time.time()
    print("Loading data...")
    oxygen = loadData(oxygen_location)
    nitrogen = loadData(nitrogen_location)

    print(f"Loaded data in {time.time() - start} seconds")

    start = time.time()
    a, b = donor_acceptor_dist(oxygen[:,0,:], nitrogen, np.arange(200000))
    print(f"Query processing time {time.time() - start} seconds")
    np.savetxt(r"data/original trajectory without compression.csv", a, delimiter=",")
    print("Output saved.")
