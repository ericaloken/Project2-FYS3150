"""
Problem 5
Plotting iterations as function of N.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

data = np.loadtxt("../data/problem5.txt", skiprows=1)
N = data[:,0]
iterations = data[:,1]

plt.figure()
plt.plot(N, iterations, "o-", label="Similarity transformations", color="mediumvioletred")
plt.xlabel(r"Matrix size $N$", fontsize=14)
plt.ylabel("Iterations", fontsize=14)
plt.title("Similarity transformations required for matrix (NxN)")
plt.legend()
plt.grid(True)
plt.show()