"""
Problem 6
Plotting the three eigenvectors corresponding to the three lowest 
eigenvalues for n = 10 and n = 100 (discretization steps) from Jacobi algorithm
and analytical solution.
The plot shows vector element v_i against x_i including boundary points. 
"""

import numpy as np
import matplotlib.pyplot as plt
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

for n in [10, 100]:
    data = np.loadtxt(f"../data/problem6_n{n}.txt", skiprows=1)
    x = data[:,0]
    v1_num, v1_ana = data[:,1], data[:,2]
    v2_num, v2_ana = data[:,3], data[:,4]
    v3_num, v3_ana = data[:,5], data[:,6]
    
    plt.figure(figsize=(10,6))
    plt.plot(x, v1_num, color="cornflowerblue", label=r"$v_1$ numerical")
    plt.plot(x, v1_ana, "--", color="midnightblue", label=r"$v_1$ analytical")
    plt.plot(x, v2_num, color="mediumorchid", label=r"$v_2$ numerical")
    plt.plot(x, v2_ana, "--", color="darkmagenta", label=r"$v_2$ analytical")
    plt.plot(x, v3_num, color="pink", label=r"$v_3$ numerical")
    plt.plot(x, v3_ana, "--", color="deeppink",label=r"$v_3$ analytical")

    plt.xlabel(r"$x_i$", fontsize=14)
    plt.ylabel(r"$v_i$", fontsize=14)
    plt.title(f"First 3 eigenvectors for n = {n}")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))    
    plt.grid(True)
    plt.tight_layout()
    plt.show()
