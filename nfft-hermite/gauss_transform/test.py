import numpy as np
import matplotlib.pyplot as plt
X = np.array([k for k in open("data/X192.mat", "r").read().split("\t")], dtype="float64")
H = np.array([[k for k in j.split("\t")] for j in open("data/H192_0001.mat", "r").read().split("\n")], dtype="float64")
N = np.array([a for a in range(1, 193)])
from solve_weights import solve_weights
W, R = solve_weights(X, H, N)
print(W)
print(R)

plt.scatter(N, R)
plt.show()