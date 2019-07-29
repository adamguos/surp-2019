import os
import sys
import time
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

if not len(sys.argv) == 3:
	raise Exception("gauss_transform.py expects 3 arguments: file of X, directory of H")

with open(sys.argv[1], "r") as f:
	X = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")

def gram_sum(H, W, m, l):
	y = np.float64(0)
	for j in range(W.shape[0]):
		y += H[m, j] * H[l, j] * W[j]
	return y

for file in sorted(os.listdir(sys.argv[2]))[3:8]:
	degree = int(file.split(".")[0][1:])
	x = X[int(degree / 32 - 1)]

	start = time.time()

	with open(os.path.join(sys.argv[2], file), "r") as f:
		H = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")
		
		n = H.shape[0]
		j_max = x.shape[0]

		print("Computing Psi for degree {}...".format(degree), end=" ", flush=True)
		Psi = np.fromfunction(lambda k, j: H[0,j] * H[k,j], (n, j_max - 1), dtype=int)
		print("done.")

		R = np.empty([H.shape[0]])
		R[0] = np.nan

		for n in range(1, H.shape[0]):
			e = np.zeros(n)
			e[0] = 1

			print("Computing W for n = {} (degree {})...".format(n, degree), end=" ", flush=True)

			W = np.linalg.lstsq(Psi[:n], e)[0]
			print("done. Checking Gram/identity norm...", end=" ", flush=True)

			G = np.fromfunction(lambda m, l: gram_sum(H, W, m, l), (int(n/2), int(n/2)), dtype=int)
			try:
				r = LA.norm(G - np.identity(G.shape[0]), np.inf)
			except ValueError:
				r = np.nan
			print(r)

			R[n] = r

			if n == 127:
				plt.imshow(G)
				plt.show()
		
		plt.scatter(range(H.shape[0]), R)
		plt.xlabel("Degree of Hermite function")
		plt.ylabel("Infinity norm of Gram matrix subtract identity matrix")
		plt.title("Gram matrix residual of Hermite function up to n = {}".format(degree))
		plt.show()