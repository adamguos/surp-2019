'''
refactor:
compute independent of experiments -- algorithms are constant and reusable
	gauss_transform takes X (1xs), H (mxs), N (1xn), m=max degree, s=num of samples, N=degrees to compute W for
	returns W (nx(s-1)), each row is a single set of W_j that corresponds to degree in N, R (nx1) Gram residuals corresponding to sets of W_j
compute can take vectors or single values
experiments independent of data
'''

import os
import sys
import time
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

'''
n = 2

with open("../../../../Documents/X.mat", "r") as f:
	X = np.array([a for a in f.read().split("\n")[n-1].split("\t")])

with open("../../../../Documents/H512_{}.mat".format(str(n).zfill(4)), "r") as f:
	H = np.array([a for a in f.read().split("\n")[-1].split("\t")])

plt.scatter(X, H)
plt.show()
'''

if not len(sys.argv) == 3:
	raise Exception("gauss_transform.py expects 3 arguments: file of X, direcctory of H")

with open(sys.argv[1], "r") as f:
	X = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")

def gram_sum(H, W, m, l):
	y = np.float64(0)
	for j in range(W.shape[0]):
		y += H[m, j] * H[l, j] * W[j]
	return y

degree = len(open(os.path.join(sys.argv[2], os.listdir(sys.argv[2])[0]), "r").read().split("\n"))

R = np.empty([len(os.listdir(sys.argv[2])), degree])

for file in sorted(os.listdir(sys.argv[2])):
	trial = int(file.split("_")[1].split(".")[0])
	x = X[trial - 1]

	start = time.time()

	with open(os.path.join(sys.argv[2], file), "r") as f:
		H = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")
		j_max = x.shape[0]

		print("Computing Psi for trial {}...".format(trial), end=" ", flush=True)
		Psi = np.fromfunction(lambda k, j: H[0,j] * H[k,j], (degree, j_max - 1), dtype=int)
		print("done.")

		for n in range(H.shape[0]):
			if n % 2 == 0:
				R[trial - 1, n] = np.nan
			else:
				e = np.zeros(n)
				e[0] = 1

				print("Computing W for n = {} (trial {})...".format(n, trial), end=" ", flush=True)

				W = np.linalg.lstsq(Psi[:n], e)[0]
				print("done. Checking Gram/identity norm...", end=" ", flush=True)

				G = np.fromfunction(lambda m, l: gram_sum(H, W, m, l), (int(n/2), int(n/2)), dtype=int)
				try:
					r = LA.norm(G - np.identity(G.shape[0]), np.inf)
				except ValueError:
					r = np.nan
				print(r)

				R[trial - 1, n] = r

				if n == 191:
					plt.imshow(G)
					plt.show()
	
	print("Time elapsed for trial {}: {}".format(trial, time.time() - start))

np.savetxt("../data/gauss_transform/R.mat", R)

Rm = np.array([np.nanmean(R[:, i]) for i in range(R.shape[1])])
plt.scatter(np.linspace(0, Rm.shape[0] - 1, Rm.shape[0]), Rm)
plt.show()

'''
file = "H512_0001.mat"
with open(os.path.join(sys.argv[2], file), "r") as f:
	trial = int(file.split("_")[1].split(".")[0])
	x = X[trial - 1]

	H = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")
	j_max = x.shape[0]

	R = np.array([])
	
	for n in range(H.shape[0]):
		if n % 2 == 0:
			R = np.append(R, np.nan)
		else:
			print("Computing W for n={}...".format(n), end=" ", flush=True)
			# Psi = np.zeros([n, j_max-1])
			# for k in range(Psi.shape[0]):
			# 	for j in range(Psi.shape[1]):
			# 		Psi[k,j] = (x[j] - x[j+1]) * H[0,j] * H[k,j]
			# Psi = np.fromfunction(lambda k, j: (x[j] - x[j+1]) * H[0,j] * H[k,j], (n, j_max - 1), dtype=int)
			Psi = np.fromfunction(lambda k, j: H[0,j] * H[k,j], (n, j_max - 1), dtype=int)
			e = np.zeros(n)
			e[0] = 1

			W = np.linalg.lstsq(Psi, e)[0]
			print("done. Checking accuracy...", end=" ", flush=True)
			print("{}. Checking Gram/identity norm...".format(np.allclose(np.matmul(Psi, W), e)), end=" ", flush=True)

			# G = np.zeros([int(n/2), int(n/2)])
			# for m in range(G.shape[0]):
			# 	for l in range(G.shape[1]):
			# 		G[m,l] = gram_sum(H, W, m, l)
			G = np.fromfunction(lambda m, l: gram_sum(H, W, m, l), (int(n/2), int(n/2)), dtype=int)
			r = LA.norm(G - np.identity(G.shape[0]))
			# r = LA.norm(np.matmul(Psi, W) - e)
			print(r)

			if n == 33:
				plt.imshow(G)
				plt.show()
			
			R = np.append(R, r)
	
	np.savetxt("../data/gauss_transform/R.mat", R)
	
	plt.scatter(range(H.shape[0]), R)
	plt.show()
'''

'''
with open(sys.argv[2], "r") as f:
	H = np.array([[k for k in j.split("\t")] for j in f.read().split("\n")], dtype="float64")

n = H.shape[0]
j_max = X.shape[0]

Psi = np.fromfunction(lambda k, j: (X[j] - X[j+1]) * H[0,j] * H[k,j], (n, j_max - 1), dtype=int)
e = np.zeros(n)
e[0] = 1

print("Computing W...")
W = np.linalg.lstsq(Psi, e)[0]
print("W computed.")
# np.savetxt("../data/gauss_transform/W.mat", W)
# print("W saved.")

print(np.allclose(np.dot(Psi, W), e))
'''