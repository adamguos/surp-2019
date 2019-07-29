'''
Goal:	compute W_j such that $\\sum_j W_j \\psi_0(x_j) \\psi_k(x_j) = \\delta_{k, 0}$,
		and Gram matrix residuals
Input:	X (1, s): x-coordinates of data samples
		H (m, s): H_{m, s} = \\psi_m(X_s)
		N (n, 1): list of degrees to compute W for
Output:	W (n, s): rows of W_j, corresponding to degrees in N
		R (n, 1): Gram residuals of W_j, corresponding to degrees in N
'''

import numpy as np
import matplotlib.pyplot as plt

# compute sum for each element of Gram matrix
def gram_sum(H, w, m, l):
	y = 0
	for j in range(w.shape[0]):
		y += H[m, j] * H[l, j] * w[j]
	return y

# main function
def solve_weights(X, H, N):
	dims = {
		"s": X.shape[0],
		"m": H.shape[0],
		"n": N.shape[0]
	}

	# initialise output arrays
	W = np.empty([dims["n"], dims["s"]])
	R = np.empty([dims["n"]])

	# compute coefficient matrix
	Psi = np.fromfunction(lambda k, j: H[0, j] * H[k, j], (dims["m"], dims["s"]), dtype=int)

	# iterate through N
	for i, n in enumerate(N):
		e = np.zeros(n)
		e[0] = 1

		# Solve for w
		w = np.linalg.lstsq(Psi[:n], e)[0]

		# Compute Gram matrix
		G = np.fromfunction(lambda m, l: gram_sum(H, w, m, l), (int(n / 2), int(n / 2)), dtype=int)

		# Attempt to compute r (fails if n=0)
		try:
			r = np.linalg.norm(G - np.identity(G.shape[0]), np.inf)	# infinity norm
		except ValueError:
			r = np.nan
		
		# Save computed values
		W[i] = w
		R[i] = r

		# Display Gram matrix (temp)
		# if (n % 64 == 0):
		# 	plt.imshow(G)
		# 	plt.title("n = {}".format(n))
		# 	plt.show()
	
	return W, R