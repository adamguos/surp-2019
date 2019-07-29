'''
Reference:	M. Gu and S. C. Eisenstat, A Divide-and-Conquer Algorithm For The Symmetric Tridiagonal Eigenproblem, 1992
Goal:		evaluate the eigendecomposition of a symmetric tridiagonal matrix in O(n log n) time
			T = U*D*Ut

Input:		vecD0 = main diagonal of symmetric tridiagonal matrix T
			vecD1 = one-off diagonal of T
			x = vector applied to computation of Ut*x
Output:		y = Ut*x
'''

# Imports
from tree import Tree
import numpy as np
from scipy import sparse

# ADC main function
def adc(vecD0, vecD1, x):
	T = sparse.diags([vecD0, vecD1, vecD1], [0, -1, 1], format="csr")	# generate T = sym. tri. diag. matrix with main diagonal vecD0 and one-off diagonal vecD1
																		# note: "csr" is used to enable subscripting but may be less efficient than the default, "dia"
	ADC_tree = Tree(T)													# instantiate main ADC tree; each node contains the (sub)matrix T at that level
	min_n = 4															# size of submatrix to be processed directly using scipy.sparse.linalg.eigs

	return construct_node(T, min_n)

# given T, returns a tuple of two nodes: (left, right)
# recursively creates nodes for each of left and right until size of submatrix <= min_n
def construct_node(T, min_n):
	N = T.shape[0]		# get size of matrix at current node

	# if N > min_n, split T into two submatrices and construct their nodes
	# else, just return T itself with no children, i.e. (None, None)
	if N > min_n:
		# split T into submatrices T1 and T2
		k = int(N / 2)		# k is the index of the row & column that separates T1 and T2
		T1 = T[:k, :k]		# T1.shape = (k, k), T2.shape = (N - k - 1, N - k - 1)
		T2 = T[-(N-k-1):, -(N-k-1):]

		# get elements in the dividing row and column between T1 and T2
		alpha = T[k, k]
		beta1 = T[k-1, k]
		beta2 = T[k+1, k]

		return (T1.shape, T2.shape)
	else:
		return 0