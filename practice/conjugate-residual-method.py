# https://en.wikipedia.org/wiki/Conjugate_residual_method
# question: complexity of matrix-vector multiplication is O(n^2)? is this a problem?
# slower and requires more storage than conjugate gradient method?

from math import sqrt
import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

# compute root-mean-square
def rms(num):
    return sqrt(sum(n*n for n in num)/len(num))

# N = 10  # dimensions of A
N = int(sys.argv[1])    # number of iterations

A = 0.1 * np.random.randn(N, N) + 0.75
# A = (A + A.conjugate()) / 2 # force A to be Hermitian
A = np.diag(np.diag(A))

print(A)

A0 = np.random.randn(N, N)
# A0 = A0.dot(A0.transpose())
A0 = (A0 + A0.transpose()) / 2

values, vectors = LA.eig(A0)

newmatrix = vectors.dot(A).dot(vectors.transpose())
print(newmatrix)
print(LA.cond(newmatrix))

A = newmatrix

# replace diagonal matrix from decomposition with generated diagonal matrix

b = np.random.randn(N, 1)
x0 = np.random.randn(N, 1)   # arbitrary initial guess

# print(LA.cond(A))   # minimise condition number to converge faster/more effectively

r = b - A.dot(x0)  # error of initial guess
p = r

x = x0  # initialise iterated solution
X = []    # list of all RMSs

for k in range(N):
    alpha = np.asscalar(r.transpose().dot(A).dot(r)) / np.asscalar((A.dot(p)).transpose().dot(A).dot(p))
    x2 = x + alpha * p
    r2 = r - alpha * A.dot(p)
    beta = np.asscalar(r2.transpose().dot(A).dot(r2)) / np.asscalar(r.transpose().dot(A).dot(r))
    p2 = r2 + beta * p

    # iterate variables
    x = x2
    r = r2
    p = p2

    # add current solution to list of solutions
    X.append(rms(b - A.dot(x)))

    print(k, LA.norm(b - A.dot(x)))

# r0 = b - A.dot(x0)
# r = b - A.dot(x)
# print(r0)
# print(r)
# print(rms(r0))
# print(rms(r))

plt.plot(X)
plt.ylabel("rms error")
plt.xlabel("iteration")
plt.show()