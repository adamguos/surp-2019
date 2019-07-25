# https://en.wikipedia.org/wiki/Kaczmarz_method

from math import sqrt
import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

# compute root-mean-square
def rms(num):
    return sqrt(sum(n*n for n in num)/len(num))

N = 5  # dimensions of A
n = int(sys.argv[1])    # number of iterations

A = np.random.randn(N, N)
xa = np.random.randn(N, 1)  # "correct" answer
b = A.dot(xa)
x0 = np.random.randn(N, 1)   # arbitrary initial guess
# x0 = np.zeros((N, 1))

x = x0  # initialise iterated solution
X = []  # list of all iterated solution RMSs

for k in range(n):
    i = k % N

    numerator = b[i,0] - np.inner(A[i], x.transpose())[0]
    denominator = np.power(LA.norm(A[i]), 2)
    conjugate = A[i].conj()
    term = np.reshape(numerator / denominator * conjugate, (-1, 1))

    x2 = x + term
    
    x = x2
    X.append(rms(xa - x))

    # x2 = x + ((np.asscalar(b[i]) - np.asscalar(np.inner(A[i], x.transpose()))) / np.power(LA.norm(A[i]), 2) * A[i].conj()).transpose()
    
    # x = x2
    # print(x)
    # X.append(rms(b - A.dot(x)))

plt.plot(X)
plt.ylabel("rms error")
plt.xlabel("iteration")
plt.show()