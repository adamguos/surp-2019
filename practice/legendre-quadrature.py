# https://www.cs.purdue.edu/archives/2001/wxg/codes/ORTHPOL

import numpy as np
from numpy import linalg as LA
from math import sqrt

# def Legendre(n, x):
#     x = np.array(x)
#     if (n == 0):
#         return x * 0 + 1.0
#     elif (n == 1):
#         return x
#     else:
#         return ((2.0 * n - 1.0) * x * Legendre(n - 1, x) - (n - 1) * Legendre(n - 2, x)) / n

# def DLegendre(n, x):
#     x = np.array(x)
#     if (n == 0):
#         return x * 0
#     else:
#         return n / (x ** 2 - 1) * (x * Legendre(n, x) - Legendre(n - 1, x))

def Legendre(n):
    if (n == 0):
        return np.poly1d([1])
    elif (n == 1):
        return np.poly1d([1, 0])
    else:
        return ((2.0 * n - 1.0) * np.poly1d([1, 0]) * Legendre(n - 1) - (n - 1) * Legendre(n - 2)) / n

def DLegendre(n):
    return Legendre(n).deriv()

def nodes(n):
    return Legendre(n).r

def weights(n):
    W = []
    X = nodes(n)
    P = DLegendre(n)
    for i in range(n):
        W.append(2 / ((1 - X[i]**2) * P(X[i])**2))
    return W

def integrate(a, b, n, f):
    X = nodes(n)
    print(X)
    W = weights(n)
    print(W)
    r = 0

    for i in range(n):
        r += W[i] * f((b-a)/2 * X[i] + (a+b)/2)
        print(i, r)
    r *= (b-a)/2

    return r

def monicise(p):
    return p / p.c[0]

def companion(p):
    p = monicise(p)
    n = len(p.c) - 1
    m = np.insert(np.eye(n-1), 0, 0, axis=1)
    coefs = []
    for c in list(reversed(p.c))[:-1]:
        coefs.append(c * -1.0)
    return np.insert(m, n-1, coefs, axis=0)

def rms(num):
    return sqrt(sum(n*n for n in num)/len(num))

def cjm(A):
    # initialise
    N = A[0].size
    b = np.zeros(N)
    xi = np.random.randn(N)
    r = b - A.dot(xi)
    p = r
    x = xi
    Xr = []

    # iterate
    for k in range(10000):
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
        Xr.append(rms(b - A.dot(x)))

        print(k, LA.norm(b - A.dot(x)))
    
    return x

