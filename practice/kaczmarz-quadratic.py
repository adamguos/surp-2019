# https://lu.seas.harvard.edu/files/yuelu/files/kaczmarz_quadratic_final.pdf

n = 1000
m = 1000

from timeit import default_timer as timer

def compute_cpu():

    import numpy as np
    import math

    # generate a positive semi-definite matrix
    def gen_PSD(size):
        A = np.random.rand(size, size)
        return A.dot(A.transpose())

    # compute an observation of sigma with some random measurement matrix
    def compute_observation(sigma):
        n = sigma.shape[0]
        a = np.random.rand(n, 1)
        return a.transpose().dot(sigma).dot(a)[0][0], a

    # generate m observations of sigma
    def gen_observations(m, sigma):
        out = [compute_observation(sigma) for i in range(m)]
        obvs = [o[0] for o in out]
        A = [o[1] for o in out]
        return obvs, A

    print("Generating initial values...")
    sigma = gen_PSD(n)
    obvs, A = gen_observations(m, sigma)
    Y = [math.sqrt(a) for a in obvs]
    r = np.linalg.matrix_rank(np.linalg.cholesky(sigma))
    T = m

    # intial guess
    U = np.random.rand(n, r)
    print("Beginning computation...")
    start = timer()

    for k in range(T):
        U0 = U
        U = U0 - ((np.linalg.norm(U0.transpose().dot(A[k])) - Y[k]) / (np.linalg.norm(U0.transpose().dot(A[k]))) * A[k].dot(A[k].transpose()).dot(U0) / np.power(np.linalg.norm(A[k]), 2))
        if k % 100 == 0:
            print("Iteration {0} out of {1} -- time elapsed: {2} s".format(k, T, timer() - start))
            start = timer()

    print(U)
    print("RMS error: {0}".format(np.sqrt(np.mean(np.square(sigma - U.dot(U.transpose()))))))

def compute_gpu():

    import pycuda.autoinit
    import pycuda.gpuarray as gpuarray
    import numpy as np
    import skcuda.linalg as linalg
    import skcuda.misc as misc
    import math

    linalg.init()

    # multiply matrices a and b using CUDA
    def gpu_mult(a, b):
        a_gpu = gpuarray.to_gpu(a)
        b_gpu = gpuarray.to_gpu(b)
        c_gpu = linalg.dot(a_gpu, b_gpu)
        return c_gpu.get()
    
    def gpu_transpose(a):
        a_gpu = gpuarray.to_gpu(a)
        at_gpu = linalg.transpose(a_gpu)
        return at_gpu.get()

    # generate a positive semi-definite matrix
    def gen_PSD(size):
        A = np.random.rand(size, size)
        return gpu_mult(A, gpu_transpose(A))

    # compute an observation of sigma with some random measurement matrix
    def compute_observation(sigma):
        n = sigma.shape[0]
        a = np.random.rand(n, 1)
        b = gpu_mult(gpu_transpose(a), sigma)
        return gpu_mult(b, a)[0][0], a

    # generate m observations of sigma
    def gen_observations(m, sigma):
        out = []
        start = timer()
        for i in range(m):
            out.append(compute_observation(sigma))
            if i % 1000 == 0:
                print("Observation {0} out of {1} -- time elapsed: {2} s".format(i, m, timer() - start))
                start = timer()

        obvs = [o[0] for o in out]
        A = [o[1] for o in out]
        return obvs, A
    
    print("Generating initial values...")
    sigma = gen_PSD(n)
    obvs, A = gen_observations(m, sigma)
    Y = [math.sqrt(a) for a in obvs]
    r = np.linalg.matrix_rank(np.linalg.cholesky(sigma))
    T = m

    # intial guess
    U = np.random.rand(n, r)
    print("Beginning computation...")
    start = timer()

    for k in range(T):
        U0 = U
        U = U0 - ((np.linalg.norm(gpu_mult(gpu_transpose(U0), A[k])) - Y[k]) / (np.linalg.norm(gpu_mult(gpu_transpose(U0), A[k]))) * gpu_mult(gpu_mult(A[k], gpu_transpose(A[k])), U0) / np.power(np.linalg.norm(A[k]), 2))
        if k % 100 == 0:
            print("Iteration {0} out of {1} -- time elapsed: {2} s".format(k, T, timer() - start))
            start = timer()

    print(U)
    print("RMS error: {0}".format(np.sqrt(np.mean(np.square(sigma - U.dot(U.transpose()))))))

compute_gpu()
# compute_cpu()