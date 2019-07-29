# Imports
from adc import adc
import numpy as np

# Problem size
N = 1000

# Initialise random input
vecD0 = np.random.rand(N)
vecD1 = np.random.rand(N-1)
x = np.random.rand(N)

# Run adc and print output
print(adc(vecD0, vecD1, x))