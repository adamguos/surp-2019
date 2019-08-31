import numpy as np
import matplotlib.pyplot as plt

mut = np.loadtxt("mut.mat")
mutn = np.loadtxt("mutn.mat")

plt.plot(mut, mutn)
plt.savefig("plot.svg")