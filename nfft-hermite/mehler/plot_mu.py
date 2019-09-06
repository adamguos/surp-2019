import numpy as np
import matplotlib.pyplot as plt

mut = np.loadtxt("mut.mat")
# mut = mut[:-40]
mutn = np.loadtxt("mutn.mat")
# mutn = mutn[:-40]

plt.plot(mut, mutn)
plt.savefig("plot.svg")