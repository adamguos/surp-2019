import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt("t.mat")
ft = np.loadtxt("ft.mat")

plt.plot(t, ft)
plt.savefig("plot.svg")