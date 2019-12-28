import numpy as np
import matplotlib.pyplot as plt

t = np.loadtxt("t.mat")
tn = np.loadtxt("tn.mat")
# tnl = np.loadtxt("tnl.mat")
T = np.loadtxt("T.mat")
A = np.loadtxt("A.mat")

plt.plot(t, np.abs(tn), color="r", label="direct compute")
# plt.plot(t, np.abs(tnl), color="b", label="ifft")
plt.legend()
for i in range(len(T)):
	plt.axvline(x=T[i], color="k", linestyle="--")
	plt.text(T[i] + 0.02, 0, "t = {}\na = {}".format(T[i], A[i]))
# plt.savefig("plot.png")
plt.savefig("plot.svg")