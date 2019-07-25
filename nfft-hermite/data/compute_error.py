import sys
import numpy as np
import matplotlib.pyplot as plt

if not len(sys.argv) == 4:
	raise Exception("compute_error.py expects 3 arguments: paths for x, y1, and y2 (y2 is base of relative error)")

with open(sys.argv[1], "r") as f:
	x = [float(n) for n in f.read().split("\n")]

with open(sys.argv[2], "r") as f:
	y1 = []
	for n in f.read().split("\n"):
		try:
			y1.append(float(n))
		except ValueError:
			y1.append(None)

with open(sys.argv[3], "r") as f:
	y2 = []
	for n in f.read().split("\n"):
		try:
			y2.append(float(n))
		except ValueError:
			y2.append(None)

if not len(x) == len(y1) == len(y2):
	raise Exception("compute_error.py expects all three datasets to be of equal size")

diff = []
for i in range(len(y1)):
	if not (y1[i] == None or y2[i] == None):
		diff.append(abs(y1[i] - y2[i]))
	else:
		diff.append(None)

pdiff = []
for i in range(len(y1)):
	if not (y1[i] == None or y2[i] == None):
		pdiff.append(abs((y1[i] - y2[i]) / (y2[i])) * 100)
	else:
		pdiff.append(None)

rdiff = []
for i in range(len(y1)):
	if not (y1[i] == None or y2[i] == None):
		rdiff.append(abs(y1[i]) / (abs(y2[i] + 10**(-8))))

plt.plot(x, diff)
# axes = plt.gca()
# axes.set_ylim(0, max(diff))
plt.title("Absolute error between {} and {}".format(sys.argv[2], sys.argv[3]))
plt.show()

plt.plot(x, pdiff)
plt.title("Percentage error between {} and {} relative to latter".format(sys.argv[2], sys.argv[3]))
plt.show()

plt.plot(x, rdiff)
plt.title("Relative error between {} and {}".format(sys.argv[2], sys.argv[3]))
plt.show()

plt.scatter(x, y1)
plt.show()

print(pdiff.index(max(pdiff)))
print(y1[pdiff.index(max(pdiff))], y2[pdiff.index(max(pdiff))])

'''
file1 = "hermpoly_asy_airy_out.mat"

f = open(file1, "r")
content = f.read()

content = content.split("\n")
y1 = [float(l.split(",")[1]) if len(l.split(",")) > 1 and not l.split(",")[1] == "" else None for l in content]

f.close()

###

file2 = "mathematica_Yp.txt"

f = open(file2, "r")
content = f.read()

y2 = [float(n) for n in content.split("\n")]

f.close()

###

file3 = "hermpoly_asy_airy_x.mat"

f = open(file3, "r")
content = f.read()

xs = [float(n) for n in content.split("\n")]

f.close()
'''