import sys

if not len(sys.argv) == 2:
	raise Exception("trim_octave_data.py expects 1 argument: path to file")

with open(sys.argv[1], "r") as f:
	input = f.read()

output = input.split("\n")[5][1:]

with open(sys.argv[1], "w") as f:
	f.write(output)

print("{} trimmed".format(sys.argv[1]))

'''
file1 = "hermpoly_asy_airy_new.mat"

f = open(file1, "r")
content = f.read()

ys = ""

content = content.split(" ")

for c in content:
	re = c.split(",")[0][1:]
	im = c.split(",")[1][:-1]
	
	# only keep non-complex values
	if im == "0":
		ys += re
	
	# delineate entries by newline
	ys += "\n"

###

file2 = "hermpoly_asy_airy_x.mat"

f = open(file2, "r")
content = f.read()

xs = ""

content = content.split(" ")

for c in content:
	xs += c
	xs += "\n"

###

print(ys.count("\n"))
print(xs.count("\n"))

ys = ys.split("\n")
xs = xs.split("\n")

out = ""

for i in range(len(ys) - 1):
	out += xs[i] + "," + ys[i]
	out += "\n"

print(out)

fileo = open("hermpoly_asy_airy_out.mat", "w")
fileo.write(out)
fileo.close()
'''