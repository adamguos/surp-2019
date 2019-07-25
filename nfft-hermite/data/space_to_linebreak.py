import sys

if len(sys.argv) < 2:
	raise Exception("space_to_linebreak.py expects at least 1 argument")

for file in sys.argv[1:]:
	if not sum(1 for line in open(file)) == 1:
		raise Exception("space_to_linebreak.py expects files with 1 line of data points, separated by spaces")
	
	with open(file, "r") as f:
		input = f.read()
	
	with open(file, "w") as f:
		f.write("\n".join(input.split(" ")))
	
	print("%s reformatted" % file)