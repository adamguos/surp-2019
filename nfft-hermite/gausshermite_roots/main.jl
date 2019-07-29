include("generate_nodes_weights.jl")
include("hermite_gram_residual.jl")

using Plots

# generate_weights(range(32, stop=1024, step=32), "data/weights")

weights_dir = "data/weights"
hermites_dir = "data/hermites"

weights = sort(readdir(weights_dir))
hermites = sort(readdir(hermites_dir))

iterator = size(weights)[1] < size(hermites)[1] ? weights : hermites

R = []

for i in 1:size(iterator, 1)
	W = readdlm("$weights_dir/$(weights[i])", '\t', BigFloat, '\n')
	H = readdlm("$hermites_dir/$(hermites[i])", '\t', BigFloat, '\n')
	r = hermite_gram_residual(W, H)
	push!(R, r)
end

N = []

for i in iterator
	n = parse(Int, split(split(i, "_")[2], ".")[1])
	push!(N, n)
end

display(plot(N, R, line=:scatter))