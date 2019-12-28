#=
Goal:	Generate gausshermite nodes for many degrees
Input:	N (1, n): list of degrees to generate nodes for
		dir (str): target directory to save files
Output:	{dir}/nodes_nnnn.mat: n files, each a list of nodes corresponding to N[n]
=#

using FastGaussQuadrature
using DelimitedFiles

function generate_nodes(N, dir)
	for n in N
		nodes, weights = gausshermite(n)
		writedlm("$dir/nodes_$(lpad(n, 4, '0')).mat", nodes, " ")
		println("n = $n saved")
	end
end

function generate_weights(N, dir)
	for n in N
		nodes, weights = gausshermite(n)
		writedlm("$dir/weights_$(lpad(n, 4, '0')).mat", weights, " ")
		println("n = $n saved")
	end
end