#=
Generates Gauss Hermite nodes and weights and saves them to data/nodes.mat and data/weights.mat respectively
Input: n = degree of Gauss quadrature
=#

using DelimitedFiles
using FastGaussQuadrature

nodes, weights = gausshermite(parse(Int, ARGS[1]))

writedlm("data/nodes.mat", nodes, '\n')
writedlm("data/weights.mat", weights, '\n')