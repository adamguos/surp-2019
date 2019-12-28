#=
cluster.jl

Main function:	cluster(n = 256)
Goal:			Cluster points of |T_n|
Input:			n = degree of Gauss quadrature
Output:			tbd
=#

using DelimitedFiles
using FastGaussQuadrature
using LinearAlgebra
using Plots

include("hfuncinfty.jl")

global w = (sqrt(5) - 1) / 2

function cluster(n = 128)
	# Generate Gauss quadrature nodes and weights of degree n using Townsend library
	# flush(stdout) used to ensure that messages are printed in a timely manner
	print("Generating Gauss quadrature rule... "); flush(stdout)
	x, lambda = gausshermite(2n)
	println("done"); flush(stdout)

	# for further computation, each lambda is multiplied by exp((x_l)^2)
	# enumerate(lambda) returns array of tuples (index, value)
	W = map(l -> l[2] * exp(x[l[1]] ^ 2), enumerate(lambda))

	psi = read_hermite_at_nodes(2n)

	# pick some test values for A and T for use by F()
	A = [1]
	T = [0]

	# compute mu_j by generating array of elements over l and taking the sum
	# mu = map(j -> sum([W[i] * F(A, T, x[i]) * psi[j, i] for i in 1:2n]) * w^(-j), 1:n)
	mu = map(j -> sum([W[i] * F(A, T, x[i]) * psi[j, i] for i in 1:2n]), 1:n)

	t = read_t()
	h = read_hermite_at_t()

	return mu, h, t

	# compute T_n by passing all indices of t
	# T() takes care of matching indices with the correct values of hermite_t
	Tn = map(ti -> T_n(n, ti, h, mu), 1:size(t, 1))

	return t, abs.(Tn)
end

function F(A, T, t)
	# F(t) = (\pi * w)^(-1/2) * exp(-(sqrt(5) / 2)t^2) * \sum_{j=1}^K a_j e^(2tt_j)
	# A = [a_1, ..., a_K]
	# T = [t_1, ..., t_K]
	# F is a function of t; A and T are other parameters

	if size(A) != size(T)
		error("F expects parameters A and T to be of equal size K")
	end

	K = size(A, 1)

	return map(tj -> (pi * w)^(-1/2) * exp(-(sqrt(5) / 2) * tj^2) * sum([A[j] * exp(2 * tj * T[j]) for j in 1:K]), t)
end

function read_hermite_at_nodes(n)
	# read values of hermite polynomials at nodes computed by mathematica from "data/hermites/H_nnnn.mat"
	try
		return readdlm("data/hermites/H_$(lpad(string(n), 4, '0')).mat", '\t', BigFloat, '\n')
	catch ArgumentError
		error("Hermite polynomial data missing for n = $n")
	end
end

function read_t()
	return readdlm("data/t.mat", '\t', BigFloat, '\n')
end

function read_hermite_at_t()
	return readdlm("data/hermite_t.mat", '\t', BigFloat, '\n')
end

function T_n(n, xi, psi, mu)
	# T_n(x) = \sum_j H(j/n) \psi_j(x) \mu_j
	# xi denotes the index of x
	print("Computing point $(xi)... "); flush(stdout)
	s = sum([hfuncinfty(j / n) * psi[j, xi] * mu[j] for j in 1:n])
	println("done"); flush(stdout)
	return s
end

function run_mathematica(n = 128, points = 10000)
	run(`./generate_hermite_at_t.wls $(n) $(points)`)
end

function export_quadrature(n = 256)
	nodes, weights = gausshermite(n)
	writedlm("data/w.mat", weights, '\n')
	writedlm("data/x.mat", nodes, '\n')
end

# x, Tn = cluster()

# # Plot absolute value of T_n
# pyplot()
# plot(x, Tn, seriestype=:scatter)