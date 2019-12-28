#=
Compute discrete Hermite polynomial transforms on the quadrature nodes
Input:
	poly_c = column vector of transformation coefficients
Output:
	poly_x = column vector of function evaluated at Hermite polynomial roots
=#

using FastGaussQuadrature
using LinearAlgebra

include("hermitefn.jl")

function fft_hermite(poly_c)
	quad_n = size(poly_c)[1]

	nodes, weights = gausshermite(quad_n)	# nodes is a vector containing the roots of the quad_nth Hermite polynomial/function

	bigU = zeros(quad_n, quad_n)
	for k = 1:quad_n
		bigU[:, k] = hermitefn(quad_n, nodes[k])
	end

	bigD = Diagonal(nodes)

	return bigU' * poly_c
end