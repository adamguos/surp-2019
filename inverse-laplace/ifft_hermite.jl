#=
Compute inverse discrete Hermite polynomial transforms on the quadrature nodes
Input:
	poly_x = column vector of input function evaluated at Hermite polynomial roots
Output:
	poly_c = column vector of transformation coefficients
=#

using FastGaussQuadrature
using LinearAlgebra

include("hermitefn.jl")
include("fft_hermite.jl")

function ifft_hermite(poly_x)
	poly_x = BigFloat.(poly_x)

	quad_n = size(poly_x)[1]

	nodes = BigFloat.(gausshermite(quad_n)[1])	# nodes is a vector containing the roots of the quad_nth Hermite polynomial/function

	bigU = hermitefn(quad_n, nodes)'

	return inv(bigU) * poly_x
end