#=
Compute inverse discrete Hermite polynomial transforms on the quadrature nodes
Input:
	poly_x = column vector of input function evaluated at Hermite polynomial roots
Output:
	poly_c = column vector of transformation coefficients
=#

setprecision(332)

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

# test if fft_hermite(ifft_hermite(poly)) returns poly, for poly: y = x^2
poly = BigFloat.(zeros(64, 1))
nodes = BigFloat.(gausshermite(size(poly)[1])[1])
for k in 1:size(poly)[1]
	poly[k] = nodes[k]^2
end

println(norm(poly - fft_hermite(ifft_hermite(poly))))