#=
Compute inverse discrete Hermite polynomial transforms on scattered data
Input:
	poly_x = column vector of input function evaluated at the scattered points
Output:
	poly_c = column vector of transformation coefficients
=#

using FastGaussQuadrature
using LinearAlgebra

function ifftsc_hermite(poly_x)
	
end