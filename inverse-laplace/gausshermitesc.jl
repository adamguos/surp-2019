#=
Compute Gauss-Hermite nodes and weights on scattered data
Input:
	n = quadrature degree
	x = column vector of scattered data
Output:
	weights = W such that \sum_j W_j \psi_k(x_j) \psi_0(x_j) = \delta_{k, 0}, k = 0, ..., n
=#

using LinearAlgebra
using Printf

include("hermitefn.jl")

function gausshermitesc(n, x)
	x_len = size(x, 1)

	# Compute matrix of products of Hermite functions
	# V_{j, k} = \psi_j(x_k) \psi_0(x_k) (x_k - x_{k+1}), j = 1...n, k = 1...x_len-1
	V = hermitefn(n, x[1:end-1])
	psi0 = V[1, :]	# first row of V, which is \psi_0(x_k)
	for j = 1:n
		V[j, :] = V[j, :] .* psi0
	end

	# Multiply each column of V by (x_k - x_{k+1})
	for k in 1:x_len-1
		V[:, k] = V[:, k] * (x[k] - x[k+1])
		# V[:, k] = V[:, k] / exp(-x[k]^2 / 2)
	end

	# Initialise column vector of 0s except for the first element, which is 1
	# e_k = 1 if k=0, 0 otherwise, k = 1...n
	e = zeros(n)
	e[1] = 1

	# Compute W using Julia backslash notation (QR)
	W = V \ e
	println("matrix solve done")

	return W
end

# function to check accuracy of Gauss-Hermite weights
function gramcheck(w, x, n)
	quad_n = n ÷ 2
	num_points = size(w, 1)
	g = fill(BigFloat(0), (quad_n, quad_n))	# force g to be a BigFloat type matrix

	psi = hermitefn(quad_n, x)
	for m = 1:quad_n
		for l = 1:quad_n
			for j = 1:num_points
				g[m, l] = g[m, l] + (w[j] * psi[m, j] * psi[l, j] * (x[j] - x[j+1]))
			end
		end
	end

	println("gramcheck done")
	return norm(g - I(quad_n))
end

function test()
	# good numbers: n=128, precision=600
	# n=256, precision=1400
	# random needs higher n than linspaced

	precision = 1400
	setprecision(precision)

	n = 256
	x_len = n * 2
	# x_len = n
	interval_max = sqrt(2) * n

	# x = (rand(x_len, 1) * 2 .- 1) * interval_max	# randomly distributed points between -interval_max and interval_max
	x = collect(range(-interval_max, stop=interval_max, length=x_len))	# equispaced points between -interval_max and interval_max

	println("n = ", n, ", x_len = ", x_len, ", precision = ", precision)

	@time w = gausshermitesc(n, x)	# use Julia built-in time recording
	@time gram = gramcheck(w, x, n)
	@printf "%.5E" gram
end
