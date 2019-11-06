#=
Compute Gauss-Hermite nodes and weights on scattered data
Input:
	n = quadrature degree
	x = column vector of scattered data
Output:
	weights = W such that \sum_j W_j \psi_k(x_j) \psi_0(x_j) = \delta_{k, 0}, k = 0, ..., n
=#

using LinearAlgebra

include("hermitefn.jl")

function gausshermitesc(n, x)
	x_len = size(x, 1)

	# V_{j, k} = \psi_j(x_k) \psi_0(x_k), j = 1...n, k = 1...x_len
	V = hermitefn(n, x)
	psi0 = V[1, :]	# first row of V, which is \psi_0(x_j)
	for k = 1:n
		V[k, :] = V[k, :] .* psi0
	end

	# e_k = 1 if k=0, 0 otherwise, k = 1...n
	e = zeros(n)
	e[1] = 1

	# compute W
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
				g[m, l] = g[m, l] + (w[j] * psi[m, j] * psi[l, j])
			end
		end
	end

	println("gramcheck done")
	return norm(g - I(quad_n))
end

# good numbers: n=128, precision=440
# n=256, precision=1200-1500
# random needs higher n than linspaced

n = 256
x_len = 2048
# x = rand(x_len, 1) * 4 .- 2	# randomly distributed points between -2 and 2
x = collect(range(-2, stop=2, length=x_len))	# equispaced points between -2 and 2
@time w = gausshermitesc(n, x)	# use Julia built-in time recording
@time println(gramcheck(w, x, n))