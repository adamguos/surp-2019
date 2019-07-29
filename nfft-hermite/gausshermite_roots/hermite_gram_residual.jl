#=
Goal:	Verify accuracy of FastGaussQuadrature.gausshermite using Gram matrix residuals
		G_{m, l} = \sum_j W_j \psi_m(y_j) \psi_l(y_j), m, l = 0, ..., n/2
Input:	W (1, n): weights from output of gausshermite(n)
		H (m, n): H_{j, k} = \psi_j(X_k)
Output: r: Gram residual
=#

using FastGaussQuadrature
using DelimitedFiles
using LinearAlgebra

function hermite_gram_residual(W, H)
	G = zeros(Float64, size(H, 1) ÷ 2, size(H, 1) ÷ 2)

	for m in 1:(size(G, 1))
		for l in 1:(size(G, 2))
			for j in 1:size(W, 1)
				G[m, l] += (W[j] * H[m, j] * H[l, j])
			end
		end
	end

	return norm(G - I, Inf)
end