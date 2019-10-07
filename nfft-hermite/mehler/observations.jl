#=
Input:
	a = 1-D vector, magnitude parameter
	t = 1-D vector, location parameter
	x = 1-D vector
Output:
	y = length(x) vector whose j-th entry is F(x_j)
=#

function observations(a, t, x)
	K = (length(a) == length(t)) ? length(a) : throw(ArgumentError("t and a must be of equal length"))
	w = (sqrt(5) - 1) / 2
	s = size(x, 1)

	y = zeros(s)

	for k = 1:s
		y[k] = (pi * w)^(-1/2) * exp(-sqrt(5) / 2 * x[k]^2)
		obs = zeros(K)
		for j = 1:K
			obs[j] = a[j] * exp(2 * x[k] * t[j])
		end
		y[k] *= sum(obs)
	end

	return y
end