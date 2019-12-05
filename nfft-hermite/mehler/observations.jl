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
	w = BigFloat((sqrt(5) - 1) / 2)
	s = size(x, 1)

	y = fill(BigFloat(0), s)

	for k = 1:s
		y[k] = (BigFloat(pi) * w)^(-1/2) * exp(BigFloat(-sqrt(5) / 2) * x[k]^2)
		obs = fill(BigFloat(0), K)
		for j = 1:K
			obs[j] = a[j] * exp(2 * x[k] * t[j])
		end
		y[k] *= sum(obs)
	end

	return y
end