#=
Compute Hermite functions by recurrence relation
Input:
	n = positive integer
	x = 1-D vector
Output:
	y = n * length(x) matrix whose (j, k)-th entry is h_{j-1}(x_k)\exp(-x_k^2/2)
=#

# compute hermitefn according to Zhuang's ifft_hermite.m
function hermitefn(n, x)
	x = BigFloat.(x)

	s = size(x, 1)
	y = BigFloat.(zeros(s, n))

	weight = BigFloat.(exp.(-x .^ 2 / 2))

	y[:, 1] = ones(1, s) .* weight' / BigFloat(pi^(1/4))
	y[:, 2] = BigFloat(sqrt(2)) .* x .* y[:, 1]

	for k = 3:n
		y[:, k] = BigFloat(1 / sqrt((k-1)/2)) .* (x .* y[:, k-1] - BigFloat(sqrt((k-2) / 2)) .* y[:, k-2])
	end

	return y'
end