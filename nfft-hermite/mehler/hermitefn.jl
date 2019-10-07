#=
Input:
	n = positive integer
	x = 1-D vector
Output:
	y = n * length(x) matrix whose (j, k)-th entry is h_{j-1}(x_k)\exp(-x_k^2/2)
=#

function hermitefn(n, x)
	s = size(x, 1);
	y = zeros(n, s);

	weight = exp.(-x .^ 2 / 2)

	y[1, :] = ones(1, s) .* weight' / pi^(1/4)
	y[2, :] = sqrt(2) .* x .* y[1, :]

	for k = 3:n
		y[k, :] = sqrt(2 / (k-1)) .* x .* y[k-1, :] - sqrt((k-2) / (k-1)) .* y[k-2, :]
		print("hermitefn: $k")
	end

	return y
end