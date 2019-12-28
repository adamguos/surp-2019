function hfuncinfty(x)
	# y = zeros(size(x))
	y = fill(BigFloat(0), size(x))
	x = abs.(x)

	ind = findall(x .<= 1/2)
	for i in ind
		y[i] = 1
	end

	ind = findall(1/2 .< x .< 1)
	for i in ind
		y[i] = exp(-exp(BigFloat(2 / (1 - 2 * x[i]))) / (1 - x[i]))
	end

	return y
end