function hfuncinfty(x)
	y = zeros(size(x))
	x = abs.(x)

	ind = findall(x .<= 1/2)
	for i in ind
		y[i] = 1
	end

	ind = findall(1/2 .< x .< 1)
	for i in ind
		y[i] = exp(-exp(2 / (1 - 2 * x[i])) / (1 - x[i]))
	end

	return y
end