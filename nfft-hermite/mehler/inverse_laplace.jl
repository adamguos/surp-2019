using FastGaussQuadrature
using LinearAlgebra
using Plots

include("hermitefn.jl")
include("hfuncinfty.jl")
include("observations.jl")

function expmult(a, b)
	s = sign(b)
	return exp(log(a) + log(abs(b))) * s
end

# Set precision of BigFloat to 300 bits (~90 digits)
setprecision(300)

# Set parameters/constants
npoints = 3072
samples = collect(range(-3, stop=3, length=npoints))
quad_n = 128
poly_n = quad_n ÷ 2
w = (sqrt(5) - 1) / 2
# a = [1, 3, 2, 1]
a = [1, 3, 1]
# a = [1, 1, 1, 1, 1]
# t = [-1, -0.25, 0.5, 1]
t = [-1, 0, 0.75]
# t = [-1, -0.5, 0, 0.5, 1]

# Compute normalised Gauss-Hermite nodes and weights
nodes, weights = FastGaussQuadrature.unweightedgausshermite(quad_n)

obs = observations(a, t, nodes)
hermites_nodes = hermitefn(poly_n, nodes)
println(size(hermites_nodes))

moments = zeros(poly_n)
for k = 1:poly_n
	arr = zeros(quad_n)
	for l = 1:quad_n
		arr[l] = weights[l] * obs[l] * hermites_nodes[k, l]
	end

	# moments[k] = w^(-k) * sum(arr)
	moments[k] = expmult(w^(-k), sum(arr))
	print("moments: $k")
end

hermites_samples = hermitefn(poly_n, samples)

kernelsum = zeros(npoints)
for i = 1:npoints
	arr = zeros(poly_n)
	for k = 1:poly_n
		arr[k] = hfuncinfty(k / poly_n) * moments[k] * hermites_samples[k, i]
	end
	kernelsum[i] = abs(sum(arr))
	print("kernelsum: $i")
end

plot!(samples, kernelsum)
plot!(t, seriestype="vline")
savefig("kernelsum.svg")
plot()