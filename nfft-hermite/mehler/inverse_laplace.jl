using FastGaussQuadrature
using LinearAlgebra
using Plots

include("hermitefn.jl")
include("hfuncinfty.jl")
include("observations.jl")
include("gausshermitesc.jl")

function expmult(a, b)
	r = sign(a)
	s = sign(b)
	return exp(log(abs(a)) + log(abs(b))) * r * s
end

# Set precision of BigFloat to 300 bits (~90 digits)
setprecision(3000)

# Set parameters/constants
npoints = 2048
samples = BigFloat.(collect(range(-3, stop=3, length=npoints)))
quad_n = 128
poly_n = quad_n ÷ 2
w = BigFloat((sqrt(5) - 1) / 2)
# a = [1, 3, 2, 1]
a = BigFloat.([1, 3, 1])
# a = [1, 1, 1, 1, 1]
# t = [-1, -0.25, 0.5, 1]
t = BigFloat.([-1, 0, 0.75])
# t = [-1, -0.5, 0, 0.5, 1]

# Compute normalised Gauss-Hermite nodes and weights
# nodes, weights = FastGaussQuadrature.unweightedgausshermite(quad_n)
nodes = BigFloat.(collect(range(-sqrt(2) * quad_n, stop=sqrt(2) * quad_n, length=quad_n*2)))
weights = gausshermitesc(quad_n, nodes)

obs = observations(a, t, nodes)
hermites_nodes = hermitefn(poly_n, nodes)
println(size(hermites_nodes))
println(size(weights))

moments = fill(BigFloat(0), poly_n)
for k = 1:poly_n
	arr = fill(BigFloat(0), quad_n)
	for l = 1:quad_n
		# arr[l] = weights[l] * obs[l] * hermites_nodes[k, l]
		arr[l] = expmult(expmult(weights[l], obs[l]), hermites_nodes[k, l])
	end

	# moments[k] = w^(-k) * sum(arr)
	moments[k] = expmult(w^(-k), sum(arr))
	println("moments: $k")
end

hermites_samples = hermitefn(poly_n, samples)

# kernelsum = zeros(npoints)
kernelsum = fill(BigFloat(0), npoints)
for i = 1:npoints
	arr = fill(BigFloat(0), poly_n)
	for k = 1:poly_n
		# arr[k] = hfuncinfty([k / poly_n])[1] * moments[k] * hermites_samples[k, i]
		arr[k] = expmult(expmult(hfuncinfty([k / poly_n])[1], moments[k]), hermites_samples[k, i])
	end
	kernelsum[i] = abs(sum(arr))
	println("kernelsum: $i")
end

plot!(samples, kernelsum)
plot!(t, seriestype="vline")
savefig("kernelsum.svg")
plot()