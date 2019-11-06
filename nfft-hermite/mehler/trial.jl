setprecision(332)

include("hermitefn.jl")
include("fft_hermite.jl")
include("ifft_hermite.jl")

# test if fft_hermite(ifft_hermite(poly)) returns poly, for poly: y = x^2
poly = BigFloat.(zeros(64, 1))
nodes = BigFloat.(gausshermite(size(poly)[1])[1])
for k in 1:size(poly)[1]
	poly[k] = nodes[k]^2
end

println(norm(poly - fft_hermite(ifft_hermite(poly))))