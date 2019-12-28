using DelimitedFiles
using Plots

pyplot()

t = collect(Iterators.flatten(readdlm("t.mat", '\t', BigFloat, '\n')))
tn = collect(Iterators.flatten(readdlm("tn.mat", '\t', BigFloat, '\n')))

plot(t, abs.(tn), xlims=(-2, 2))
savefig("plot.png")