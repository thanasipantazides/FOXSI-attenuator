using LinearAlgebra
using Plots
using DataFrames
using CSV

# import design space data:
designpath = joinpath(@__DIR__, "../results/designspace_1e5.csv")
designspace = CSV.read(designpath, DataFrame)
E = designspace.energy
v1 = designspace.v1
v2 = designspace.v2
v3 = designspace.v3
r01 = designspace.r01
r02 = designspace.r02
r03 = designspace.r03
absorbprob = designspace.absorbprob

v = [v1'; v2'; v3']
r0 = [r01'; r02'; r03']
slice = NaN*r0
plotdata = permutedims(cat(cat(r0, v, dims=3), slice, dims=3), (1,3,2))

plotlyjs()
plot(vec(plotdata[1, :, :]), vec(plotdata[2, :, :]), vec(plotdata[3, :, :]), legend=false, color=:violet)