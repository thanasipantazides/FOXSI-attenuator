using LinearAlgebra
using Plots
using DataFrames
using CSV

# import design space data:
designpath = joinpath(@__DIR__, "../results/designspace_1e5csv")
designspace = CSV.read(designpath, DataFrame)
E = designspace.energy
v1 = designspace.v1
v2 = designspace.v2
v3 = designspace.v3
r01 = designspace.r01
r02 = designspace.r02
r03 = designspace.r03
absorbprob = designspace.absorbprob

