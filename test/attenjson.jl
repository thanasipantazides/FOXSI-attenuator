using Distributed

# add a process for each CPU core:
if nprocs() < length(Sys.cpu_info()) + 1
    addprocs(1 + length(Sys.cpu_info()) - nprocs())
end

using LinearAlgebra
using Plots
using StatsBase
using CSV
using DataFrames
using Serialization

# include the ray tracing library on all cores:
@everywhere include("../src/Attenuator3D.jl")
using .Attenuator3D

att = Attenuator3D.importattenuator(joinpath(@__DIR__, "../data/post_optic_attenuator.json"))

plotlyjs()
holes = att.holes[99:102,99:102,:]

ptlist = zeros(3,length(holes))

plot(0,0,0)
for i in 1:length(holes)
    ptlist[:,i] = holes[i].c
    Attenuator3D.plotcyl!(holes[i])
end
current()