using LinearAlgebra
using Plots
using LaTeXStrings
using Measures
using DataFrames
using CSV

include("../src/Attenuator3D.jl")
using .Attenuator3D

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

n = Int(length(E)^(1/4))

θ = zeros(size(E))
ϕ = zeros(size(E))
x = zeros(size(E))
γ = Array{Attenuator3D.Particle, 1}(undef, length(E))
for i = 1:length(E)
    v = [v1[i]; v2[i]; v3[i]]
    r0 = [r01[i]; r02[i]; r03[i]]
    θ[i] = acos((v'*[0;0;-1])/norm(v))
    
    vprojxy = v - (v'*[0;0;1])*[0;0;1]
    if norm(vprojxy) == 0
        ϕ[i] = π
    else
        ϕ[i] = acos((vprojxy'*[1;0;0])/norm(vprojxy))
    end

    # WRITE THIS:
    x[i] = 0
    
    γ[i] = Attenuator3D.Particle(r0, v, E[i])
end

E = reshape(E, (n,n,n,n))
θ = reshape(θ, (n,n,n,n))
ϕ = reshape(ϕ, (n,n,n,n))
x = reshape(x, (n,n,n,n))
γ = reshape(γ, (n,n,n,n))
α = reshape(absorbprob, (n,n,n,n))

θunique = unique(x->round(x, digits=8), θ)
ϕunique = unique(x->round(x, digits=8), ϕ)

θgrid = repeat(θunique, outer=[1,n])
ϕgrid = repeat(ϕunique, outer=[1,n])

gr()

# surface(180*θunique/π, 180*ϕunique/π, 1 .- α[Ei, :, :, xi], color=:viridis, xaxis=(L"\theta,\quad [\mathrm{deg}]"), yaxis=(L"\phi,\quad [\mathrm{deg}]"), zaxis=(L"\tau,\quad [-]")

xi = 1

nplots = 10

# surfs = Array{Plots.Plot{Plots.GRBackend}, 1}(undef, nplots)
surfs = Array{Plots.Plot{Plots.GRBackend}, 1}(undef, nplots)

surface()

for Ei = 1:nplots
    surfs[Ei] = surface(
        180*θunique[:]/π, 
        180*ϕunique[:]/π, 
        1 .- vec(α[Ei, :, :, xi]), 
        color=:viridis,
        marker_z=1 .- vec(α[Ei,:,:,xi]),
        colorbar=:none,
        clims=(0.6,0.8),
        xaxis=(L"\theta,\quad [\mathrm{deg}]"), 
        yaxis=(L"\phi,\quad [\mathrm{deg}]"), 
        zaxis=(L"\tau,\quad [-]"),
        include_mathjax="cdn"
    )
end
# current()
plot(surfs..., layout=length(surfs), size=(2000,1000), margin=2.0mm, clims=(0.6,0.8))