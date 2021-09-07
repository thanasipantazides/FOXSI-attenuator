using LinearAlgebra
using Plots
using LaTeXStrings
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
    ϕ[i] = acos((v'*[1;0;0])/norm(v))
    
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

plotlyjs()
Ei = 6
xi = 3
surface(unique(θ), unique(ϕ), α[Ei, :, :, xi], xaxis=(L"\theta"), yaxis=(L"\phi"), zaxis=(L"\alpha"))
