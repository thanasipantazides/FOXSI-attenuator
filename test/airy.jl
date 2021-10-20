using LinearAlgebra
using Plots
using LaTeXStrings
using Measures
# using DataFrames
# using CSV

include("../src/Attenuator3D.jl")
using .Attenuator3D


# Planck
h = 6.62607e-34
# radius of attenuator hole
R = 10e-2
# energies
E = 5e-3*(1.602180000067e-19)

λ = E/h
airycoeff = 2*π*R/λ
display(airycoeff)

# angles to check
θs = LinRange(-10*π/180, 10*π/180, 1000)
I = [Attenuator3D.airyintensity(R, θ, λ) for θ in θs]

plot(θs, I)