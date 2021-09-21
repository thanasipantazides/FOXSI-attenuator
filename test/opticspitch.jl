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

include("../src/Rotations.jl")
using .Rotations

# include the ray tracing library on all cores:
@everywhere include("../src/Attenuator3D.jl")
using .Attenuator3D


# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength.*1e-6
massattenuation = [attenenergy attenlength]

# import Milo's photons:
milopath = joinpath(@__DIR__, "../data/milo_input.csv")
indata = CSV.read(milopath, DataFrame)
inenergy = indata.energy.*1000
inangle = indata.angle.*π/180

# discard photons with energy above 30 keV or below 0.1 keV (can't interpolate):
keepI = findall(min(attenenergy...) .< inenergy .< max(attenenergy...))
inenergy = inenergy[keepI]
inangle = inangle[keepI]

# pick how many photons to sample:
# samplesize = length(keepI)
# sampleI = sample(1:length(inenergy), samplesize, replace=false)

# inenergy = inenergy[sampleI]
# inangle = inangle[sampleI]



opticalaxis = [0;0;1]

# build photons
nγ = length(inenergy)

sourcez = 5e-3
v0 = [0;0;-1]
ϕs = 2*π*rand(nγ)

γ = Array{Attenuator3D.Particle, 1}(undef, nγ)

for i = 1:nγ
    r0 = [0; 0; sourcez]
    v = Rotations.dcm3(ϕs[i])'*Rotations.dcm1(inangle[i])'*v0
    E = inenergy[i]

    γ[i] = Attenuator3D.Particle(r0, v, E)
end



# build attenuator:
# attenuator pitch angles
θrange = [0, 10*π/180]

nθ = 10


Rtop = 21e-6/2
Rbottom = 10e-6/2
htop = 120e-6
hbottom = 175e-6

sandwichheight = 105e-6
bottomz = 0.0
topz = bottomz + hbottom + sandwichheight

pitch = 60e-6
nholes = 200
attenuatornormal = [0.0; 0.0; 1.0]
attenuatorbottompoint = [0.0; 0.0; bottomz]
attenuatortoppoint = attenuatorbottompoint + attenuatornormal*(topz + htop - bottomz)

attenuatordensity = 2.3296e3 # kg/m^3, silicon

attenuatortopholes = Array{Attenuator3D.Cylinder, 2}(undef, nholes, nholes)
attenuatorbottomholes = Array{Attenuator3D.Cylinder, 2}(undef, nholes, nholes)

for i = 1:nholes
    for j = 1:nholes
        attenuatortopholes[i,j] = Attenuator3D.Cylinder(Rtop, htop, [(i-1)*pitch - (nholes-1)/2*pitch; (j-1)*pitch - (nholes-1)/2*pitch; topz], attenuatornormal)

        attenuatorbottomholes[i,j] = Attenuator3D.Cylinder(Rbottom, hbottom, [(i-1)*pitch - (nholes-1)/2*pitch; (j-1)*pitch - (nholes-1)/2*pitch; bottomz], attenuatornormal)
    end
end

attenuator = Attenuator3D.PixelatedAttenuator(
    cat(attenuatortopholes[:], attenuatorbottomholes[:], dims=1), 
    attenuatordensity,
    attenuatortoppoint,
    attenuatorbottompoint,
    attenuatornormal,
    max(Rtop, Rbottom),
    massattenuation
)



# plot
plotlyjs()
Attenuator3D.plotparticles(γ[1:100:end],false)