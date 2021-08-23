using LinearAlgebra
using StatsBase
using Plots
using CSV
using DataFrames

include("../src/Attenuator3D.jl")
using .Attenuator3D

# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength.*1e-6
massattenuation = [attenenergy attenlength]

# import Milo's post-optics photons:
milopath = joinpath(@__DIR__, "../data/milo_input.csv")
indata = CSV.read(milopath, DataFrame)
inenergy = indata.energy.*1000
inangle = indata.angle.*π/180

# discard photons with energy above 30 keV or below 0.1 keV (can't interpolate):
keepI = findall(min(attenenergy...) .< inenergy .< max(attenenergy...))
inenergy = inenergy[keepI]
inangle = inangle[keepI]

# pick how many photons to sample:
samplesize = 100
sampleI = sample(1:length(inenergy), samplesize, replace=false)

inenergy = inenergy[sampleI]
inangle = inangle[sampleI]

# build focal plane:
focalnormal = [0.0;0.0;1.0]
focaloffset = 4.0



# build attenuator:
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



# probability of hitting each part of attenuator:
Pthru = π*Rbottom^2/pitch^2
Pannulus = π/pitch^2*(Rtop^2 - Rbottom^2)
Psurface = 1 - π*Rtop^2/pitch^2

# average attenuator thickness:
meanthickness = Pthru*sandwichheight + Pannulus*(sandwichheight + htop) + Psurface*(sandwichheight + htop + hbottom)

# build photons:
nphotons = length(inenergy)
photons = Array{Attenuator3D.Particle, 1}(undef, nphotons)

# find max angle to bound distance from attenuator:
maxangle = max(inangle...)
# maximum allowable offset from attenuator to ensure all photons hit attenuator:
dzmax = pitch*nholes/2/tan(maxangle)

safety = 1.5
photonsource = attenuatortoppoint + attenuatornormal*dzmax*(2 - safety)

directions = zeros(3,nphotons)

for i = 1:nphotons
    v0 = [sin(inangle[i]); 0; cos(inangle[i])]
    ϕ = 2*π*rand()
    directions[:,i] = [cos(ϕ) -sin(ϕ) 0;
        sin(ϕ) cos(ϕ) 0;
        0 0 1]*v0
    photons[i] = Attenuator3D.Particle(photonsource, directions[:,i], inenergy[i])
end

# process photons:
absorbprob = Attenuator3D.batchphotons(photons, attenuator)

outenergy = (1 .- absorbprob).*inenergy


sortangle = sortperm(inangle)

plotlyjs()

hin = histogram(inenergy./1000, ylabel="Counts", xlabel="Energy[keV]", legend=false)
hout = histogram(outenergy./1000, xlabel="Energy[keV]", legend=false)

plot(hin, hout)

# xaxis!(:log10)
# yaxis!(:log10)
