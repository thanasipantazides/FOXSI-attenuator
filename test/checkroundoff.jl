using LinearAlgebra
using Plots
using Printf
using StatsBase
using CSV
using DataFrames

# @everywhere include("../src/Attenuator3D.jl")
include("../src/Attenuator3D.jl")
using .Attenuator3D

include("../src/Rotations.jl")
using .Rotations

# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si_test.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength
massattenuation = [attenenergy attenlength]

θrange = [0, 0.9*π/180]

nθ = 10

θs = LinRange(θrange[1], θrange[2], nθ)

Rtop = 21e-6/2
Rbottom = 10e-6/2
htop = 120e-6
hbottom = 175e-6

sandwichheight = 105e-6
bottomz = 0.0
topz = bottomz + hbottom + sandwichheight

pitch = 60e-6
nholes = 200

attenuators = Array{Attenuator3D.PixelatedAttenuator, 1}(undef, nθ)

for k = 1:nθ
    rotation = Rotations.dcm1(θs[k])'
    attenuatornormal = rotation*[0.0; 0.0; 1.0]
    attenuatorbottompoint = [0.0; 0.0; bottomz]
    attenuatortoppoint = attenuatorbottompoint + attenuatornormal*(topz + htop - bottomz)

    attenuatordensity = 2.3296e3 # kg/m^3, silicon

    attenuatortopholes = Array{Attenuator3D.Cylinder, 2}(undef, nholes, nholes)
    attenuatorbottomholes = Array{Attenuator3D.Cylinder, 2}(undef, nholes, nholes)

    for i = 1:nholes
        for j = 1:nholes

            topholecenter = rotation*[(i-1)*pitch - (nholes-1)/2*pitch; (j-1)*pitch - (nholes-1)/2*pitch; topz]
            attenuatortopholes[i,j] = Attenuator3D.Cylinder(Rtop, htop, topholecenter, attenuatornormal)

            bottomholecenter = rotation*[(i-1)*pitch - (nholes-1)/2*pitch; (j-1)*pitch - (nholes-1)/2*pitch; bottomz]
            attenuatorbottomholes[i,j] = Attenuator3D.Cylinder(Rbottom, hbottom, bottomholecenter, attenuatornormal)
        end
    end

    attenuators[k] = Attenuator3D.PixelatedAttenuator(
        cat(attenuatortopholes, attenuatorbottomholes, dims=3), 
        attenuatordensity,
        attenuatortoppoint,
        attenuatorbottompoint,
        attenuatornormal,
        max(Rtop, Rbottom),
        massattenuation
    )
end
attenuator = attenuators[1]

# make photons:
nγ = 10
r0 = [0;0;5e-2]
v0 = [0;0;-1]
vscale = 1
θ = 0*π/180
v = vscale*Rotations.dcm1(θ)'*v0

E = 131.95

γ = Attenuator3D.Particle(r0, v, E)

h1 = Attenuator3D.Cylinder(Rtop, htop, attenuator.toppoint, attenuator.normal)
h2 = Attenuator3D.Cylinder(Rbottom, hbottom, attenuator.bottompoint, attenuator.normal)
(ttop1, ttop2) = Attenuator3D.cylinderentryexit(γ, h1)
(tbot1, tbot2) = Attenuator3D.cylinderentryexit(γ, h2)
lengths = Attenuator3D.lengthsinattenuator(γ,attenuator)

# α = Attenuator3D.batchphotons(γ, attenuator)
α = Attenuator3D.absorptionprobability(γ, attenuator)
τ = Attenuator3D.transmissionprobability(γ, attenuator)

print("\nnonzero results:")
print("\n\ttop cylinder times:\t\t", (ttop2 - ttop1) != 0)
print("\n\tbottom cylinder times:\t\t", (tbot2 - tbot1) != 0)
print("\n\tlength in attenuator:\t\t", any(lengths .!= 0))
print("\n\ttransmission probability:\t", all(τ .!= 0.0))

plot(attenenergy, exp.(-attenlength.*0.2e-6))
xaxis!(:log10)
# df = DataFrame(alpha=α[:])
# tempname = "tempfile.csv"
# CSV.write(joinpath(@__DIR__,"../results/"*tempname),df)
# readdf = CSV.read(joinpath(@__DIR__,"../results/"*tempname),DataFrame)

# readα = readdf.alpha

# print("\n\tCSV read/write OK:\t\t", all(readα .== α[:]))




# plotlyjs()
# Attenuator3D.plotattenuator(attenuator)
# Attenuator3D.plotparticles!(γ,1e-1,0)


# # read save file
# data = CSV.read(joinpath(@__DIR__, "../results/pre_opt_pitch_atten_0.0_rad.csv"), DataFrame)

# energy = data.energy
# x = data.x
# y = data.y
# absorbprob = data.absorbprob

# # zeroI = findall(absorbprob .== 1.0)
# zeroI = 1:100

# nγ = length(zeroI)
# γ = Array{Attenuator3D.Particle, 1}(undef, nγ)

# for i = 1:nγ
#     γ[i] = Attenuator3D.Particle([x[zeroI[i]]; y[zeroI[i]]; 5e-2], [0;0;-1], energy[zeroI[i]])
# end

# print("\n\nRunning batch...")
# αbatch = Attenuator3D.batchphotons(γ, attenuator)
# print("\nRunning loop...")
# αloop = zeros(nγ)
# for i = 1:nγ
#     αloop[i] = Attenuator3D.absorptionprobability(γ[i], attenuator)
# end

# print("\nWriting to DataFrame...")
# dfbatch = DataFrame(energy=energy[zeroI], x=x[zeroI], y=y[zeroI], absorbprob=αbatch)
# dfloop = DataFrame(energy=energy[zeroI], x=x[zeroI], y=y[zeroI], absorbprob=αloop)

# print("\nWriting CSV...")
# CSV.write(joinpath(@__DIR__, "../results/tempbatch.csv"), dfbatch)
# CSV.write(joinpath(@__DIR__, "../results/temploop.csv"), dfloop)

# print("\nReading CSV...")
# batchdata = CSV.read(joinpath(@__DIR__, "../results/tempbatch.csv"), DataFrame)
# loopdata = CSV.read(joinpath(@__DIR__, "../results/temploop.csv"), DataFrame)

# print("\nChecking...")
# print("\n\tNo 1s in batch:\t\t\t", all(αbatch .!= 1.0))
# print("\n\tNo 1s in loop:\t\t\t", all(αloop .!= 1.0))
# print("\n\tLoop matches batch:\t\t", all(αbatch .== αloop))
# print("\n\tBatch matches CSV batch:\t", all(αbatch .== batchdata.absorbprob))
# print("\n\tLoop matches CSV loop:\t\t", all(αloop .== loopdata.absorbprob))
# print("\n\tBatch matches original:\t\t", all(αbatch .== absorbprob[zeroI]))

# sortI = sortperm(energy[zeroI])
# plot(energy[sortI]/1000, absorbprob[sortI])
# plot!(energy[sortI]/1000, αbatch[sortI])