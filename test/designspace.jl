using Distributed

# add a process for each CPU core:
if nprocs() < length(Sys.cpu_info()) + 1
    addprocs(1 + length(Sys.cpu_info()) - nprocs())
end

using LinearAlgebra
using Plots
using CSV
using DataFrames

include("../src/Rotations.jl")
using .Rotations

# include the library on all cores:
@everywhere include("../src/Attenuator3D.jl")
using .Attenuator3D


# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength.*1e-6
massattenuation = [attenenergy attenlength]


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



# build photons:

# initial x-position
x0 = attenuatorbottomholes[Int(round(nholes/2)), Int(round(nholes/2))].c[1] + pitch/2

# photons: span energy, impact position, angle, optical distance
energyrange = [0.1e3, 30e3]
pitchrange = [0, 2*π/180]
twistrange = [0, π/2]
positionrange = [0, 60e-6] .+ x0

# number of steps in each range:
n = 5

Es = LinRange(energyrange[1], energyrange[2], n)
θs = LinRange(pitchrange[1], pitchrange[2], n)
ϕs = LinRange(twistrange[1], twistrange[2], n)
xs = LinRange(positionrange[1], positionrange[2], n) .+ x0

# basic velocity
v0 = [0;0;-1]
# photon source z-plane
z0 = 0.0005

γ = Array{Attenuator3D.Particle, 4}(undef, n, n, n, n)
for iE = 1:n
    for iθ = 1:n
        for iϕ = 1:n
            for ix = 1:n
                v = Rotations.dcm3(ϕs[iϕ])'*Rotations.dcm2(θs[iθ])'*v0
                rp = [xs[ix]; 0; z0]
                dt0 = (z0 - attenuator.normal'*rp)/(attenuator.normal'*v)
                r0 = rp + dt0*v
                E = Es[iE]

                γ[iE, iθ, iϕ, ix] = Attenuator3D.Particle(r0, v, E)
            end
        end
    end
end

# process photons:
@time absorbprob = Attenuator3D.batchphotons(γ, attenuator)

plotlyjs()
# plot([0;0],[0;0],[0;0])

# for i = 1:4:length(γ)
#     pts = [γ[i].r0 γ[i].v]'
#     plot!(pts[:,1], pts[:,2], pts[:,3], color=:blue, legend=false)
# end

# current()