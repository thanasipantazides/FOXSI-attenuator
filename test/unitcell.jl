using LinearAlgebra
using Distributed
using Plots
using CSV
using DataFrames

# add a process for each CPU core:
if nprocs() < length(Sys.cpu_info()) + 1
    addprocs(1 + length(Sys.cpu_info()) - nprocs())
end

# include the library on all cores:
@everywhere include("../src/Attenuator3D.jl")
using .Attenuator3D



# build attenuator:
Rtop = 21e-6/2
Rbottom = 10e-6/2
htop = 120e-6
hbottom = 175e-6

sandwichheight = 105e-6
bottomz = 0.0
topz = bottomz + hbottom + sandwichheight

pitch = 60e-6

attenuatornormal = [0.0; 0.0; 1.0]
attenuatorbottompoint = [0.0; 0.0; bottomz]
attenuatortoppoint = attenuatorbottompoint + attenuatornormal*(topz + htop - bottomz)

attenuatordensity = 2.3296e3 # kg/m^3, silicon

# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength
massattenuation = [attenenergy attenlength]

tophole = Attenuator3D.Cylinder(Rtop, htop, attenuatortoppoint - [0;0;htop], attenuatornormal)
bottomhole = Attenuator3D.Cylinder(Rbottom, hbottom, attenuatorbottompoint, attenuatornormal)

attenuator = Attenuator3D.PixelatedAttenuator(
    cat(tophole, bottomhole, dims=1), 
    attenuatordensity,
    attenuatortoppoint,
    attenuatorbottompoint,
    attenuatornormal,
    max(Rtop, Rbottom),
    massattenuation
)



# build photons:
γz = 0.0005
γv = [0; 0; -1]
xrange = [-pitch/2, pitch/2]
yrange = [-pitch/2, pitch/2]
Erange = [10e3, 10e3]
nγ = 20
γ = Array{Attenuator3D.Particle, 2}(undef, nγ, nγ)
x = zeros(nγ,nγ)
y = zeros(nγ,nγ)
E = zeros(nγ,nγ)

for i = 1:nγ
    for j = 1:nγ
        # x[i,j] = xrange[1] + (xrange[2] - xrange[1])/(nγ - 1)*(i - 1)
        # y[i,j] = yrange[1] + (yrange[2] - yrange[1])/(nγ - 1)*(j - 1)

        x[i,j] = pitch/2*(i - 1)/(nγ - 1)*cos((j - 1)/(nγ - 1)*2*π)
        y[i,j] = pitch/2*(i - 1)/(nγ - 1)*sin((j - 1)/(nγ - 1)*2*π)

        E[i,j] = Erange[1] + (Erange[2] - Erange[1])*rand()

        γ[i,j] = Attenuator3D.Particle([x[i,j]; y[i,j]; γz], γv, E[i,j]) 
    end
end

@time absorbprob = Attenuator3D.batchphotons(γ, attenuator)
absorbprob = reshape(absorbprob, (nγ, nγ))

# plot:
upz = attenuatortoppoint[3]
downz = attenuatorbottompoint[3]
boundboxpoints = [-pitch/2  pitch/2 pitch/2 -pitch/2 -pitch/2  pitch/2 pitch/2 -pitch/2;
                  -pitch/2 -pitch/2 pitch/2  pitch/2 -pitch/2 -pitch/2 pitch/2  pitch/2;
                  downz    downz    downz    downz    upz      upz     upz      upz]
draworder = [1, 2, 3, 4, 1, 5, 6, 2, 6, 7, 3, 7, 8, 4, 8, 5]

boundx = boundboxpoints[1,draworder]
boundy = boundboxpoints[2,draworder]
boundz = boundboxpoints[3,draworder]

plotlyjs()
plot(boundx, boundy, boundz, color=:blue, legend=false)

scale = upz/max(absorbprob...)
offset = max(scale.*absorbprob...)
surf = surface!(x, y, scale.*absorbprob .- offset, color=:thermal, alpha=1, legend=false)

scatter!(x, y, γz.*ones(nγ,nγ), color=:yellow, markersize=1, legend=false)

Attenuator3D.plotcyl!(tophole)
Attenuator3D.plotcyl!(bottomhole)
