using Plots
using LinearAlgebra
using Distributed
using SharedArrays
using DataFrames
using CSV


include("../src/Attenuator3D.jl")
using .Attenuator3D



# import Yixian's reference data:
yixianpath = joinpath(@__DIR__, "../data/yixian_baseline.csv")
reference = CSV.read(yixianpath, DataFrame)
refenergy = parse.(Float64, reference.Energy[2:end])
reftransmissivity = reference.Transmission[2:end]

# import attenuation length data from LBNL:
attenpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
attendata = CSV.read(attenpath, DataFrame)
attenenergy = attendata.energy
attenlength = attendata.attenlength
massattenuation = [attenenergy attenlength]

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
sourcez = -0.01
nphotons = 100
photonpitch = pitch*nholes/nphotons
focus = [0; 0; 2e-3]
velocityscale = 1; # scale for numerics
photons = Array{Attenuator3D.Particle, 2}(undef, nphotons, nphotons)

for i = 1:nphotons
    for j = 1:nphotons
        # source of this photon:
        thisposition = [(i-1)*photonpitch - photonpitch/2*nphotons; (j-1)*photonpitch - photonpitch/2*nphotons; sourcez]
        # shoot all photons towards focus:
        thisvelocity = focus - thisposition;
        photons[i,j] = Attenuator3D.Particle(thisposition, velocityscale*thisvelocity./norm(thisvelocity),15e3)
    end
end

@time prob = Attenuator3D.absorptionprobability(photons[1], attenuator)



# plotting
plotlyjs()



# collisions
ax1 = -[1;1;1]
ax1 = ax1./norm(ax1)
ax2 = [1;1;1]
ax2 = ax2./norm(ax2)

cyl1 = Attenuator3D.Cylinder(1, 2, [4;2;4], ax1)
cyl2 = Attenuator3D.Cylinder(1, 2, [0;0;0], ax2)

ray = [-3 6;0 3;0 4]
γ = Attenuator3D.Particle(ray[:,1], (ray[:,2] - ray[:,1])./norm(ray[:,2] - ray[:,1]), 0)

(t1in, t1out) = Attenuator3D.cylinderentryexit(γ, cyl1)
(t2in, t2out) = Attenuator3D.cylinderentryexit(γ, cyl2)

plot(ray[1,:], ray[2,:], ray[3,:], linewidth=3)
Attenuator3D.plotcyl!(cyl1)
Attenuator3D.plotcyl!(cyl2)
xaxis!("x")
yaxis!("y")

current()



# flat spectrum
Erange = [1e3, 30e3]
xrange = [-pitch*(nholes - 1)/2, pitch*(nholes - 1)/2]
yrange = [-pitch*(nholes - 1)/2, pitch*(nholes - 1)/2]
z = -10e-2
dir = [0;0;1]

nγ = 20
γs = Array{Attenuator3D.Particle, 2}(undef, nγ, nγ)
absorblikelihood = zeros(nγ, nγ)
# absorblikelihood = SharedArray{Float64}(nγ*nγ)
xs = xrange[1] .+ (xrange[2] - xrange[1]).*rand(nγ, nγ)
ys = yrange[1] .+ (yrange[2] - yrange[1]).*rand(nγ, nγ)
# Es = Erange[1] .+ (Erange[2] - Erange[1]).*rand(nγ, nγ)
Es = Erange[1]:(Erange[2] - Erange[1])./(nγ^2 - 1):Erange[2]
Es = Matrix(reshape(Es, (nγ,nγ)))

@distributed for i = 1:1
end

# write function in attenuator3D.jl that takes batch of photons and processes in parallel. Call it from here.

for i = 1:nγ*nγ
    γs[i] = Attenuator3D.Particle([xs[i];ys[i];z], dir, Es[i])
    # absorblikelihood[i] = Attenuator3D.absorptionprobability(γs[i], attenuator)
    # println("solved ", 100*(i)/nγ/nγ, "%")
end

@time absorblikelihood = Attenuator3D.batchphotons(γs, attenuator)


# l = @layout [grid(2,1)]
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

bins = 20
inhist = histogram(
    vec(Es)./1000,
    xlabel="Energy [keV]",
    ylabel="Counts",
    legend=false
)

sortI = sortperm(vec(Es))
outhist = plot(vec(Es[sortI])./1000,
    vec(absorblikelihood[sortI]).*100,
    linewidth=2,
    xlabel="Energy [keV]",
    ylabel="Probability of absorption [%]",
    legend=false
)

plot(inhist, outhist)
current()

# compare to Yixian's spectrum:
plot(vec(Es[sortI])./1000,
    (1 .- vec(absorblikelihood[sortI])).*100,
    color=:blue,
    linewidth=2,
    legend=false
)

# must be odd:
nsmooth = 51
smoothprob = moving_average(absorblikelihood[sortI], nsmooth)
dn = Int64((nsmooth - 1)/2)

I1 = dn
I2 = length(absorblikelihood) - dn - 1

plot!(vec(Es[sortI])[I1:I2]./1000,
    (1 .- vec(smoothprob)).*100,
    linewidth=2,
    color=:green,
    legend=false
)

plot!(refenergy, reftransmissivity.*100,
    color=:red,
    linewidth=2,
    xlabel="Energy [keV]",
    ylabel="Probability of transmission [%]",
    legend=false
)

# savefig("../results/uniformspectrum.pdf")
