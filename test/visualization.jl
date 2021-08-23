using LinearAlgebra
using Plots
using DataFrames
using CSV

# import Yixian's reference data:
yixianpath = joinpath(@__DIR__, "../data/yixian_baseline.csv")
reference = CSV.read(yixianpath, DataFrame)
refenergy = parse.(Float64, reference.Energy[2:end])
reftransmissivity = reference.Transmission[2:end]

# flat spectrum data:
flatpath = joinpath(@__DIR__, "../results/absorblikelihood_LBNL.csv")
flatdata = CSV.read(flatpath, DataFrame)
energy = flatdata.energy
absorbprob = flatdata.absorbprob
x = flatdata.x
y = flatdata.y

# post-optic data:
opticpath = joinpath(@__DIR__, "../results/post_optic_atten.csv")
opticdata = CSV.read(opticpath, DataFrame)
inenergy = opticdata.inenergy
inangle = opticdata.angle
opticabsprob = opticdata.absorbprob

sortI = sortperm(vec(energy))



moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

# must be odd:
nsmooth = 411
smoothprob = moving_average(absorbprob[sortI], nsmooth)
dn = Int64((nsmooth - 1)/2)

I1 = dn
I2 = length(absorbprob) - dn - 1



plotlyjs()

plot(vec(energy[sortI])./1000,
    (1 .- vec(absorbprob)).*100,
    color=:chartreuse,
    legend=false
)

plot!(vec(energy[sortI])[I1:I2]./1000,
    (1 .- vec(smoothprob)).*100,
    linewidth=3,
    color=:violet,
    legend=false,
    xlabel="Energy [keV]",
    ylabel="Probability of transmission [%]"
)

plot!(refenergy,
    reftransmissivity.*100,
    color=:black,
    linestyle=:dash,
    linewidth=3,
    legend=false
)
yaxis!(:log10)
current()

transmitprob = 1 .- opticabsprob
outsample = rand(length(transmitprob))
transmitmask = outsample .< transmitprob
outenergy = transmitmask.*inenergy
nbins = 30
binticks = LinRange(0,30,nbins+1)

histogram(inenergy./1000, bins=binticks, color=:blue, alpha=0.6)
histogram!(outsample./1000, bins=binticks, color=:red, alpha=0.6, xaxis=("Energy [keV]", binticks), yaxis=("Counts"))

# scatter(x, y,
#     marker_z=energy, 
#     markersize=2,
#     markerstrokewidth=0,
#     legend=false, 
#     aspect_ratio=1.0
# )

# current()

# savefig("/Users/Thanasi/Documents/University/FOXSI/Rays/results/imageprojection.pdf")

# plot!(refenergy, reftransmissivity.*100,
#     color=:red,
#     linewidth=2,
#     xlabel="Energy [keV]",
#     ylabel="Probability of transmission [%]",
#     legend=false
# )
