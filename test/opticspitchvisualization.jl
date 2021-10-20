using LinearAlgebra
using Plots
using Printf
using StatsBase
using CSV
using DataFrames
using Measures
using LaTeXStrings

include("../src/Rotations.jl")
using .Rotations

resultfiles = readdir(joinpath(@__DIR__,"../results/"))
ncases = 0
θs = []
energy = []
angle = []
transmitprob = []
for file in resultfiles
    if endswith(file, ".serial")
        ncases = ncases + 1
        
        # extract pitch angle (θ) from the file name 
        starti = findnext("0", file, 1)[1]
        endi = findlast("_", file)[1] - 1
        append!(θs, parse(Float64, file[starti:endi]))
        
        datamat = deserialize(joinpath(@__DIR__, "../results/"*file))

        # store file data
        if ncases == 1
            # store energy, angle (should be same for all runs)
            energy = datamat[:,1]
            angle = datamat[:,2]
            
            # start transmitprob
            transmitprob = datamat[:,3]
        else
            # horzcat to transmitprob
            transmitprob = cat(transmitprob, datamat[:,3], dims=2)
        end
    end
end

# sort data
sortθ = sortperm(θs)
θs = θs[sortθ]
transmitprob = transmitprob[:,sortθ]

sortE = sortperm(energy)
energy = energy[sortE]
angle = angle[sortE]

for i = 1:ncases
    transmitprob[:,i] = transmitprob[sortE,i]
end



# make energy bins
binwidth = 0.5e3
bins = 0:binwidth:29e3

# sample probability of falling in each bin
binleft = []
for i = 1:length(bins)
    rh = findfirst(energy .> bins[i])
    push!(binleft,rh)
end

# make counts per bin
baseprob = ones(length(bins), 1)
baseweights = zeros(length(bins), 1)
binnedprob = zeros(length(bins), ncases)
nonzfrac = zeros(length(bins), ncases)

for i = 1:ncases
    for j = 1:length(binleft)
        if j < length(binleft)
            thisprobset = transmitprob[binleft[j]:(binleft[j+1] - 1), i]
            bincount = binleft[j+1] - binleft[j]

            baseweights[j] = bincount
            binnedprob[j,i] = sum(thisprobset)/bincount
            nonzfrac[j,i] = length(findall(thisprobset .!= 0))/bincount
            
        else
            thisprobset = transmitprob[binleft[j]:length(energy), i]
            bincount = length(energy) - binleft[j] + 1
            
            baseweights[j] = bincount
            binnedprob[j,i] = sum(thisprobset)/bincount
            nonzfrac[j,i] = length(findall(thisprobset .!= 0))/bincount

        end
    end
end

# total probability of absorption for each attenuator pitch:
totalprobability = sum(transmitprob, dims=1)/length(energy)
print("\ntotal absorption probability:\n")
display(totalprobability)
print("\n")



# plot
gr()
# plotlyjs()

binskev = bins/1000
# binlabels = [string(Int(bin)) for bin in binskev]
binlabels = [string(round(bin,digits=1)) for bin in binskev]

h = bar(
            binskev,
            baseweights.*baseprob,
            nbins=length(binskev),
            orientation=:vertical,
            linecolor=nothing,
            xaxis=("Energy [keV]"),
            yaxis=("Counts"),
            label="Spectrum from optics"
)

anglerange = [1,3,10]
for i = anglerange
    nonemptys = findall(binnedprob[:,i] .!= 0)
    fullbins = binskev[nonemptys]
    weights = baseweights[nonemptys].*binnedprob[nonemptys,i]

    h = bar!(
        fullbins,
        weights,
        nbins=length(weights),
        bar_edges=true,
        orientation=:vertical,
        linecolor=nothing,
        xticks=([binskev .- 0.5;],binlabels),
        xaxis=("Energy [keV]"),
        yaxis=("Counts"),
        label="attenuator angle: "*@sprintf("%.1f deg",θs[i]*180/π)
    )
end
current()



p = plot(
    binskev,
    baseweights.*baseprob,
    xaxis=("Energy [keV]"),
    yaxis=("Counts",:log10),
    label=nothing
)

for i = anglerange
    nonemptys = findall(binnedprob[:,i] .!= 0)
    fullbins = binskev[nonemptys]
    weights = baseweights[nonemptys].*binnedprob[nonemptys,i]

    p = plot!(
        fullbins,
        weights,
        xticks=([binskev .- 0.5;],binlabels),
        ylims=[10^0,10^4],
        yticks=(10 .^ LinRange(0,4,5), ["1", "10", "100", "1,000", "10,000"]),
        xaxis=("Energy [keV]"),
        yaxis=("Counts",:log10),
        label=nothing,
        legend=true
        # legend=false#"$(@sprintf("%.2f", θs[i]*180/π))"
    )
end
current()

plot(h,p, layout=(1,2), size=(1280,400), margins=6mm, show=true)
# savefig(joinpath(@__DIR__, "../results/attenuator_angle_optic.pdf"))


plotlyjs()

cscheme = cgrad(:viridis, ncases-1, categorical=true, rev=true)
plot(0,0)
for i = 3:ncases
    nonemptys0 = findall(binnedprob[:,2] .!= 0 )
    thisnonempty = findall(binnedprob[:,i] .!= 0)

    intersectI = intersect(nonemptys0, thisnonempty)

    # fullbins = binskev[intersectI]
    # weights0 = baseweights[intersectI].*binnedprob[intersectI,1] 
    # thisweight = baseweights[intersectI].*binnedprob[intersectI,i]
    fullbins = binskev
    weights0 = baseweights.*binnedprob[:,2]
    thisweight = baseweights.*binnedprob[:,i]

    display(i)
    plot!(
            fullbins,
            (thisweight - weights0)./ weights0,
            xticks=([binskev;],binlabels),
            xlim=[0,15],
            color=cscheme[i-1],
            xaxis=("Energy [keV]"),
            yaxis=(L"\frac{T(\theta) - T_0}{T_0}"),
            label=@sprintf("%.2f deg", θs[i]*180/π),
            legend=:bottomright,
            margins=6mm,
            title="Normalized transmission, after optics"
    )
end
current()

# savefig(joinpath(@__DIR__, "../results/attenuator_angle_optic_normalize.pdf"))