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
absorbprob = []
dfs = Array{DataFrame, 1}(undef, 1)
for file in resultfiles
    if startswith(file, "pre_opt_pitch_atten_")
        ncases = ncases + 1
        
        # extract pitch angle (θ) from the file name 
        starti = length("pre_opt_pitch_atten_") + 1
        endi = findnext("_", file, starti)[1] - 1
        append!(θs, parse(Float64, file[starti:endi]))
        
        df = CSV.read(joinpath(@__DIR__, "../results/"*file), DataFrame)

        # store file data
        if ncases == 1
            # store energy, angle (should be same for all runs)
            energy = df.energy
            
            # start absorbprob
            absorbprob = df.absorbprob
        else
            # horzcat to absorbprob
            absorbprob = cat(absorbprob, df.absorbprob, dims=2)
        end
    end
end

# sort data
sortθ = sortperm(θs)
θs = θs[sortθ]
absorbprob = absorbprob[:,sortθ]

sortE = sortperm(energy)
energy = energy[sortE]

for i = 1:ncases
    absorbprob[:,i] = absorbprob[sortE,i]
end

transmitprob = 1 .- absorbprob



# make energy bins
binwidth = 0.1e3
bins = 0.1:binwidth:29e3

# sample probability of falling in each bin
binleft = []
for i = 1:length(bins)
    rh = findfirst(energy .> bins[i])
    push!(binleft,rh)
end

# make counts per bin
baseprob = ones(length(bins), 1)
binnedprob = zeros(length(bins), ncases)
for i = 1:ncases
    for j = 1:length(binleft)
        # WATCH FOR OFF-BY-ONE

        if j < length(binleft)
            binnedprob[j,i] = sum(transmitprob[binleft[j]:(binleft[j+1] - 1), i])/(binleft[j+1] - binleft[j])
        else
            binnedprob[j,i] = sum(transmitprob[binleft[j]:length(energy), i])/(length(energy) - binleft[j] + 1)
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
binlabels = [string(round(bin)) for bin in binskev]

histin = fit(Histogram, energy/1000, nbins=length(bins))

# h = bar(
#             binskev,
#             histin.weights.*baseprob,
#             nbins=length(binskev),
#             orientation=:vertical,
#             linecolor=nothing,
#             xaxis=("Energy [keV]"),
#             yaxis=("Counts"),
#             label="Input spectrum"
# )

# anglerange = [1,3,10]
# # FIX THIS TO MAKE VISIBLE
# for i = anglerange
#     display(i)
#     nonemptys = findall(binnedprob[:,i] .!= 0)
#     fullbins = binskev[nonemptys]
#     weights = histin.weights[nonemptys].*binnedprob[nonemptys,i]

#     display(length(nonemptys)/length(bins))
#     h = bar!(
#         fullbins,
#         weights,
#         nbins=length(weights),
#         bar_edges=true,
#         orientation=:vertical,
#         linecolor=nothing,
#         xticks=([binskev .- 0.5;],binlabels),
#         xaxis=("Energy [keV]"),
#         yaxis=("Counts"),
#         label="attenuator angle: "*@sprintf("%.1f deg",θs[i]*180/π),
#         legend=:bottomright
#     )
# end
# current()



# p = plot(
#     binskev,
#     histin.weights.*baseprob,
#     xaxis=("Energy [keV]"),
#     yaxis=("Counts",:log10),
#     label=nothing
# )

# for i = anglerange
#     nonemptys = findall(binnedprob[:,i] .!= 0)
#     fullbins = binskev[nonemptys]
#     weights = histin.weights[nonemptys].*binnedprob[nonemptys,i]

#     p = plot!(
#         fullbins,
#         weights,
#         xticks=([binskev .- 0.5;],binlabels),
#         ylims=[10^0,10^4],
#         yticks=(10 .^ LinRange(0,4,5), ["1", "10", "100", "1,000", "10,000"]),
#         xaxis=("Energy [keV]"),
#         yaxis=("Counts",:log10),
#         label=nothing
#         # legend=false#"$(@sprintf("%.2f", θs[i]*180/π))"
#     )
# end
# current()

# plot(h,p, layout=(1,2), size=(1280,400), margins=6mm, show=true)
# savefig(joinpath(@__DIR__, "../results/attenuator_angle.pdf"))



# cscheme = cgrad(:roma, ncases-1, categorical=true, rev=true)
# plot(0,0)
# for i = 2:ncases
#     nonemptys0 = findall(binnedprob[:,1] .!= 0 )
#     thisnonempty = findall(binnedprob[:,i] .!= 0)

#     intersectI = intersect(nonemptys0, thisnonempty)

#     fullbins = binskev[intersectI]
#     weights0 = histin.weights[intersectI].*binnedprob[intersectI,1] 
#     thisweight = histin.weights[intersectI].*binnedprob[intersectI,i]

#     display(i)
#     plot!(
#             fullbins,
#             (thisweight - weights0)./ weights0,
#             xticks=([binskev .- 0.5;],binlabels),
#             xlim=[0,15],
#             color=cscheme[i-1],
#             xaxis=("Energy [keV]"),
#             yaxis=(L"\frac{T(\theta) - T_0}{T_0}"),
#             label=@sprintf("%.2f deg", θs[i]*180/π),
#             legend=:bottomright,
#             margins=6mm,
#             title="Normalized transmission, before optics"
#     )
# end
# current()

# savefig(joinpath(@__DIR__, "../results/attenuator_angle_normalize.pdf"))


plotlyjs()
plot(0,0)
for i = 2:ncases
    nonemptys0 = findall(binnedprob[:,1] .!= 0 )
    thisnonempty = findall(binnedprob[:,i] .!= 0)

    intersectI = intersect(nonemptys0, thisnonempty)

    fullbins = binskev[intersectI]

    weights0 = histin.weights[intersectI].*binnedprob[intersectI,1]
    thisweight = histin.weights[intersectI].*binnedprob[intersectI,i]


    # plot!(
    #         fullbins,
    #         maxvalue,
    #         fillrange=minvalue,
    #         xticks=([binskev .- 0.5;],binlabels),
    #         xlim=[0,15],
    #         color=cscheme[i-1],
    #         alpha=0.4,
    #         xaxis=("Energy [keV]"),
    #         yaxis=(L"\frac{T(\theta) - T_0}{T_0}"),
    #         label=nothing,
    #         margins=6mm,
    # )

    plot!(
            fullbins,
            (thisweight .- weights0)./weights0,
            # xticks=([binskev .- 0.5;],binlabels),
            xlim=[0,29],
            color=cscheme[i-1],
            xaxis=("Energy [keV]"),
            yaxis=(L"\frac{T(\theta) - T_0}{T_0}"),
            label=@sprintf("%.2f deg", θs[i]*180/π),
            legend=false,
            margins=6mm,
            title="Normalized transmission, after optics"
    )
end
current()