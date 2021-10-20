# using Interpolations
using CSV
using DataFrames
using LinearAlgebra
using Plots

include("../src/Attenuator3D.jl")
using .Attenuator3D

# massattenuation = [1.00000e+03  1.570e+02  1.567e+02;
#                     1.50000e+03  5.355e+01  5.331e+01;
#                     1.83890e+03  3.092e+01  3.070e+01;
#                     1.83890e+03  3.192e+02  3.059e+02;
#                     2.00000e+03  2.777e+02  2.669e+02;
#                     3.00000e+03  9.784e+01  9.516e+01;
#                     4.00000e+03  4.529e+01  4.427e+01;
#                     5.00000e+03  2.450e+01  2.400e+01;
#                     6.00000e+03  1.470e+01  1.439e+01;
#                     8.00000e+03  6.468e+00  6.313e+00;
#                     1.00000e+04  3.389e+00  3.289e+00;
#                     1.50000e+04  1.034e+00  9.794e-01;
#                     2.00000e+04  4.464e-01  4.076e-01;
#                     3.00000e+04  1.436e-01  1.164e-01;
#                     4.00000e+04  7.012e-02  4.782e-02;
#                     5.00000e+04  4.385e-02  2.430e-02;
#                     6.00000e+04  3.207e-02  1.434e-02;
#                     8.00000e+04  2.228e-02  6.896e-03;
#                     1.00000e+05  1.835e-02  4.513e-03;
#                     1.50000e+05  1.448e-02  3.086e-03;
#                     2.00000e+05  1.275e-02  2.905e-03;
#                     3.00000e+05  1.082e-02  2.932e-03;
#                     4.00000e+05  9.614e-03  2.968e-03;
#                     5.00000e+05  8.748e-03  2.971e-03;
#                     6.00000e+05  8.077e-03  2.951e-03;
#                     8.00000e+05  7.082e-03  2.875e-03;
#                     1.00000e+06  6.361e-03  2.778e-03;
#                     1.25000e+06  5.688e-03  2.652e-03;
#                     1.50000e+06  5.183e-03  2.535e-03;
#                     2.00000e+06  4.480e-03  2.345e-03;
#                     3.00000e+06  3.678e-03  2.101e-03;
#                     4.00000e+06  3.240e-03  1.963e-03;
#                     5.00000e+06  2.967e-03  1.878e-03;
#                     6.00000e+06  2.788e-03  1.827e-03;
#                     8.00000e+06  2.574e-03  1.773e-03;
#                     1.00000e+07  2.462e-03  1.753e-03;
#                     1.50000e+07  2.352e-03  1.746e-03; 
#                     2.00000e+07  2.338e-03  1.757e-03]



lbnlpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si.csv")
lbnldata = CSV.read(lbnlpath, DataFrame)
lbnlenergies = lbnldata.energy
lbnlattenlength = lbnldata.attenlength
lbnlmassattenuation = [lbnlenergies lbnlattenlength]

nistpath = joinpath(@__DIR__, "../data/NIST_attenlength_Si.csv")
nistdata = CSV.read(nistpath, DataFrame)
nistenergies = nistdata.energy
nistattenlength = nistdata.attenlength
nistmassattenuation = [nistenergies nistattenlength]

testpath = joinpath(@__DIR__, "../data/LBNL_attenlength_Si_test.csv")
testdata = CSV.read(testpath, DataFrame)
testenergies = testdata.energy
testattenlength = testdata.attenlength
testmassattenuation = [testenergies testattenlength]


# density = 2.3296
# density = 1

plotlyjs()


# attenuation data
# energies = massattenuation[:,1]
# coeffs = massattenuation[:,2]

# fineenergies = min(energies...):(max(energies...) - min(energies...))/10000:max(energies...)

# sample energy space logorithmically
expspace = log10(30.01):0.01:log10(3.0e4)
fineenergies = 10.0.^expspace

finecoeffs = zeros(size(fineenergies))

for i = 1:length(finecoeffs)
    finecoeffs[i] = Attenuator3D.interpolateattenuation(fineenergies[i], lbnlmassattenuation)
end



plot(
    1e-3.*fineenergies, finecoeffs,
    grid=true,
    legend=false,
    color=:red,
    linestyle=:dash
)

plot!(
    1e-3*testenergies, 1 ./(1e-6*testattenlength),
    grid=true,
    legend=false,
    color=:purple
)

plot!(
    1e-3.*lbnlenergies, lbnlattenlength,
    grid=true,
    legend=false,
    color=:green
)

plot!(
    1e-3.*nistenergies, nistattenlength,
    grid=true,
    legend=false,
    color=:blue
)

yaxis!("attenuation length [m]", :log10, minorgrid=true)
xaxis!("E [keV]", :log10, minorgrid=true)

current()