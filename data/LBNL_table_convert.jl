using LinearAlgebra
using DataFrames
using CSV

files = readdir()

rawtable = CSV.read("data/LBNL_raw.csv", DataFrame, header=2)
energy = rawtable.energy
transmission = rawtable.transmission

sidensity = 2.3296e3
thickness = 0.2e-6      # from the table

# attenuation length is multiplying coefficient of length an x-ray photon spends in side

attenlength = -log.(transmission)/thickness

outtable = DataFrame(energy=energy, attenlength=attenlength)

CSV.write("data/LBNL_attenlength_Si_test.csv", outtable)