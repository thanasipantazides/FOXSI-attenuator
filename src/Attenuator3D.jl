module Attenuator3D
using Distributed, DistributedArrays
using SharedArrays
using LinearAlgebra
using Random
import JSON, CSV, DataFrames
using Plots
# using ArbNumerics

import SpecialFunctions

export Particle
struct Particle
    r0  # initial position
    v   # propagation direction (unit)
    E   # energy
end

export Cylinder
struct Cylinder
    R   # radius
    h   # height
    c   # center point
    a   # axis (unit) (NOTE: c + h*a yields a point on the opposite cap => this is an inward normal axis)
end

export PixelatedAttenuator
struct PixelatedAttenuator
    holes
    density     # in kg/m^3
    toppoint    # a point (x,y,z) on the one side of the attenuator (in m)
    bottompoint # a point on the other side of the attenuator
    normal      # normal vector to attenuator surface
    largeradius # largest hole radius in attenuator 
    massattenuation # table of mass attenuation coefficients per energy
end

export so3
"""
    so3(u)

Get skew-symmetric matrix for crossing `u` with another 3-vector. I.e.: so3(u)*v == cross(u,v)
"""
function so3(u)
    return [0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0]
end

export sliceindices
"""
    sliceindices(array, n)

Returns n linear index pairs such that `array` is split into n (roughly) equal pieces.  
"""
function sliceindices(array, n)
    N = length(array)
    k = N ÷ n
    if k == 0
        error("slicing too much")
    end
    m = rem(N, n)

    # [array[(i*k+min(i, m)):((i+1)*k+min(i+1, m))] for i = 1:n]
    starts = [1 + i*k + min(i,m) for i = 0:n-1]
    ends = ends = [i*k + min(i,m) for i = 1:n]

    return [starts ends]
end

export shufflesample
"""
    shufflesample(array)

Returns indeces sampling `array` so that the returned array has length which is the largest multiple of `nworkers()` less than or equal to the length of the original array.
"""
function shufflesample(array)
    N = length(array[:])
    n = nworkers()
    maxI = floor(Int, N / n)*n
    shuff = randperm(N)

    return shuff[1:maxI]
end

export interpolateattenuation
"""
    interpolateattenuation(E, massattenuation)

Interpolate in array `massattenuation` by energy E (in eV). First column of `massattenuation` is energy in MeV, second column is mass attenuation coefficient μ/ρ in m^2/kg. Returns mass attenuation coefficient corresponding to input energy E.
"""
function interpolateattenuation(E, massattenuation)
    energies = massattenuation[:,1]
    coefficients = massattenuation[:,2]

    diff = E .- energies

    if all(diff .> 0) || all(diff .< 0)
        error("interpolant outside range")
    end

    upper = findfirst(diff .≤ 0)
    lower = findlast(diff .≥ 0)

    if lower == upper
        return coefficients[lower]
    else
        return coefficients[lower] + (E - energies[lower])/(energies[upper] - energies[lower])*(coefficients[upper] - coefficients[lower])
    end
end

export cylinderentryexit
"""
    cylinderentryexit(particle, cylinder)

Compute times at which a ray passes through a cylinder. The ray has initial position `particle.r0` and velocity `particle.v`. The cylinder has radius `cylinder.R`, height `cylinder.h`, a point one one of its caps `cylinder.c`, and a unit-length axis pointing from `cylinder.c` towards its other cap `cylinder.a`. 

Two scalars are returned: the time at which the ray enters the cylinder and the time at which it leaves the cylinder. If the ray misses the cylinder, zeros are returned for both entry and exit time.
"""
function cylinderentryexit(particle, cylinder)
    r0 = particle.r0
    v = particle.v
    c = cylinder.c
    a = cylinder.a
    R = cylinder.R
    h = cylinder.h

    ttop = (c'*a - r0'*a)/(v'*a)
    tbottom = (c'*a + h - r0'*a)/(v'*a)

    topcapradius = norm(r0 + v*ttop - c)
    bottomcapradius = norm(r0 + v*tbottom - c - a*h)

    qa = v'*v - (v'*a)^2
    qb = 2*((r0'*v - c'*v) - ((v'*a)*(r0'*a - c'*a)))
    qc = norm(r0 - c)^2 - (r0'*a - c'*a)^2 - R^2

    # particle coaxial with cylinder:
    if qa == 0
        # check if it hits caps, if so return caps times
        if topcapradius ≤ R
            if bottomcapradius ≤ R
                # rectify to zero if particle starts within cylinder
                ttop = max(ttop, 0)
                tbottom = max(tbottom, 0)
                return (min(ttop, tbottom), max(ttop, tbottom))
            else
                error("particle supposedly parallel to cylinder axis")
            end
        else
            # outside cap radius, miss
            if bottomcapradius ≤ R
                error("particle supposedly parallel to cylinder axis")
            end

            return (0, 0)
        end
        
    # photon oblique to cylinder axis:
    else
        root1 = (-qb + sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)
        root2 = (-qb - sqrt(Complex(qb^2 - 4*qa*qc)))/(2*qa)

        if imag(root1) == 0
            root1 = real(root1)
            root2 = real(root2)
            # hits infinite cylinder somewhere
            t1 = min(root1, root2)
            t2 = max(root1, root2)

            # list of all intersection times with cap planes and cylinder
            dt = [tbottom, ttop, t1, t2]
            # order intersection times
            I = sortperm(dt)

            # check if first two intersection times are with cylinder (t1, t2)
            if length(intersect(I[1:2], [3, 4])) == 2
                # miss
                return (0, 0)
            end
            # check if last two intersection times are with cylinder (t1, t2)
            if length(intersect(I[3:4], [3, 4])) == 2
                # miss
                return (0, 0)
            end

            # if you've made it this far, cylinder is hit somewhere
            tin = dt[I[2]]
            tout = dt[I[3]]

            tin = max(tin, 0)
            tout = max(tout, 0)
            return (tin, tout)

        # miss altogether (complex roots for intersection with cylinder) 
        else
            return (0, 0)
        end
    end
end

export lengthsinattenuator
"""
    lengthsinattenuator(photon, attenuator)

Get distances a `photon::Particle` travels within an `attenuator::PixelatedAttenuator` which has cylinders cut out of it. 
"""
function lengthsinattenuator(photon, attenuator)
    topcylinderstimes = Array{Float64}(undef, 0, 2)
    # filter cylinders along path
    hits = falses(length(attenuator.holes[:]))
    for c = 1:length(attenuator.holes[:])
        if norm(cross(photon.v,attenuator.holes[c].a)) == 0
            # if cylinder parallel to ray
            hits[c] = norm(attenuator.holes[c].c - photon.r0 - photon.v*((photon.v'*(attenuator.holes[c].c - photon.r0))./(photon.v'*photon.v))) ≤ attenuator.holes[c].R
        else
            # if cylinder oblique to ray
            hits[c] = cross(photon.v, attenuator.holes[c].a)'*(photon.r0 - attenuator.holes[c].c)/norm(cross(photon.v, attenuator.holes[c].a)) ≤ attenuator.holes[c].R
        end
        
        # if this photon flies within this cylinder's radius, check for collision
        if hits[c]
            (tin, tout) = cylinderentryexit(photon, attenuator.holes[c])

            # make sure it actually hits cylinder
            if tout - tin > 0
                topcylinderstimes = cat(topcylinderstimes, [tin tout], dims=1)
            end
        end
    end
    
    # sort cylinders visited in order of visit time (increasing distance from r0 along v)
    sortedcylinderstimes = zeros(size(topcylinderstimes))
    if length(topcylinderstimes[:,1]) > 1
        topcylinderorder = sortperm(topcylinderstimes[:,1])
        sortedcylinderstimes = topcylinderstimes[topcylinderorder,:]
    else
        sortedcylinderstimes = topcylinderstimes
    end
    
    # times at which photon hits attenuator top surface and bottom surface
    tslabbottom = (attenuator.bottompoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)
    tslabtop = (attenuator.toppoint'*attenuator.normal - photon.r0'*attenuator.normal)/(photon.v'*attenuator.normal)

    tslabin = min(tslabbottom, tslabtop)
    tslabout = max(tslabbottom, tslabtop)

    if length(sortedcylinderstimes) > 0
        # times spent within attenuator material (complement of time spent in cylinders)
        timediffs = [sortedcylinderstimes[:,1]; tslabout] - [tslabin; sortedcylinderstimes[:,2]]
    else
        timediffs = tslabout - tslabin
    end

    # path lengths in material
    lengths = norm(photon.v).*timediffs

    return lengths
end

export absorptionprobability
"""
    absorptionprobability(photon, attenuator)

Return the probability of a photon::Particle being absorbed in attenuator::PixelatedAttenuator.
"""
function absorptionprobability(photon, attenuator)
    lengths = lengthsinattenuator(photon, attenuator)
    if length(lengths) != 0
        transmissionlikelihoods = zeros(BigFloat, length(lengths))
        for i = 1:length(lengths)
            # interpolate attenuation coefficients from table data:
            attenuationcoeff = interpolateattenuation(photon.E, attenuator.massattenuation)
            
            # probability of passing through material:
            transmissionlikelihoods[i] = exp(BigFloat(-lengths[i]*attenuationcoeff))
        end

        # complement of total transmission likelihood is absorption likelihood
        absorptionlikelihood = 1 - prod(transmissionlikelihoods)
        
        return absorptionlikelihood
    else
        return 0.0
    end
end

export transmissionprobability
"""
    transmissionprobability(photon, attenuator)

Return the probability of a photon::Particle being transmitted through attenuator::PixelatedAttenuator.
"""
function transmissionprobability(photon, attenuator)
    lengths = lengthsinattenuator(photon, attenuator)
    if length(lengths) != 0
        transmissionlikelihoods = zeros(BigFloat, length(lengths))
        for i = 1:length(lengths)
            # interpolate attenuation coefficients from table data:
            attenuationcoeff = interpolateattenuation(photon.E, attenuator.massattenuation)
            
            # probability of passing through material:
            transmissionlikelihoods[i] = exp(BigFloat(-lengths[i]*attenuationcoeff))
        end

        # get total transmission likelihood
        return prod(transmissionlikelihoods)
        
    else
        return 1.0
    end
end

export batchphotons
"""
    batchphotons(photons, attenuator)

Compute absorption likelihood for a set of `photons` passing through `attenuator`. If multiple processes are available, computation will be parallelized. To execute in parallel, use the `Distributed` module in the calling context, add processes using `addprocs(n)`, and include this source with the `@everywhere` macro before `using` the module:
    
    @everywhere include("Attenuator3D.jl")
    using .Attenuator3D
"""
function batchphotons(photons, attenuator)
    # check if multiple processes available:
    if nprocs() > 1
        # need to factor jobsize by nworkers()
        spares = nworkers() - rem(length(photons), nworkers())
        # populate extra "dummy" photons to make array correct size
        appendix = Array{Attenuator3D.Particle}(undef, spares, 1)
        fill!(appendix, photons[1])
        inphotons = [photons[:]; appendix]
        
        # size of batch to run on each worker
        batchsize = length(inphotons[:]) ÷ nworkers()

        if rem(length(inphotons[:]), nworkers()) != 0
            error("Cannot distribute job to number of workers")
        end

        # shape array so each worker gets one column
        inphotons = reshape(inphotons, batchsize, nworkers())
        # make distributed array for result on all workers
        transmitlikelihood = dzeros(BigFloat, size(inphotons), workers(), [1, nworkers()])

        @sync @distributed for j in 1:nworkers()
            localtransmitlikelihood = localpart(transmitlikelihood)
            locali = DistributedArrays.localindices(transmitlikelihood)

            for i = 1:length(localtransmitlikelihood)
                localtransmitlikelihood[i] = transmissionprobability(inphotons[i,locali[2][1]], attenuator)
            end
        end

        return Array(transmitlikelihood[1:(end - spares)])

    else
        transmitlikelihood = zeros(size(photons))

        @inbounds for i = 1:length(photons)
            transmitlikelihood[i] = transmissionprobability(photons[i], attenuator)
        end
        return transmitlikelihood
    end
end

export airyintensity
"""
    airyintensity(radius, angle, wavenumber)

Get intensity fraction for diffracted light of given `wavenumber` projected at `angle` through hole of `radius` (following Airy's law for large distances).
"""
function airyintensity(radius, angle, wavelength)
    x = 2*π*radius/wavelength*sin(angle)
    I = (2*SpecialFunctions.besselj1(x)/x)^2
    return I
end

export importattenuator
"""
    importattenuator(file)

Returns a `PixelatedAttenuator` object described by the JSON file.
"""
function importattenuator(file)
    data = JSON.parsefile(file)
    holedata = data["holes"]

    # order: increasing from bottom to top
    Rs = convert(Array{Float64,1}, holedata["radii"])
    hs = convert(Array{Float64,1}, holedata["heights"])
    order = convert(Array{Int64,1}, holedata["order"])

    Rs = Rs[order]
    hs = hs[order]

    hsandwich = Float64(data["sandwich"])
    density = Float64(data["density"])

    nx = Int64(data["nx"])
    ny = Int64(data["ny"])
    pitch = Float64(data["pitch"])
    normal = convert(Array{Float64,1}, data["normal"])
    normal = normal./norm(normal)
    pivot = convert(Array{Float64,1}, data["bottomcenter"])
    attenuationfile = data["attenuationdata"]
    attendata = CSV.read(attenuationfile, DataFrames.DataFrame)
    massattenuation = [attendata.energy attendata.attenlength]

    # transform based on normal direction
    cθ = dot(normal, [0,0,1])
    sθ = -sqrt(1 - cθ^2)
    transform = [1  0   0;
                 0  cθ  sθ;
                 0  -sθ cθ]'

    topholes = Array{Cylinder, 2}(undef, nx, ny)
    botholes = Array{Cylinder, 2}(undef, nx, ny)
    for i = 1:nx
        for j = 1:ny
            a = normal

            x = (i - 1)*pitch - (nx - 1)*pitch/2
            y = (j - 1)*pitch - (ny - 1)*pitch/2
            zbot = 0.0
            ztop = hs[1] + hsandwich
            cbot = transform*[x;y;zbot] + pivot
            ctop = transform*[x;y;ztop] + pivot
            
            # make hole list
            botholes[i,j] = Cylinder(Rs[1], hs[1], cbot, a)
            topholes[i,j] = Cylinder(Rs[2], hs[2], ctop, a)
        end
    end

    attenuator = PixelatedAttenuator(
        cat(topholes, botholes, dims=3), 
        density,
        topholes[1].c + topholes[1].a*topholes[1].h,
        botholes[1].c,
        normal,
        max(Rs...),
        massattenuation
    )
    return attenuator
end

# plotting utilities
export plotcyl
"""
    plotcyl(cylinder::Cylinder)

Plot a transparent cylinder described by in Cylinder struct
"""
function plotcyl(cylinder)
    R = cylinder.R
    h = cylinder.h
    c = cylinder.c
    a = cylinder.a
    a = a./norm(a)

    # angular resolution:
    n = 100

    t = range(0, stop=2π, length=n)
    
    X = R.*[cos.(t) cos.(t)]
    Y = R.*[sin.(t) sin.(t)]
    Z = [zeros(n,1) h*ones(n,1)]

    transform = zeros(3,3)
    if a'*[0;0;1] == 1
        # if cylinder axis parallel to z-axis:
        transform = I
    elseif a'*[0;0;1] == -1
        # if cylinder axis antiparallel to z-axis:
        transform = [1 0 0; 0 -1 0; 0 0 -1]
    else
        # cylinder axis oblique to z-axis:
        ax = cross([0;0;1], a)
        ax = ax./norm(ax)
        cosang = a'*[0;0;1]
        # Rodrigues's formula/axis-angle representation:
        transform = I*cosang + (1 - cosang)*(ax*ax') + [0 -ax[3] ax[2]; ax[3] 0 -ax[1]; -ax[2] ax[1] 0]*sqrt(1 - cosang^2)
    end

    # transform all cylinder points:
    for i = 1:n
        posl = [X[i,1]; Y[i,1]; Z[i,1]]
        posr = [X[i,2]; Y[i,2]; Z[i,2]]

        newposl = transform*posl
        newposr = transform*posr

        X[i,:] = [newposl[1] newposr[1]] + [c[1] c[1]]
        Y[i,:] = [newposl[2] newposr[2]] + [c[2] c[2]]
        Z[i,:] = [newposl[3] newposr[3]] + [c[3] c[3]]
    end

    # draw surface:
    s = surface(X,Y,Z,
            color=:green,
            alpha=0.75,
            legend=false)
    return s
end

export plotcyl!
"""
    plotcyl!(cylinder::Cylinder)

Plot a transparent cylinder described by in Cylinder struct
"""
function plotcyl!(cylinder)
    R = cylinder.R
    h = cylinder.h
    c = cylinder.c
    a = cylinder.a
    a = a./norm(a)

    # angular resolution:
    n = 100

    t = range(0, stop=2π, length=n)

    X = R.*[cos.(t) cos.(t)]
    Y = R.*[sin.(t) sin.(t)]
    Z = [zeros(n,1) h*ones(n,1)]

    transform = zeros(3,3)
    if a'*[0;0;1] == 1
        # if cylinder axis parallel to z-axis:
        transform = I
    elseif a'*[0;0;1] == -1
        # if cylinder axis antiparallel to z-axis:
        transform = [1 0 0; 0 -1 0; 0 0 -1]
    else
        # cylinder axis oblique to z-axis:
        ax = cross([0;0;1], a)
        ax = ax./norm(ax)
        cosang = a'*[0;0;1]
        # Rodrigues's formula/axis-angle representation:
        transform = I*cosang + (1 - cosang)*(ax*ax') + so3(ax)*sqrt(1 - cosang^2)
    end

    # transform all cylinder points:
    for i = 1:n
        posl = [X[i,1]; Y[i,1]; Z[i,1]]
        posr = [X[i,2]; Y[i,2]; Z[i,2]]

        newposl = transform*posl
        newposr = transform*posr

        X[i,:] = [newposl[1] newposr[1]] + [c[1] c[1]]
        Y[i,:] = [newposl[2] newposr[2]] + [c[2] c[2]]
        Z[i,:] = [newposl[3] newposr[3]] + [c[3] c[3]]
    end

    # draw surface:
    s = surface!(X,Y,Z,
            color=:green,
            alpha=0.75,
            legend=false)
    return s
end

export plotparticles
"""
    plotparticles(particles::Array{Particle}, vscale, energycode)

Plot particles as rays with option (set by `energycode` bit) to colorcode. Option `vscale` scales length of velocity.
"""
function plotparticles(particles, vscale, energycode)
    
    n = length(particles)

    r0 = zeros(3, n)
    v = zeros(3, n)

    for i = 1:n
        r0[:, i] = particles[i].r0
        v[:, i] = particles[i].r0 .+ vscale*particles[i].v
    end

    data = cat(r0, v, dims=3)
    data = cat(data, fill(NaN,3,n), dims=3)
    x = vec(data[1,:,:]')
    y = vec(data[2,:,:]')
    z = vec(data[3,:,:]')

    p = plot(x, y, z, color=:yellow, legend=false)
    return p
end

export plotparticles!
"""
    plotparticles!(particles::Array{Particle}, energycode)

Plot particles as rays with option (set by `energycode` bit) to colorcode. Option `vscale` scales length of velocity.
"""
function plotparticles!(particles, vscale, energycode)
    
    n = length(particles)

    r0 = zeros(3, n)
    v = zeros(3, n)

    for i = 1:n
        r0[:, i] = particles[i].r0
        v[:, i] = particles[i].r0 .+ vscale*particles[i].v
    end

    data = cat(r0, v, dims=3)
    data = cat(data, fill(NaN,3,n), dims=3)
    x = vec(data[1,:,:]')
    y = vec(data[2,:,:]')
    z = vec(data[3,:,:]')

    p = plot!(x, y, z, color=:yellow, legend=false)
    return p
end

export plotattenuator
"""
    plotattenuator(attenuator::PixelatedAttenuator)

Plot bounding box of attenuator.
"""
function plotattenuator(attenuator)
    cornerstop = [attenuator.holes[1,1,1] attenuator.holes[1,end,1]; 
                  attenuator.holes[end,1,1] attenuator.holes[end,end,1]]
    cornersbottom = [attenuator.holes[1,1,end] attenuator.holes[1,end,end]; 
                    attenuator.holes[end,1,end] attenuator.holes[end,end,end]]
    corners = cat(cornerstop, cornersbottom, dims=3)

    a = corners[1,1,1].c'*attenuator.normal
    b = corners[1,1,2].c'*attenuator.normal

    extremes = zeros(2,2,2,3)

    if a > b
        # front = a
        if corners[1,1,1].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                    extremes[i,j,2,:] = corners[i,j,2].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c
                    extremes[i,j,2,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                end
            end
        end
    else
        # front = b
        if corners[1,1,2].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                    extremes[i,j,2,:] = corners[i,j,1].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c
                    extremes[i,j,2,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                end
            end
        end
    end

    s = surface(0,0,0,legend=false)
    for i = 1:3
        for j = 1:2
            if i == 1
                surf = extremes[j,:,:,:]
            elseif i == 2
                surf = extremes[:,j,:,:]
            elseif i == 3
                surf = extremes[:,:,j,:]
            end 
            surface!(surf[:,:,1], surf[:,:,2], surf[:,:,3], color=:purple, alpha=0.5, legend=false)
        end
    end
    
    current()
    return s
end

export plotattenuator
"""
    plotattenuator!(attenuator::PixelatedAttenuator)

Plot bounding box of attenuator.
"""
function plotattenuator!(attenuator)
    cornerstop = [attenuator.holes[1,1,1] attenuator.holes[1,end,1]; 
                  attenuator.holes[end,1,1] attenuator.holes[end,end,1]]
    cornersbottom = [attenuator.holes[1,1,end] attenuator.holes[1,end,end]; 
                    attenuator.holes[end,1,end] attenuator.holes[end,end,end]]
    corners = cat(cornerstop, cornersbottom, dims=3)

    a = corners[1,1,1].c'*attenuator.normal
    b = corners[1,1,2].c'*attenuator.normal

    extremes = zeros(2,2,2,3)

    if a > b
        # front = a
        if corners[1,1,1].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                    extremes[i,j,2,:] = corners[i,j,2].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,1].c
                    extremes[i,j,2,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                end
            end
        end
    else
        # front = b
        if corners[1,1,2].a'*attenuator.normal > 0
            # axis points along normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c + corners[i,j,2].a*corners[i,j,2].h
                    extremes[i,j,2,:] = corners[i,j,1].c
                end
            end
        else
            # axis points opposite normal
            for i = 1:2
                for j = 1:2
                    extremes[i,j,1,:] = corners[i,j,2].c
                    extremes[i,j,2,:] = corners[i,j,1].c + corners[i,j,1].a*corners[i,j,1].h
                end
            end
        end
    end

    for i = 1:3
        for j = 1:2
            if i == 1
                surf = extremes[j,:,:,:]
            elseif i == 2
                surf = extremes[:,j,:,:]
            elseif i == 3
                surf = extremes[:,:,j,:]
            end 
            surface!(surf[:,:,1], surf[:,:,2], surf[:,:,3], color=:purple, alpha=0.5, legend=false)
        end
    end
    
    current()
    return s
end

end # /module